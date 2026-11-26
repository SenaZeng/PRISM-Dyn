# -*- coding: utf-8 -*-
"""
V3_clean/core/Node01_finalize.py
================================
Node 01 finalization — 把已冻结的指认结果（v0.8.0 master_microstates.csv）
转成下游可用的最终微观态表。

【本模块不做化学判断】——所有电荷/质子化/立体计数的决定已由
Node01_microstate_enumerate.py (v0.8.0) 完成并冻结。finalize 只做机械的
下游整形：排除、立体展开、ID 方案、标记。这是"化学身份 → 空间结构身份"的过渡。

策略（三方商定：user / checkagain28 / Claude）：
  1. 排除 C 类（review_tier == 'C_exclude'；数量随 scope 变化，见 manifest
     的 n_excluded_classC —— 当前 1606 分子 scope 下为 0，不再硬编码 14）
     → node01_excluded.csv
  2. 立体枚举：
       k == 0                → 单行，无后缀，parent_id = 自身
       1 <= k <= 3           → 自动枚举 2^k 个立体异构体，
                               name = <parent>_s01, _s02, ...（补零），
                               parent_id = <parent>，R/S 记于 isomeric_smiles
       k >= 4                → 方案 B：保留【单行】、needs_stereo_review=True、
                               stereo_unresolved_count=k、held out of docking（不枚举）
  3. 极端电荷 |q| >= 3      → extreme_charge_flag=True（不排除，供 pre-docking 审视）
  4. possible_overprotonation → 记录标志，供 pre-docking 单质子化审视
     （注：v0.8.0 已算好该标志；本模块不改电荷，仅传递标志）

ID 方案：
  非枚举:  name=pool_mol_015,      parent_id=pool_mol_015
  枚举:    name=pool_mol_027_s01,  parent_id=pool_mol_027
  后缀【只】给被枚举的分子；parent_id【始终】填充。

输出：
  V3_microstate_table_final.csv  —— Node 02 的冻结输入
  node01_excluded.csv            —— 14 个 C 类排除
  Node01_finalize_manifest.json  —— 统计 + provenance
"""

import json
import hashlib
from datetime import datetime, timezone
from pathlib import Path
from itertools import product

FINALIZE_VERSION = "0.2.2"

# k <= STEREO_AUTO_MAX 自动枚举；k >= 该值+1 走方案 B（保留+标记，不枚举）
STEREO_AUTO_MAX = 3
EXTREME_CHARGE_ABS = 3          # |q| >= 此值 → extreme_charge_flag

# v0.2.0: 枚举安全上限（防组合爆炸）。与 STEREO_AUTO_MAX 无关——它只是硬天花板，
# 真实预算按【每个分子自己的】潜在中心数计算（修复 v0.1.1 的 maxIsomers 恒为 8）。
ENUM_SAFETY_CEILING = 256

# v0.2.0: 门控与操作必须使用【同一套】排除原子集合，故直接复用 enumerate v0.9.0。
# 版本耦合是刻意的：若 enumerate 未升级到 v0.9.0，finalize 拒绝运行而不是静默
# 退回旧口径（那正是 v0.1.1 幻影立体中心问题的成因）。
try:
    from core import Node01_microstate_enumerate as _ENUM
except ImportError:  # 允许与模块同目录直接导入
    try:
        import Node01_microstate_enumerate as _ENUM
    except ImportError as _e:
        raise ImportError(
            "Node01_finalize v0.2.0 requires Node01_microstate_enumerate v0.9.0 "
            "(for excluded_stereo_atoms / count_undefined_stereo_corrected). "
            f"Import failed: {_e}")
if not hasattr(_ENUM, "excluded_stereo_atoms"):
    raise ImportError(
        "Node01_microstate_enumerate is older than v0.9.0 — it lacks "
        "excluded_stereo_atoms(). Upgrade it before running finalize v0.2.0.")


def utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ---------- v0.1.1 hardening helpers (lineage / eligibility / stable labels) ----------
def identity_hash(iso_smiles, net_charge, Chem):
    """Stable per-row identity: sha256 of RDKit-canonical isomeric SMILES + charge.
    Downstream nodes key on THIS, not on names/strings (Risk 1)."""
    canon = str(iso_smiles)
    m = Chem.MolFromSmiles(str(iso_smiles))
    if m is not None:
        canon = Chem.MolToSmiles(m, isomericSmiles=True)
    return hashlib.sha256(f"{canon}|{int(net_charge)}".encode("utf-8")).hexdigest()


def charge_risk_level(q):
    q = abs(int(q))
    return "extreme" if q >= 5 else "high" if q >= 3 else "normal"


def overprotonation_reason(review_reason, possible_overprotonation):
    if not possible_overprotonation:
        return "none"
    rr = str(review_reason)
    # v0.9.1: labels corrected to match true connectivity (were off by one).
    if "1,2-diammonium" in rr:
        return "adjacent_1_2_diamine"        # ethylenediamine / piperazine spacing
    if "1,1-geminal_diammonium" in rr:
        return "geminal_1_1_diamine"
    if "hydrazinium" in rr:
        return "hydrazinium"
    return "adjacent_amine_unspecified"


def exclusion_reason_code(review_reason, net_charge):
    """Stable machine label for excluded (C) molecules (Risk 6)."""
    rr = str(review_reason)
    if "hydrazinium" in rr:
        return "HYDRAZINIUM_NONPHYSIOLOGICAL"
    if abs(int(net_charge)) >= 5:
        return "EXTREME_CHARGE"
    return "CLASS_C_UNSPECIFIED"


def downstream_eligibility(needs_stereo_review, crisk, possible_overprotonation):
    """Node 01 decides eligibility so Node 02 never re-interprets flags (Risk 3).

    v0.9.1 SCOPE FIX — charge_risk no longer holds a molecule back.
    Folding `crisk in ("high","extreme")` into HOLD conflated two questions:
      "what IS this molecule at pH 7.4?"      -> Node 01's job, and it answered
      "is a trication a credible lead?"       -> druglikeness, Node 04+
    A correctly-assigned +3 (e.g. three independent alkylammonium arms) is not an
    identity defect. Such molecules are released as READY and carry
    extreme_charge_flag / charge_risk_level downstream, where the filtering
    decision properly belongs. Only genuine identity doubt (over-protonation)
    still holds.
    """
    if needs_stereo_review:
        return "HOLD_STEREO"
    if possible_overprotonation:
        return "HOLD_CHARGE_REVIEW"
    return "READY"


# ---------- 立体枚举 ----------
def _clear_phantom_stereo(mol, Chem):
    """把【不该被枚举】的原子(可变质子化氮 / 桥连双环桥头)的手性标记清除。
    这样 RDKit 即便枚举了它们，产物也会在去重时折叠回同一个分子。"""
    try:
        ex = _ENUM.excluded_stereo_atoms(mol, Chem)
        for idx in ex:
            mol.GetAtomWithIdx(idx).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    except Exception:
        pass
    return mol


def enumerate_stereoisomers(smiles, k_undefined, Chem, max_k=STEREO_AUTO_MAX,
                            ceiling=ENUM_SAFETY_CEILING):
    """
    枚举【真实】未定义立体中心的异构体，返回 (isomeric SMILES 列表, info)。

    v0.2.0 相对 v0.1.1 的三处修复：
      1) maxIsomers 按【每个分子自己】的潜在中心数计算，不再恒为 2**max_k(=8)。
         v0.1.1 收下了 k_undefined 却从未使用它 —— 静默截断无任何标志。
      2) 产物上清除幻影立体(质子化氮 / 桥头)后再规范化去重，故仅在【真实】中心
         上展开；仅在氮上不同的重复异构体自动折叠。
      3) 截断被显式检测并上报(info["truncated"])，不再冒充"完整枚举"。

    info = {n_raw_generated, n_unique, k_used, k_raw, budget, truncated, note}
    """
    from rdkit.Chem.EnumerateStereoisomers import (
        EnumerateStereoisomers, StereoEnumerationOptions)
    info = {"n_raw_generated": 0, "n_unique": 0, "k_used": int(k_undefined),
            "k_raw": None, "budget": 0, "truncated": False, "note": ""}
    m = Chem.MolFromSmiles(str(smiles))
    if m is None:
        info["note"] = "parse_failed"
        return [], info

    # 预算：按 RDKit 实际可能看到的中心数(k_raw)给足空间，避免"去重前就被截断"
    try:
        k_corr, k_raw, _, _ = _ENUM.count_undefined_stereo_corrected(m, Chem)
    except Exception:
        k_corr, k_raw = int(k_undefined), int(k_undefined)
    info["k_raw"] = k_raw
    exp = max(int(k_raw), int(k_undefined), 1)
    budget = min(2 ** exp, ceiling)
    info["budget"] = budget

    opts = StereoEnumerationOptions(onlyUnassigned=True, maxIsomers=budget,
                                    unique=True)
    out = []
    try:
        for iso in EnumerateStereoisomers(m, options=opts):
            info["n_raw_generated"] += 1
            iso = _clear_phantom_stereo(iso, Chem)          # 去幻影
            Chem.AssignStereochemistry(iso, cleanIt=True, force=True)
            out.append(Chem.MolToSmiles(iso, isomericSmiles=True, canonical=True))
    except Exception as e:
        info["note"] = f"enumeration_error: {type(e).__name__}"
        return [], info

    uniq = sorted(set(out))
    info["n_unique"] = len(uniq)
    # 截断检测：生成数触到预算，且预算小于理论完整数 → 结果不完整
    if info["n_raw_generated"] >= budget and budget < 2 ** exp:
        info["truncated"] = True
        info["note"] = (f"TRUNCATED at safety ceiling {ceiling} "
                        f"(potential centres={exp})")
    # 一致性检查：去幻影后的唯一数不应超过 2**k_corrected
    if k_corr is not None and len(uniq) > 2 ** max(int(k_corr), 0):
        info["note"] = (info["note"] + "; " if info["note"] else "") + \
            f"unique({len(uniq)}) > 2^k_corrected({2 ** max(int(k_corr), 0)})"
    return uniq, info


# ---------- 主 finalize ----------
def finalize(df, Chem, pandas_mod,
             smiles_col="microstate_smiles",
             name_col="name",
             charge_col="net_charge",
             stereo_col="stereo_undefined_centers",
             tier_col="review_tier",
             overprot_col="possible_overprotonation"):
    """消费 v0.8.0 指认表，产出 (final_df, excluded_df, stats)。不改任何化学判断。"""
    pd = pandas_mod

    def _int(x, default=0):
        try:
            return int(float(x))
        except (TypeError, ValueError):
            return default

    # ---- 1. 排除 C 类（尊重 py 的 tier 决定，不重新分类）----
    is_excluded = df[tier_col].astype(str) == "C_exclude"
    excluded_df = df[is_excluded].copy()
    kept = df[~is_excluded].copy()

    final_rows = []
    stats = {
        "n_input": len(df),
        "n_excluded_classC": int(is_excluded.sum()),
        "n_kept_parents": len(kept),
        "n_stereo_auto_parents": 0,
        "n_stereo_isomers_generated": 0,
        "n_stereo_held_k_ge_4": 0,
        "n_no_stereo": 0,
        "n_extreme_charge": 0,
        "n_overprotonation_flagged": 0,
        "stereo_enumeration_failures": [],
        "n_stereo_collapsed_single": 0,
        "n_stereo_truncated": 0,
        "stereo_truncated_parents": [],
    }

    for _, r in kept.iterrows():
        name = str(r[name_col])
        q = _int(r[charge_col])
        k = _int(r[stereo_col])
        micro = str(r[smiles_col])
        overprot = str(r.get(overprot_col, "False")) in ("True", "true", "1")
        extreme = abs(q) >= EXTREME_CHARGE_ABS
        if extreme:
            stats["n_extreme_charge"] += 1
        if overprot:
            stats["n_overprotonation_flagged"] += 1

        base = {
            "parent_id": name,
            "canonical_smiles": r.get("canonical_smiles"),
            "microstate_smiles": micro,   # v0.2.0: 母体未枚举的指认态(回退/溯源)
            "net_charge": q,
            "extreme_charge_flag": extreme,
            "possible_overprotonation": overprot,
            "microstate_status": r.get("protonation_status"),
            "resolved_by": r.get("resolved_by"),
            "_review_reason": str(r.get("review_reason", "")),
            # v0.2.0 血缘：保留母体原始 k(修复 D3),并记录枚举是否被截断
            "stereo_k_original": k,
            "stereo_enum_truncated": False,
            "stereo_enum_note": "",
        }

        if k == 0:
            # 无未定义立体中心 → 单行，无后缀
            stats["n_no_stereo"] += 1
            final_rows.append({
                **base, "name": name,
                "isomeric_smiles": micro,
                "stereo_status": "defined_or_none",
                "needs_stereo_review": False,
                "stereo_unresolved_count": 0,
            })

        elif 1 <= k <= STEREO_AUTO_MAX:
            # 自动枚举 2^k 个立体异构体
            isomers, einfo = enumerate_stereoisomers(micro, k, Chem, STEREO_AUTO_MAX)
            base["stereo_enum_truncated"] = bool(einfo.get("truncated"))
            base["stereo_enum_note"] = str(einfo.get("note", ""))
            if einfo.get("truncated"):
                stats["n_stereo_truncated"] += 1
                stats["stereo_truncated_parents"].append(name)
            if not isomers:
                stats["stereo_enumeration_failures"].append(name)
                # 枚举失败 → 退回方案 B（保留 + 标记），不丢
                final_rows.append({
                    **base, "name": name, "isomeric_smiles": micro,
                    "stereo_status": "enumeration_failed_held",
                    "needs_stereo_review": True, "stereo_unresolved_count": k,
                })
                continue

            if len(isomers) == 1:
                # v0.2.0: 去掉幻影中心后只剩一个异构体 → 本就没有真实立体待定。
                # 不加 _sNN 后缀（避免制造无意义的“立体异构体”身份）。
                stats["n_stereo_collapsed_single"] += 1
                final_rows.append({
                    **base, "name": name, "isomeric_smiles": isomers[0],
                    "stereo_status": "no_real_stereo_after_correction",
                    "needs_stereo_review": False,
                    "stereo_unresolved_count": 0,
                })
                continue

            stats["n_stereo_auto_parents"] += 1
            stats["n_stereo_isomers_generated"] += len(isomers)
            for i, iso in enumerate(isomers, start=1):
                final_rows.append({
                    **base,
                    "name": f"{name}_s{i:02d}",     # 补零
                    "isomeric_smiles": iso,
                    "stereo_status": "auto_enumerated",
                    "needs_stereo_review": False,
                    "stereo_unresolved_count": 0,
                })

        else:  # k >= 4 → 方案 B：保留 + 标记 + 不枚举 + 不进 docking
            stats["n_stereo_held_k_ge_4"] += 1
            final_rows.append({
                **base, "name": name, "isomeric_smiles": micro,
                "stereo_status": "held_high_stereo",
                "needs_stereo_review": True,
                "stereo_unresolved_count": k,
            })

    # ---- v0.1.1: enrich every row with lineage / eligibility / stable labels ----
    for row in final_rows:
        iso = row["isomeric_smiles"]
        q = row["net_charge"]
        nm = row["name"]
        crisk = charge_risk_level(q)
        row["identity_hash"] = identity_hash(iso, q, Chem)
        row["charge_risk_level"] = crisk
        row["overprotonation_reason"] = overprotonation_reason(
            row.get("_review_reason", ""), row["possible_overprotonation"])
        # stereo lineage
        ss = row["stereo_status"]
        if ss == "auto_enumerated":
            row["stereo_generation_method"] = "RDKit_EnumerateStereoisomers"
            row["stereo_parent"] = row["parent_id"]
            row["stereo_index"] = nm.split("_s")[-1] if "_s" in nm else ""
        elif ss == "no_real_stereo_after_correction":
            row["stereo_generation_method"] = "phantom_centres_removed_v090"
            row["stereo_parent"] = row["parent_id"]
            row["stereo_index"] = ""
        elif ss in ("held_high_stereo", "enumeration_failed_held"):
            row["stereo_generation_method"] = "held_not_enumerated"
            row["stereo_parent"] = row["parent_id"]
            row["stereo_index"] = ""
        else:  # defined_or_none
            row["stereo_generation_method"] = "none"
            row["stereo_parent"] = row["parent_id"]
            row["stereo_index"] = ""
        row["downstream_eligibility"] = downstream_eligibility(
            row["needs_stereo_review"], crisk, row["possible_overprotonation"])
        row.pop("_review_reason", None)

    cols = ["name", "parent_id", "identity_hash",
            "canonical_smiles", "microstate_smiles", "isomeric_smiles",
            "net_charge", "charge_risk_level",
            "microstate_status", "resolved_by",
            "stereo_status", "stereo_generation_method", "stereo_parent",
            "stereo_index", "needs_stereo_review", "stereo_unresolved_count",
            "stereo_k_original", "stereo_enum_truncated", "stereo_enum_note",
            "extreme_charge_flag", "possible_overprotonation",
            "overprotonation_reason", "downstream_eligibility"]
    final_df = pd.DataFrame(final_rows)[cols]
    stats["n_final_rows"] = len(final_df)
    stats["downstream_eligibility_counts"] = (
        final_df["downstream_eligibility"].value_counts().to_dict())
    stats["charge_risk_counts"] = final_df["charge_risk_level"].value_counts().to_dict()

    # 排除表：保留原因供审计 + 稳定标签码
    exc_cols = [name_col, charge_col, stereo_col, tier_col, "review_reason",
                smiles_col]
    exc_cols = [c for c in exc_cols if c in excluded_df.columns]
    excluded_df = excluded_df[exc_cols].copy()
    if len(excluded_df):
        excluded_df["exclusion_reason_code"] = [
            exclusion_reason_code(rr, q) for rr, q in
            zip(excluded_df.get("review_reason", ""), excluded_df[charge_col])]

    return final_df, excluded_df, stats


def write_outputs(final_df, excluded_df, stats, out_dir, input_path,
                  dry_run=True):
    """写 V3_microstate_table_final.csv + node01_excluded.csv + manifest。"""
    out_dir = Path(out_dir)
    final_path = out_dir / "V3_microstate_table_final.csv"
    exc_path = out_dir / "node01_excluded.csv"
    mani_path = out_dir / "Node01_finalize_manifest.json"

    manifest = {
        "module": "Node01_finalize",
        "version": FINALIZE_VERSION,
        "created_at_utc": utc_now(),
        "dry_run": dry_run,
        "input": str(input_path),
        "input_sha256": file_sha256(input_path) if Path(input_path).exists() else None,
        "policy": {
            "exclude": "review_tier == 'C_exclude' (respects enumerate tier decision)",
            "stereo_auto_max_k": STEREO_AUTO_MAX,
            "stereo_k_ge_4": "Option B: keep + needs_stereo_review + hold from docking",
            "extreme_charge_abs": EXTREME_CHARGE_ABS,
            "id_scheme": "enumerated -> <parent>_sNN (zero-padded); parent_id always filled",
            "chemistry": "NOT re-decided here; frozen from enumerate v0.9.0",
            "stereo_perception": ("undefined stereocentres counted on the corrected "
                                  "atom set (enumerate v0.9.0): protonation-created "
                                  "ammonium N and bridged-bicyclic bridgeheads are "
                                  "excluded; fused-ring and quaternary-ammonium "
                                  "centres are retained"),
            "enum_budget": "maxIsomers = 2**(per-molecule potential centres), "
                           f"capped at {ENUM_SAFETY_CEILING}; truncation is flagged",
            "requires": "Node01_microstate_enumerate >= v0.9.0",
        },
        "stats": stats,
    }

    if not dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)
        final_df.to_csv(final_path, index=False)
        excluded_df.to_csv(exc_path, index=False)
        mani_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False),
                             encoding="utf-8")
    return manifest, {"final": final_path, "excluded": exc_path, "manifest": mani_path}
