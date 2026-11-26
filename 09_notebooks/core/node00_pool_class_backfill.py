# -*- coding: utf-8 -*-
"""
V3_clean/core/node00_pool_class_backfill.py
===========================================
Node 00 追加标注 —— 给候选池每个分子打一个【与打分无关】的 pool_class 标签。

背景（checkagain33 + 溯源）：
  Node 00 的并集里混入了约 76 个远小于候选药的碎片（苯、丁烷、甲醇、单氨基酸等），
  以及个别价态不合理的物种（O[Cl+2](O)O、N[SH](=O)=O）。溯源显示它们不是某个
  坏文件的污染，而是【生成器早期未充分训练批次】的小产物，散布在多个 VAE run 中
  （primary_source 集中在 07VAE / 02VAE 等早期运行）。

  这些分子若与真实候选一同对接，会在 BBS 分布里堆出一个虚假的低亲和力簇，而 MD
  选择阈值恰恰从该分布的聚类读取 —— 于是碎片会悄悄改变哪些真分子进入 MD。

原则（与反循环原则一致）：
  - 【不删除】任何分子。全部对接。
  - 标签【只依化学】：分子尺寸 + 价态合理性。【绝不】使用任何 docking / BBS / reward 列。
  - 阈值写进 Methods，可复现、可辩护。
  - 下游聚类只用 pool_class == 'candidate'；fragment / implausible_valence 在 SI 单列报告。

pool_class 取值：
  'implausible_valence' —— RDKit 严格 sanitize 失败（价态/化学图不成立）。最高优先。
  'fragment'            —— 能 sanitize，但重原子数 < HEAVY_ATOM_MIN。
  'candidate'           —— 其余（真实候选药分子）。

关于 implausible_valence（诚实说明）：
  本类别用 RDKit 严格 sanitize 作判据。RDKit 在形式电荷平衡后对价态很宽容 ——
  形式上"怪"但合法的写法（如 O[Cl+2](O)O 含氧酸、[SH](=O)=O 磺酸互变体）都能
  通过 sanitize，因而按【尺寸】落入 fragment/candidate，而非 implausible_valence。
  故本类别在当前池上【预期计数为 0】—— 这是正确的，不是 bug。
  我们【刻意不】加启发式电荷规则去强行捕获这些物种：那会重新引入假阳性风险
  （把 pool_mol_154/577 这类合法大磺酸误标），只为把两个已被判为 fragment 的
  4-重原子微小分子改个名。两种标法下聚类结果完全一致。
  若需在 SI 点名这些物种，那是叙事选择，而非要构建的检测器。

阈值 HEAVY_ATOM_MIN = 12：
  依据重原子数分布的自然断层 —— 池中 [9,12) 仅 15 个（低谷），[12,15) 起回升，
  全池中位数 24。12 是分布本身给出的分界，非任意选取。

输入：Node 00 的 all_candidates_v3_start.csv（含 canonical_smiles + name + source 列）
输出：同表 + 三列 pool_class / heavy_atom_count / pool_class_reason
      —— 回填，不重建；下游节点继承该标签。
"""

import json
import hashlib
from datetime import datetime, timezone
from pathlib import Path

BACKFILL_VERSION = "0.1.0"

HEAVY_ATOM_MIN = 12          # < 此值 → fragment（依分布断层，见 docstring）
SMILES_COL = "canonical_smiles"
NAME_COL = "name"


def utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def valence_ok(smiles, Chem):
    """RDKit 严格核验：能否在【完整 sanitize】下解析。
    这是价态合理性的权威判据 —— 启发式 SMARTS 会误判（如磺酸的 [SH](=O)=O
    互变写法在大分子里是合法的），故必须用 RDKit 而非字符串规则。
    返回 (ok: bool, reason: str)。"""
    try:
        m = Chem.MolFromSmiles(str(smiles))          # 默认 sanitize=True
        if m is None:
            # 关闭 sanitize 再解析，区分"根本非法"与"价态问题"
            m2 = Chem.MolFromSmiles(str(smiles), sanitize=False)
            if m2 is None:
                return False, "unparseable"
            try:
                Chem.SanitizeMol(m2)
                return True, "ok"                     # 理论上不会到这
            except Exception as e:
                return False, f"sanitize_failed:{type(e).__name__}"
        return True, "ok"
    except Exception as e:
        return False, f"error:{type(e).__name__}"


def heavy_atom_count(smiles, Chem):
    """重原子数（RDKit 权威计数；解析失败返回 -1）。"""
    m = Chem.MolFromSmiles(str(smiles), sanitize=False)
    if m is None:
        return -1
    return sum(1 for a in m.GetAtoms() if a.GetSymbol() != "H")


def classify_one(smiles, Chem, heavy_min=HEAVY_ATOM_MIN):
    """返回 (pool_class, heavy_atom_count, reason)。仅依化学，不看任何分数。"""
    ok, vreason = valence_ok(smiles, Chem)
    if not ok:
        ha = heavy_atom_count(smiles, Chem)
        return "implausible_valence", ha, vreason
    ha = heavy_atom_count(smiles, Chem)
    if 0 <= ha < heavy_min:
        return "fragment", ha, f"heavy_atoms={ha} < {heavy_min}"
    return "candidate", ha, f"heavy_atoms={ha}"


def backfill(df, Chem, pandas_mod, smiles_col=SMILES_COL, heavy_min=HEAVY_ATOM_MIN):
    """给 df 追加 pool_class / heavy_atom_count / pool_class_reason 三列。
    返回 (annotated_df, stats)。不改任何既有列、不删任何行。"""
    pd = pandas_mod
    classes, has, reasons = [], [], []
    for smi in df[smiles_col]:
        c, ha, r = classify_one(smi, Chem, heavy_min)
        classes.append(c); has.append(ha); reasons.append(r)
    out = df.copy()
    out["pool_class"] = classes
    out["heavy_atom_count"] = has
    out["pool_class_reason"] = reasons

    vc = out["pool_class"].value_counts().to_dict()
    stats = {
        "n_total": len(out),
        "heavy_atom_min_threshold": heavy_min,
        "pool_class_counts": vc,
        "n_candidate": int((out["pool_class"] == "candidate").sum()),
        "n_fragment": int((out["pool_class"] == "fragment").sum()),
        "n_implausible_valence": int((out["pool_class"] == "implausible_valence").sum()),
        "implausible_examples": out.loc[out["pool_class"] == "implausible_valence",
                                        [NAME_COL, smiles_col]].values.tolist()[:20],
        "criterion": "chemistry only (heavy-atom count + RDKit strict sanitize); "
                     "NO docking/BBS/reward columns used",
    }
    return out, stats


def write_outputs(annotated_df, stats, out_csv, manifest_json, input_path,
                  dry_run=True):
    """回填写盘 + manifest。默认 dry_run（fail-safe）。"""
    out_csv = Path(out_csv)
    manifest = {
        "module": "node00_pool_class_backfill",
        "version": BACKFILL_VERSION,
        "created_at_utc": utc_now(),
        "dry_run": dry_run,
        "input": str(input_path),
        "input_sha256": file_sha256(input_path) if Path(input_path).exists() else None,
        "policy": {
            "pool_class": "candidate / fragment / implausible_valence",
            "criterion": "chemistry only — heavy-atom count + RDKit strict sanitize",
            "heavy_atom_min": HEAVY_ATOM_MIN,
            "threshold_rationale": "natural gap in heavy-atom distribution: [9,12) is a "
                                   "trough (15 mols) before [12,15) rises; pool median 24",
            "no_deletion": "no molecule removed; all remain for docking",
            "no_scoring": "docking/BBS/reward columns never consulted",
            "downstream_use": "cluster only pool_class=='candidate'; fragments & "
                              "implausible_valence reported separately in SI",
        },
        "stats": stats,
    }
    if not dry_run:
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        annotated_df.to_csv(out_csv, index=False)
        Path(manifest_json).write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return manifest


# ---- 溯源辅助（活1）：碎片按来源分布，供 Methods/SI 报告 ----
def trace_fragment_sources(annotated_df, source_col="primary_source"):
    """把 fragment / implausible_valence 按来源 run 分布，返回可报告的表。
    仅【读取】既有 source 列，不做任何筛选。"""
    if source_col not in annotated_df.columns:
        return {"error": f"no column {source_col}"}
    flagged = annotated_df[annotated_df["pool_class"] != "candidate"]
    by_src = flagged[source_col].fillna("(none)").value_counts().to_dict()
    return {
        "n_flagged": int(len(flagged)),
        "flagged_by_source": by_src,
        "note": "fragments are spread across early under-trained VAE runs, "
                "not a single contaminated file",
    }
