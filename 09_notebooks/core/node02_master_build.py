# -*- coding: utf-8 -*-
"""
V3_clean/core/node02_master_build.py
====================================
Node 02 — 3D 权威分子构造（化学身份 → 空间结构身份）。

消费 Node 01 finalize 的 `V3_microstate_table_final.csv`，为每行构造【一个权威
3D 主分子对象】（master SDF），并强制核心不变量：

    CanonicalIsoSMILES(built master SDF) == CanonicalIsoSMILES(input isomeric_smiles)

即：构造出的 3D 分子可证与输入的化学身份（含立体、电荷）完全一致 —— 3D 化只加
坐标，不改身份。这是"one molecule object"原则在空间层的落实。

【本模块不做化学判断】——电荷/质子化/立体已由 Node 01 冻结。Node 02 只：
  - 用 isomeric_smiles（已含 R/S 与电荷）构造 3D
  - 验证不变量（身份未变）
  - 写 master SDF + identity manifest（每分子）
  - QC 门：不变量失败 = 该分子 QC_FAIL，记录、不静默通过

策略：
  - needs_stereo_review=True（k≥4 保留）的分子【不构造 3D】（held out of docking），
    只登记、标记，等人工立体决定后再回来。
  - 其余全部构造。构造失败（嵌入失败/sanitize/不变量不符）→ QC_FAIL，进 failures。

依赖 core/master_molecule.py（build_master_molecule / identity_manifest /
canonical_iso_smiles / write_master_sdf）。
"""

import json
import hashlib
from datetime import datetime, timezone
from pathlib import Path

NODE02_VERSION = "0.2.0"

# v0.1.3:嵌入迭代上限。RDKit 默认 0 = 自动 = 10 x 原子数(约 500)。
# 原硬编码 50,对 ~50 重原子骨架约少 10 倍尝试,会把可嵌入分子误判为不可行。
EMBED_MAX_ITERATIONS = 0


# v0.2.0: 上游血缘。此前链路仅由输入 CSV 的 sha256 锚定；这里额外锚定
# Node 01 finalize 的 manifest 本身（版本 + 文件 hash），使可追溯性完整。
def _read_node01_provenance(input_dir):
    """读取同目录 Node01_finalize_manifest.json，返回 finalize 版本与 manifest hash。
    缺失时优雅返回 None（不假死、不写死数字）。"""
    p = Path(input_dir) / "Node01_finalize_manifest.json"
    if not p.is_file():
        return None
    try:
        fm = json.loads(p.read_text(encoding="utf-8"))
        return {
            "finalize_version": fm.get("version"),
            "manifest_sha256": file_sha256(p),
            "input_sha256_declared": fm.get("input_sha256"),
        }
    except Exception:
        return None


def _rdkit_version():
    try:
        import rdkit
        return rdkit.__version__
    except Exception:
        return None


def _master_module_version():
    """尽力获取 master_molecule 版本，不改 write_outputs 的调用签名。"""
    try:
        from core import master_molecule as _mm
        return getattr(_mm, "MASTER_MODULE_VERSION", None)
    except Exception:
        try:
            import master_molecule as _mm
            return getattr(_mm, "MASTER_MODULE_VERSION", None)
        except Exception:
            return None


# v0.2.0: 记录"嵌入静默固定"的立体中心 —— 即 Node 01 微观态未指定、但 3D 化
# 必须落定的中心。这不是缺陷，而是一次被记录下来的具体空间选择：
#   - labile_N     : 四配位质子化氮（溶液中快速互变；当前只表征其一）
#   - symmetric_C  : 组成对称的碳（两种"构型"其实是同一分子；仅需登记）
#   - other        : 需人工判断（S/P 等）——应始终为空，非空即需化学复核
# 判据严格：只统计【微观态未指定】的中心，已冻结中心不计入。
def _refined_stereocentres(ref, molH, Chem):
    """从 3D 坐标感知立体后，返回被嵌入固定、但输入未指定的立体中心。
    返回 [{atom_idx, symbol, class}]；class ∈ {labile_N, symmetric_C, other}。
    ref  : 输入 SMILES 解析所得 RDKit 分子（无坐标，携带已冻结 R/S）。
    molH : 已构造的 3D 主分子（含 H、含坐标）。"""
    m3 = Chem.Mol(molH)
    Chem.AssignStereochemistryFrom3D(m3)
    m3 = Chem.RemoveHs(m3)
    for a in m3.GetAtoms():
        a.SetAtomMapNum(0)                       # atom-map 会扰乱规范排序，先清零
    Chem.AssignStereochemistry(m3, cleanIt=True, force=True)
    r = Chem.Mol(ref)
    Chem.AssignStereochemistry(r, cleanIt=True, force=True)
    match = m3.GetSubstructMatch(r, useChirality=False)
    frozen = ({match[a.GetIdx()] for a in r.GetAtoms()
               if a.GetPropsAsDict().get("_CIPCode") is not None}
              if match else set())
    out = []
    for a in m3.GetAtoms():
        if a.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
            continue
        if a.GetIdx() in frozen:                 # 已冻结中心：不算"嵌入固定"
            continue
        sym = a.GetSymbol()
        q = a.GetFormalCharge()
        deg = a.GetDegree() + a.GetTotalNumHs()
        klass = ("labile_N" if (sym == "N" and q > 0 and deg == 4)
                 else "symmetric_C" if sym == "C"
                 else "other")
        out.append({"atom_idx": int(a.GetIdx()), "symbol": sym, "class": klass})
    return out


def _summarize_refined(manifests):
    """聚合 refined_stereocentres：分子数、逐类中心数、labile-N 若枚举的增量。"""
    from collections import Counter
    per_class = Counter()
    mols_with = 0
    enum_extra = 0
    for m in manifests:
        cs = m.get("refined_stereocentres") or []
        if cs:
            mols_with += 1
        for c in cs:
            per_class[c["class"]] += 1
        n_lab = sum(1 for c in cs if c["class"] == "labile_N")
        if n_lab:
            enum_extra += 2 ** n_lab - 1
    return {
        "n_molecules_with_refined": mols_with,
        "centres_by_class": dict(per_class),
        "labile_N_enumeration_extra_molecules": enum_extra,
        "note": ("refined = microstate-unspecified stereocentres fixed by ETKDG; "
                 "labile_N are physically interconverting (candidate for epimer "
                 "expansion); symmetric_C are the same molecule (record only)"),
    }


def annotate_refined_from_disk(manifests, Chem):
    """为【已在内存但缺 refined_stereocentres 字段】的 manifests 就地补齐该字段，
    数据来源为已写出且经校验的磁盘 SDF（读盘，不重新嵌入）。
    用途：内核仍存活、只想补元数据而不重跑 40 分钟全量构建时。
    返回补齐的条数。"""
    n = 0
    for m in manifests:
        if "refined_stereocentres" in m:
            continue
        sdf = m.get("sdf_path")
        ref_smi = m.get("canonical_isomeric_smiles") or m.get("built_canonical")
        if not sdf or not ref_smi or not Path(sdf).is_file():
            m["refined_stereocentres"] = []
            continue
        mol = Chem.MolFromMolFile(str(sdf), sanitize=True, removeHs=False)
        ref = Chem.MolFromSmiles(str(ref_smi))
        if mol is None or ref is None:
            m["refined_stereocentres"] = []
            continue
        m["refined_stereocentres"] = _refined_stereocentres(ref, mol, Chem)
        n += 1
    return n


def utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _write_sdf_with_props(molH, path, Chem):
    """写 SDF 并保留 molH 上的所有 properties（SDWriter 默认会写非私有 prop）。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    w = Chem.SDWriter(str(path))
    w.write(molH)
    w.close()


def _optimize_geometry(molH, Chem):
    """力场优化，返回使用的方法名。
    先查 MMFF 是否有全参数：有 → MMFF 优化；无 → 跳过（保留 ETKDG 几何），
    【不】回退 UFF（UFF 对某些电荷态不认识，会打 UFFTYPER warning 且几何未必更好）。
    ETKDG 嵌入几何对 docking 起点已足够（docking 自行搜索构象）。"""
    from rdkit.Chem import AllChem
    try:
        if AllChem.MMFFHasAllMoleculeParams(molH):
            AllChem.MMFFOptimizeMolecule(molH)
            return "MMFF94"
        else:
            return "none_ETKDG_only"    # 无 MMFF 参数 → 保留嵌入几何，不用 UFF
    except Exception:
        return "none_optimize_failed"


def _quiet_rdkit():
    """屏蔽 RDKit C++ 日志的上下文管理器（用于嵌入期间的 UFFTYPER 噪声）。
    若当前 RDKit 版本无 BlockLogs，则退化为无操作（不影响正确性）。"""
    try:
        from rdkit import rdBase
        return rdBase.BlockLogs()
    except Exception:
        import contextlib
        return contextlib.nullcontext()


def _robust_embed(smiles, Chem, seeds=(42,)):
    """有界重试嵌入（调用方传 2 种子 × 2 模式 = 4 attempts；maxIterations=0
    即 RDKit 自动 ≈ 10×原子数）。
    成功返回 ((molH, ref, opt_method, attempts), None)；失败返回 (None, reason_category)。
    reason_category ∈ {parse_failure, stereo_geometry_infeasible, embedding_failure}。
    注：stereo_geometry_infeasible 仅表示"该分子含立体标记且四次嵌入均失败"，
    并非断言立体本身不可行；has_stereo 另存为结构化字段供下游判读。"""
    from rdkit.Chem import AllChem
    ref = Chem.MolFromSmiles(str(smiles))
    if ref is None:
        return None, "parse_failure"
    has_stereo = any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
                     for a in ref.GetAtoms()) or bool(ref.GetStereoGroups())
    for atom in ref.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    attempts = 0
    for seed in seeds:
        for use_random in (False, True):
            attempts += 1
            molH = Chem.AddHs(ref)
            p = AllChem.ETKDGv3()
            p.randomSeed = seed
            p.useRandomCoords = use_random
            p.maxIterations = EMBED_MAX_ITERATIONS
            # v0.1.4: 屏蔽 RDKit 内部构象清理产生的 UFFTYPER 噪声。
            # 该警告来自 RDKit 对带电氮无 UFF 原子类型,与本模块的优化策略无关,
            # 且不影响写出的几何。用 attempts 作为结构化替代信号(见下)。
            with _quiet_rdkit():
                ok = (AllChem.EmbedMolecule(molH, p) == 0)
            if ok:
                opt_method = _optimize_geometry(molH, Chem)
                Chem.AssignStereochemistry(molH, cleanIt=True, force=True)
                return (molH, ref, opt_method, attempts), None
    reason = ("stereo_geometry_infeasible" if has_stereo else "embedding_failure")
    return None, reason


def build_one(row, mm, Chem, out_dir, embed_seed=42, dry_run=True,
              smiles_col="isomeric_smiles", name_col="name"):
    """为单行构造 3D 主分子。返回 (status, record)。

    status ∈ {"built", "held", "qc_fail"}。
    不变量：canonical_iso_smiles(built) == canonical_iso_smiles(input isomeric_smiles)。
    """
    name = str(row[name_col])
    smi = str(row[smiles_col])
    q_target = row.get("net_charge")

    # k≥4 保留：不构造 3D
    # v0.1.2: honour Node 01's eligibility decision rather than re-deriving it.
    # Node 01 finalize sets downstream_eligibility (READY / HOLD_STEREO /
    # HOLD_CHARGE_REVIEW) precisely so Node 02 never re-interprets raw flags.
    # Falls back to needs_stereo_review when the column is absent (older tables).
    elig = str(row.get("downstream_eligibility", "") or "")
    if elig.startswith("HOLD"):
        return "held", {
            "name": name, "parent_id": str(row.get("parent_id", "")),
            "reason": f"node01_eligibility:{elig}",
            "isomeric_smiles": smi,
        }
    if not elig and str(row.get("needs_stereo_review")) in ("True", "true", "1"):
        return "held", {
            "name": name, "parent_id": str(row.get("parent_id", "")),
            "reason": "needs_stereo_review (k>=4 held out of docking)",
            "isomeric_smiles": smi,
        }

    # 用有界迭代嵌入（避免几何不可行分子卡死）。稳健重试内含 maxIterations 上限。
    embed_seeds = (embed_seed, embed_seed + 1)
    rescued, fail_reason = _robust_embed(smi, Chem, seeds=embed_seeds)
    if rescued is None:
        # Q3: 区分失败原因（论文统计用）
        return "held", {
            "name": name,
            "parent_id": str(row.get("parent_id", "")),
            "reason": fail_reason,
            "isomeric_smiles": smi,
        }
    molH, ref, opt_method, embed_attempts = rescued

    # 不变量：构造后身份 == 输入身份
    try:
        built_can = mm.canonical_iso_smiles(molH)
        input_can = mm.canonical_iso_smiles(smi)
    except Exception as e:
        return "qc_fail", {"name": name, "stage": "canonicalize",
                           "error": f"{type(e).__name__}: {e}", "isomeric_smiles": smi}

    if built_can != input_can:
        return "qc_fail", {
            "name": name, "stage": "invariant",
            "error": "CanonicalIsoSMILES(built) != CanonicalIsoSMILES(input)",
            "input_canonical": input_can, "built_canonical": built_can,
            "isomeric_smiles": smi,
        }

    # 电荷一致性核查（构造电荷应等于冻结的 net_charge）
    built_charge = Chem.GetFormalCharge(molH)
    charge_ok = True
    try:
        charge_ok = (int(float(q_target)) == built_charge)
    except (TypeError, ValueError):
        charge_ok = None    # 目标缺失，跳过核查

    if charge_ok is False:
        return "qc_fail", {
            "name": name, "stage": "charge",
            "error": f"built charge {built_charge} != frozen net_charge {q_target}",
            "isomeric_smiles": smi,
        }

    # 写 master SDF + identity manifest
    sdf_path = Path(out_dir) / "sdf" / f"{name}.sdf"
    man = mm.identity_manifest(molH, smi, tag=name)
    man["name"] = name
    # Q1: 把 Node 01 关键身份字段写进 SDF properties（SDF 离开 CSV 也不失身份）
    # v0.1.1 (checkagain30): identity_hash【只由化学身份决定】，不含 name。
    #   理由：name 是数据库标签，非化学身份；改名不应改变身份 hash。
    #   行的唯一标识仍是 name（单独存为 property）。
    identity_hash = hashlib.sha256(
        f"{built_can}|{built_charge}".encode()).hexdigest()
    props = {
        "name": name,
        "parent_id": str(row.get("parent_id", "")),
        "identity_hash": identity_hash,
        "identity_hash_basis": "sha256(canonical_isomeric_smiles|formal_charge)",
        "net_charge": str(built_charge),
        "microstate_status": str(row.get("microstate_status", "")),
        "canonical_isomeric_smiles": built_can,
        # Q4 + v0.1.1: conformer / 嵌入溯源
        "conformer_id": "0",
        "generation_method": "ETKDGv3",
        "embedding_seeds": str(list(embed_seeds)),
        "embedding_modes": "standard;useRandomCoords",
        "optimization_method": opt_method,   # MMFF94 / none_ETKDG_only / none_optimize_failed
        "optimization_policy": "MMFF94 if full params else keep ETKDG (no UFF fallback)",
    }
    for k, v in props.items():
        molH.SetProp(k, v)
    # v0.1.3 修复:DRY_RUN=True 时不写 SDF。此前 build_one 无条件写盘,
    # 导致"预览"运行也会在磁盘上留下上千个 SDF,且没有对应 manifest 描述它们。
    if not dry_run:
        _write_sdf_with_props(molH, sdf_path, Chem)
    man["written"] = (not dry_run)
    man["identity_hash"] = identity_hash
    man["sdf_path"] = str(sdf_path)
    man["optimization_method"] = opt_method
    man["embed_attempts"] = embed_attempts   # >1 = 嵌入吃力(原 UFFTYPER 警告的结构化替代)
    man["sdf_sha256"] = file_sha256(sdf_path) if (not dry_run) else None
    man["built_canonical"] = built_can
    man["frozen_net_charge"] = q_target
    man["built_charge"] = built_charge
    man["parent_id"] = row.get("parent_id")
    # v0.2.0: 记录嵌入静默固定的、输入未指定的立体中心（labile_N / symmetric_C）。
    # ref 与 molH 此处均在作用域内，无需重嵌入。
    try:
        man["refined_stereocentres"] = _refined_stereocentres(ref, molH, Chem)
    except Exception as e:
        man["refined_stereocentres"] = []
        man["refined_stereo_error"] = f"{type(e).__name__}: {e}"
    return "built", man


def find_orphaned_parents(df, manifests, name_col="name", parent_col="parent_id",
                          held_records=None):
    """审计检查：找出【所有子代都未成功构建】的母体分子（"孤儿母体"）。

    背景：Node 01 finalize 把一个母体展开成 2^k 个立体异构体；Node 02 可能因
    几何不可行把部分（甚至全部）异构体判为 held。若某母体【全部】子代都没建成，
    该母体实际【整体掉出】docking —— 这是必须显式报告的信息损失，不能沉默。

    返回 (orphan_df, summary)：
      orphan_df: 每个孤儿母体一行（母体、子代总数、held/fail 数、原因集合）
      summary:   计数汇总
    """
    built_names = {m.get("name") for m in manifests}
    # v0.1.3:有意保留(node01_eligibility / needs_stereo_review)不是"意外掉出"。
    # 此前 4 个 HOLD_STEREO + 1 个 HOLD_CHARGE_REVIEW 被误报为孤儿母体。
    intentional = set()
    for h in (held_records or []):
        rr = str(h.get("reason", ""))
        if rr.startswith("node01_eligibility") or rr.startswith("needs_stereo_review"):
            intentional.add(h.get("name"))
    # 按母体分组统计子代成败
    rows = []
    for parent, grp in df.groupby(parent_col):
        children = list(grp[name_col].astype(str))
        n_children = len(children)
        n_built = sum(1 for c in children if c in built_names)
        n_int = sum(1 for c in children if c in intentional)
        if n_built == 0 and n_int < n_children:
            rows.append({
                "parent_id": parent,
                "n_children": n_children,
                "n_built": n_built,
                "n_not_built": n_children - n_built,
                "children": ";".join(children[:8]),
                "n_intentionally_held": n_int,
                "enumerated": n_children > 1,   # True=立体展开过；False=单分子
            })
    summary = {
        "n_orphaned_parents": len(rows),
        "n_orphaned_enumerated": sum(1 for r in rows if r["enumerated"]),
        "n_orphaned_single": sum(1 for r in rows if not r["enumerated"]),
        "orphaned_molecules_lost": sum(r["n_children"] for r in rows),
        "n_intentionally_held": len(intentional),
        "note": "intentionally held parents are NOT counted as orphans",
    }
    return rows, summary


def run_node02(df, mm, Chem, out_dir, embed_seed=42, limit=None, dry_run=True):
    """遍历 finalize 表，构造全部 3D 主分子。返回 (manifests, stats)。

    manifests: 成功构造的 identity manifest 列表（写入 identity_manifests.jsonl）。
    stats: 计数 + failures + held 列表。
    """
    out_dir = Path(out_dir)
    # v0.1.4: 默认 dry_run=True(fail-safe)。忘记传参不会静默写盘。
    print(f"[node02] dry_run={dry_run} -> "
          f"{'不写 SDF(预览)' if dry_run else '将写 SDF 到 ' + str(out_dir / 'sdf')}",
          flush=True)
    manifests, held, failures = [], [], []
    n = 0
    for _, row in df.iterrows():
        if limit is not None and n >= limit:
            break
        n += 1
        status, rec = build_one(row, mm, Chem, out_dir, embed_seed, dry_run=dry_run)
        if status == "built":
            manifests.append(rec)
        elif status == "held":
            held.append(rec)
        else:
            failures.append(rec)

    orphans, orphan_summary = find_orphaned_parents(
        df.head(n) if limit else df, manifests, held_records=held)
    stats = {
        "n_processed": n,
        "n_built": len(manifests),
        "n_held": len(held),
        "n_qc_fail": len(failures),
        "optimization_methods": _count_field(manifests, "optimization_method"),
        "held_reasons": _count_field(held, "reason", split_paren=True),
        "refined_stereo_summary": _summarize_refined(manifests),
        "orphaned_parents": orphan_summary,
        "dry_run_no_sdf_written": dry_run,
        "failures": failures[:50],
        "held_examples": held[:20],
        "orphan_rows": orphans,
    }
    return manifests, held, failures, stats


def _count_field(records, field, split_paren=False):
    """统计某字段的分布（held 原因/优化方法）。"""
    from collections import Counter
    c = Counter()
    for r in records:
        v = str(r.get(field, "unknown"))
        if split_paren:
            v = v.split("(")[0].strip()
        c[v] += 1
    return dict(c)


def write_outputs(manifests, held, failures, stats, out_dir, input_path,
                  dry_run=True):
    """写 identity manifests + node02 manifest。SDF 已在 build_one 中逐个写。"""
    out_dir = Path(out_dir)
    idman_path = out_dir / "identity_manifests.jsonl"
    held_path = out_dir / "node02_held.csv"
    orphan_path = out_dir / "node02_orphaned_parents.csv"
    fail_path = out_dir / "node02_qc_failures.json"
    mani_path = out_dir / "Node02_manifest.json"

    manifest = {
        "module": "node02_master_build",
        "version": NODE02_VERSION,
        "created_at_utc": utc_now(),
        "dry_run": dry_run,
        "input": str(input_path),
        "input_sha256": file_sha256(input_path) if Path(input_path).exists() else None,
        # v0.2.0: 完整可复现血缘
        "rdkit_version": _rdkit_version(),
        "master_module_version": _master_module_version(),
        "upstream_node01": _read_node01_provenance(Path(input_path).parent),
        # v0.2.0: 不变量拆成两条不同的声明 —— 混淆二者正是先前"磁盘重读=1765失配"
        # 误报的根源。internal_gate 是本节点实际强制的东西；disk_reverification 是
        # 外部从 SDF 重新核验时【必须】采用的非对称判据。
        "invariant": {
            "internal_gate": (
                "CanonicalIsoSMILES(built, in-memory) == CanonicalIsoSMILES(input); "
                "标签层，经 master_molecule.canonical_iso_smiles、剥离 atom-map。"
                "qc_fail 计的即此条。"),
            "disk_reverification": (
                "从 SDF 重解析并从 3D 感知立体，会在微观态未指定的中心上追加描述符"
                "（可翻转质子化 N；组成对称 C）。正确核验是非对称的：每个【已指定】"
                "中心恢复出一致 CIP 码；未指定中心可获得由几何决定的构型"
                "（记于 refined_stereocentres，不构成化学身份）。"),
        },
        "policy": {
            "held": "honours Node01 downstream_eligibility (HOLD_STEREO / HOLD_CHARGE_REVIEW); not built",
            "qc_gate": "invariant + charge consistency; failure -> qc_fail (not silent)",
            "chemistry": "NOT decided here; frozen from Node 01 finalize",
            # v0.1.1 (checkagain30): 嵌入与优化的完整可复现记录
            "embedding_method": "ETKDGv3",
            "embedding_seeds": "embed_seed, embed_seed+1 (default 42, 43)",
            "embedding_modes": ["standard", "useRandomCoords"],
            "embedding_max_iterations": EMBED_MAX_ITERATIONS,   # 0 = RDKit 自动 (~10×原子数)
            "optimization_policy": ("MMFF94 only when MMFFHasAllMoleculeParams; "
                                    "otherwise retain ETKDG geometry. "
                                    "No UFF fallback (avoids undefined charge-state "
                                    "parameterization)."),
            "identity_hash_basis": "sha256(canonical_isomeric_smiles|formal_charge)",
        },
        "stats": {k: v for k, v in stats.items()
                  if k not in ("failures", "held_examples", "orphan_rows")},
    }

    # v0.1.4 一致性守卫:若构建阶段是 dry_run(未写 SDF),则不允许写出
    # 声称 n_built 个分子的 manifest —— 否则记录与磁盘不符。
    build_was_dry = bool(stats.get("dry_run_no_sdf_written", False))
    if (not dry_run) and build_was_dry:
        raise RuntimeError(
            "拒绝写盘:构建阶段以 dry_run=True 运行(未写 SDF),"
            "但 write_outputs 被要求以 dry_run=False 落盘。"
            "请先用 dry_run=False 重跑 run_node02,再写 manifest。")
    if not dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "sdf").mkdir(parents=True, exist_ok=True)
        with open(idman_path, "w", encoding="utf-8") as f:
            for m in manifests:
                f.write(json.dumps(m, ensure_ascii=False) + "\n")
        import csv
        if held:
            with open(held_path, "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=list(held[0].keys()))
                w.writeheader(); w.writerows(held)
        # 孤儿母体审计表（所有子代都未构建）
        orphan_rows = stats.get("orphan_rows", [])
        if orphan_rows:
            with open(orphan_path, "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=list(orphan_rows[0].keys()))
                w.writeheader(); w.writerows(orphan_rows)
        fail_path.write_text(json.dumps(failures, indent=2, ensure_ascii=False),
                             encoding="utf-8")
        mani_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False),
                             encoding="utf-8")
    return manifest, {"identity_manifests": idman_path, "held": held_path,
                      "orphaned_parents": orphan_path,
                      "failures": fail_path, "manifest": mani_path}
