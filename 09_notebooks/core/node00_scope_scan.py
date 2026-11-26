# -*- coding: utf-8 -*-
"""
V3_clean/core/scope_scan.py
===========================
Node 00 — Scope reconstruction 的核心算法（只读）。

职责：**发现**候选 CSV、**识别** SMILES 列、按 canonical SMILES **去重**、
**追踪**来源、**审计** pre-docking 过滤列 —— 但【绝不执行任何筛选】。

原则（checkagain22）：
  - 372 是历史候选上限，不是机械送 docking 的名单。
  - N_clean = 通过 pre-docking（ADMET/QM/结构有效性）的唯一结构数。
  - 【不】用旧 Broad_Spectrum_Score / docking rank / worst-case 缩小候选池。
  - 本模块只【报告】哪些列存在、哪些是 pre-docking vs post-docking，交人决策。

不假设列名：先扫描、报告发现的列，再据规则分类。
"""

import re
from pathlib import Path

SCOPE_SCAN_VERSION = "0.2.0"   # v0.2.0: 对齐真实候选池来源（14 run 目录的 docking CSV）

# 候选池的真实来源（来自 2_PRISM_pipeline_part1 cell 14）：
# 每个 run 目录下 step5a_broadspectrum_docking.csv 含 'smiles'+'name'（+docking分数）。
# 这 14 个历史运行目录去重后 = 候选池上限（约 372）。
DEFAULT_RUN_DIRS = [
    "01.1更新第一层后的3*50*30",
    "01.5good-3*50*30",
    "01.6good200-30",
    "01VAE.0后的3*50*30",
    "02VAE.1重训后3*50*30",
    "04VAE.3重训但导航没有重训",
    "05VAE.3重训+导航重训",
    "07VAE_2+U-3_200_50",
    "08VAE_2+U-3_100_30",
    "09VAE_2+U-3_100_30",
    "11VAE_BBS_cycle1",
    "12VAE_BBS_cycle2",
    "13VAE_BBS_cycle3",
    "14VAE_BBS_cycle4",
]
# 每个 run 目录内的候选 docking CSV 文件名
RUN_DOCKING_CSV = "step5a_broadspectrum_docking.csv"
# 顶层参考 CSV（最终入选/静态基线；仅供交叉参照，非候选池来源）
TOPLEVEL_REFERENCE_CSVS = [
    "final_candidates_paper.csv",
    "step5b_final_candidates.csv",
    "step5b_master_summary.csv",
]

# SMILES 列候选名（大小写不敏感匹配）
_SMILES_COL_CANDIDATES = [
    "smiles", "canonical_smiles", "isomeric_smiles", "smi",
    "authoritative_smiles", "source_smiles", "SMILES",
]

# 已知的 post-docking 列（【不可】用于缩小候选池 —— 仅报告，不使用）
_POST_DOCKING_MARKERS = [
    "broad_spectrum_score", "bbs", "dock", "vina", "worst_case",
    "rank", "e_sars", "e_mers", "affinity", "binding",
]

# 已知的 pre-docking 过滤列（可用于定义 N_clean —— 仅报告，交人决策）
_PRE_DOCKING_MARKERS = [
    "qed", "sa_score", "sascore", "mw", "molwt", "logp", "tpsa",
    "admet", "herg", "qm", "gap", "homo", "lumo", "filter_status",
    "lipinski", "veber", "pains", "valid",
]


def find_candidate_csvs(results_dir):
    """只读扫描 results/ 下所有 .csv，返回路径列表（不读内容）。
       泛用兜底；候选池请优先用 find_run_docking_csvs（对齐真实来源）。"""
    results_dir = Path(results_dir)
    if not results_dir.is_dir():
        return []
    return sorted(results_dir.rglob("*.csv"))


def find_run_docking_csvs(results_dir, run_dirs=None, docking_csv=RUN_DOCKING_CSV):
    """候选池的【权威来源】：14 个历史 run 目录下的 docking CSV（含 smiles+name）。

    返回 [(run_name, csv_path, exists_bool), ...]。只读，不读内容。
    对齐 2_PRISM_pipeline_part1 cell 14 的候选池构建逻辑。
    """
    results_dir = Path(results_dir)
    run_dirs = run_dirs if run_dirs is not None else DEFAULT_RUN_DIRS
    out = []
    for run in run_dirs:
        p = results_dir / run / docking_csv
        out.append((run, p, p.exists()))
    return out


def find_toplevel_reference_csvs(results_dir, names=None):
    """顶层参考 CSV（最终入选名单 / 静态基线）。仅交叉参照，非候选池来源。"""
    results_dir = Path(results_dir)
    names = names if names is not None else TOPLEVEL_REFERENCE_CSVS
    out = []
    for n in names:
        p = results_dir / n
        out.append((n, p, p.exists()))
    return out


def detect_smiles_column(columns):
    """从列名列表识别 SMILES 列。返回匹配的列名或 None。"""
    lower = {c.lower(): c for c in columns}
    for cand in _SMILES_COL_CANDIDATES:
        if cand.lower() in lower:
            return lower[cand.lower()]
    # 兜底：任何含 'smiles' 的列
    for c in columns:
        if "smiles" in c.lower():
            return c
    return None


def classify_columns(columns):
    """把列分成 pre-docking / post-docking / other（仅报告，不筛选）。"""
    pre, post, other = [], [], []
    for c in columns:
        cl = c.lower()
        if any(m in cl for m in _POST_DOCKING_MARKERS):
            post.append(c)
        elif any(m in cl for m in _PRE_DOCKING_MARKERS):
            pre.append(c)
        else:
            other.append(c)
    return {"pre_docking": pre, "post_docking": post, "other": other}


def inventory_csv(csv_path, pandas_mod):
    """读一个 CSV 的表头 + 行数（只读），返回清单条目。不做任何筛选。"""
    pd = pandas_mod
    entry = {"path": str(csv_path), "name": Path(csv_path).name}
    try:
        # 只读表头判断列；再读全量算行数（候选 CSV 通常不大）
        df = pd.read_csv(csv_path, dtype=str)
        entry["n_rows"] = len(df)
        entry["columns"] = list(df.columns)
        entry["smiles_column"] = detect_smiles_column(df.columns)
        entry["column_classes"] = classify_columns(df.columns)
        entry["readable"] = True
    except Exception as e:
        entry["readable"] = False
        entry["error"] = f"{type(e).__name__}: {e}"
    return entry


def canonicalize_smiles(smiles, rdkit_chem):
    """规范 SMILES（去 atom map + 去氢），失败返回 None。"""
    try:
        m = rdkit_chem.MolFromSmiles(str(smiles))
        if m is None:
            return None
        for a in m.GetAtoms():
            a.SetAtomMapNum(0)
        return rdkit_chem.MolToSmiles(m, isomericSmiles=True, canonical=True)
    except Exception:
        return None


def dedup_from_run_csvs(run_csv_list, pandas_mod, rdkit_chem, smiles_col="smiles",
                        name_col="name"):
    """从 14 个 run 目录的 docking CSV 去重（对齐真实候选池来源）。

    run_csv_list: find_run_docking_csvs 的返回 [(run, path, exists), ...]
    返回 (unique_df, stats)。按 canonical SMILES 去重，追踪出现在哪些 run。
    ⚠️ 只读 smiles+name，【不】读或使用任何 docking 分数列（防止用 post-docking 缩池）。
    """
    pd = pandas_mod
    records = {}     # canonical -> {runs:set, raw, names:set}
    invalid = []
    total_rows = 0
    per_run_counts = {}

    for run, path, exists in run_csv_list:
        if not exists:
            per_run_counts[run] = None      # 缺失
            continue
        try:
            df = pd.read_csv(path, dtype=str)
        except Exception as e:
            per_run_counts[run] = f"read_error: {type(e).__name__}"
            continue
        if smiles_col not in df.columns:
            per_run_counts[run] = f"no_{smiles_col}_col"
            continue
        cnt = 0
        for _, row in df.iterrows():
            raw = row.get(smiles_col)
            if raw is None or (isinstance(raw, float)):
                continue
            total_rows += 1
            can = canonicalize_smiles(raw, rdkit_chem)
            if can is None:
                invalid.append((run, str(raw)))
                continue
            rec = records.setdefault(can, {"runs": set(), "raw": raw, "names": set()})
            rec["runs"].add(run)
            if name_col in df.columns and row.get(name_col) is not None:
                rec["names"].add(str(row.get(name_col)))
            cnt += 1
        per_run_counts[run] = cnt

    rows = []
    for can, rec in records.items():
        rows.append({
            "canonical_smiles": can,
            "n_runs": len(rec["runs"]),
            "first_run": sorted(rec["runs"])[0],
            "raw_smiles_example": rec["raw"],
            "source_names": ";".join(sorted(rec["names"])[:5]),
        })
    cols = ["canonical_smiles", "n_runs", "first_run", "raw_smiles_example", "source_names"]
    unique_df = (pd.DataFrame(rows).sort_values("canonical_smiles").reset_index(drop=True)
                 if rows else pd.DataFrame(columns=cols))
    stats = {
        "total_smiles_rows_scanned": total_rows,
        "unique_canonical": len(records),
        "invalid_smiles": len(invalid),
        "invalid_examples": invalid[:10],
        "per_run_counts": per_run_counts,
    }
    return unique_df, stats


def dedup_candidates(inventory, pandas_mod, rdkit_chem, source_filter=None):
    """跨所有可读 CSV 的 SMILES 去重（canonical）。
       返回 (unique_df, stats)。不做任何 pre/post-docking 筛选。

    unique_df 列：canonical_smiles, first_source_csv, source_count,
                  raw_smiles_example, invalid（能否解析）。
    """
    pd = pandas_mod
    records = {}       # canonical -> {sources:set, raw_example, valid}
    invalid = []
    total_rows = 0

    for entry in inventory:
        if not entry.get("readable"):
            continue
        smi_col = entry.get("smiles_column")
        if not smi_col:
            continue
        if source_filter and entry["name"] not in source_filter:
            continue
        df = pd.read_csv(entry["path"], dtype=str)
        for raw in df[smi_col].dropna():
            total_rows += 1
            can = canonicalize_smiles(raw, rdkit_chem)
            if can is None:
                invalid.append((entry["name"], raw))
                continue
            rec = records.setdefault(can, {"sources": set(), "raw": raw})
            rec["sources"].add(entry["name"])

    rows = []
    for can, rec in records.items():
        rows.append({
            "canonical_smiles": can,
            "source_count": len(rec["sources"]),
            "first_source_csv": sorted(rec["sources"])[0],
            "raw_smiles_example": rec["raw"],
        })
    unique_df = pd.DataFrame(rows).sort_values("canonical_smiles").reset_index(drop=True) \
        if rows else pd.DataFrame(columns=["canonical_smiles", "source_count",
                                           "first_source_csv", "raw_smiles_example"])
    stats = {
        "total_smiles_rows_scanned": total_rows,
        "unique_canonical": len(records),
        "invalid_smiles": len(invalid),
        "invalid_examples": invalid[:10],
    }
    return unique_df, stats
