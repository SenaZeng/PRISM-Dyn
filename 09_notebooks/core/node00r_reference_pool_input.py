#!/usr/bin/env python3
"""
build_reference_pool_input.py — 为参考小池构造 Node 01 的输入表

参考化合物没有"13 个历史运行"的来源，故【不经 Node 00】。
本脚本手工构造一张与 all_candidates_v3_start.csv 【同 schema】的 3 行表，
使其可从 Node 01 进入，走与 31 个候选完全相同的代码路径。

设计要点：
  · schema 从真实文件读取，不假设列数（32 还是 35）
  · 参考物无 pIC50 / admet_* / PySCF_* 等属性 -> 留空
  · provenance 列标记为 reference_panel，不伪造历史来源
  · pool_class（若存在）赋 'reference'，使其在选择逻辑中可被排除

输出：
  results/00_scope_ref/reference_pool_v3_start.csv        (3 行，同 schema)
  results/00_scope_ref/reference_pool_identity.json       (SMILES/CID/IC50/出处/预期断言)

只读 all_candidates_v3_start.csv 取 schema；只写上述两个新文件。

用法:  cd /mnt/hpc/home/wanghongmei/lang2 && python3 core/build_reference_pool_input.py
"""
import csv
import json
import hashlib
import datetime
from pathlib import Path

ROOT = Path("/mnt/hpc/home/wanghongmei/lang2")
SCHEMA_SRC = ROOT / "results/00_scope/all_candidates_v3_start.csv"
OUT_DIR = ROOT / "results/00_scope_ref"
OUT_CSV = OUT_DIR / "reference_pool_v3_start.csv"
OUT_JSON = OUT_DIR / "reference_pool_identity.json"


# ── 参考面板（方案 A）。SMILES 全部经 CID 解析并逐项核验 ────────────────────
# 注意 Ensitrelvir 的 SMILES 含两个反斜杠（C=N 的 E 构型）
# -> Python 中必须为原始字符串 r"..."，否则 \N 触发 unicode 转义错误
REFERENCES = [
    {
        "name": "ref_AVI_1084",
        "common_name": "AVI-1084",
        "enamine_id": "Z3535317212",
        "smiles": "C1CN(C(=O)N(C1=O)CC2=CC=CNC2=O)CC3=CC=CS3",
        "pubchem_cid": 167818027,
        "expected_formula": "C15H15N3O3S",
        "expected_mw": 317.4,
        "expected_net_charge": 0,
        "expected_stereocentres": 0,
        "ic50_uM": 29.0,
        "ic50_source": "Detomasi et al., Sci Adv 11, eadt7836 (2025) — 862M docking screen hit",
        "mechanism": "noncovalent",
        "role": "weak anchor",
    },
    {
        "name": "ref_X77",
        "common_name": "X77",
        "enamine_id": None,
        "smiles": "CC(C)(C)C1=CC=C(C=C1)N([C@H](C2=CN=CC=C2)C(=O)NC3CCCCC3)C(=O)C4=CN=CN4",
        "pubchem_cid": 145998279,
        "expected_formula": "C27H33N5O2",
        "expected_mw": 459.6,
        "expected_net_charge": 0,
        "expected_stereocentres": 1,          # (1R)，RDKit 实测 [(11, 'R')]
        "ic50_uM": 2.8,
        "ic50_source": "6W63 co-crystal report; 另有 4.1 μM 的报道值（J Med Chem 综述）",
        "mechanism": "noncovalent",
        "role": "mid-range; 6W63 co-crystal ligand (crystallographic pose available)",
    },
    {
        "name": "ref_ensitrelvir",
        "common_name": "Ensitrelvir (S-217622)",
        "enamine_id": None,
        "smiles": r"Cn1cnc(CN2C(=O)N(Cc3cc(F)c(F)cc3F)C(=N\c3cc4cn(C)nc4cc3Cl)\NC2=O)n1",
        "pubchem_cid": 162533924,
        "expected_formula": "C22H17ClF3N9O2",
        "expected_mw": 531.9,
        "expected_net_charge": 0,
        "expected_stereocentres": 0,          # 立体信息在 C=N 双键几何(E)，非手性中心
        "ic50_uM": 0.013,
        "ic50_source": "PubChem CID 162533924 / 多来源；另有 0.049 μM 的实测值 (Commun Biol 2023)",
        "mechanism": "noncovalent",
        "role": "potent end; clinical-grade",
    },
]

# provenance 列的填充值：不伪造历史来源
PROVENANCE_FILL = {
    "source_runs": "reference_panel",
    "source_count": 1,
    "run_types": "reference_panel",
    "is_bbs_guided": 0,
    "primary_source": "reference_panel",
    "pool_class": "reference",
    "pool_class_reason": "external_reference_compound",
    "qm_status": "none",
}


def sha256_str(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def preflight():
    """开跑前的自检：反斜杠、共价弹头、必需字段。"""
    problems = []
    for r in REFERENCES:
        smi = r["smiles"]
        if "C#C" in smi:
            problems.append(f"{r['common_name']}: SMILES 含 C#C —— 疑为共价弹头，不可用")
        if r["common_name"].startswith("Ensitrelvir") and smi.count("\\") != 2:
            problems.append(f"{r['common_name']}: 反斜杠数为 {smi.count(chr(92))}，应为 2 "
                            f"—— C=N 的 E 构型已丢失（是否忘了 r\"\" ？）")
        for k in ("name", "smiles", "pubchem_cid", "ic50_uM", "ic50_source"):
            if r.get(k) in (None, ""):
                problems.append(f"{r['common_name']}: 缺字段 {k}")
    return problems


def main():
    print("=" * 74)
    print("  参考小池 · Node 01 输入表构造")
    print("=" * 74)

    # 1) 自检
    print("\n-- 开跑前自检 --")
    problems = preflight()
    if problems:
        for p in problems:
            print("  ✗", p)
        raise SystemExit("自检未通过，已中止；未写任何文件。")
    print(f"  ✓ {len(REFERENCES)} 个化合物：无 C#C、反斜杠完整、字段齐全")

    # 2) 读真实 schema —— 不假设列数
    if not SCHEMA_SRC.is_file():
        raise SystemExit(f"找不到 schema 来源: {SCHEMA_SRC}")
    with open(SCHEMA_SRC, newline="", encoding="utf-8") as f:
        header = next(csv.reader(f))
    print(f"\n-- schema 来源 --")
    print(f"  {SCHEMA_SRC.name}: {len(header)} 列")
    print(f"  前 8 列: {header[:8]}")

    # 3) 逐列决定填什么
    id_cols = {"canonical_smiles", "name"}
    filled, blanked, unknown = [], [], []
    rows = []
    for r in REFERENCES:
        row = {}
        for col in header:
            if col == "canonical_smiles":
                row[col] = r["smiles"]
            elif col == "name":
                row[col] = r["name"]
            elif col in PROVENANCE_FILL:
                row[col] = PROVENANCE_FILL[col]
            else:
                row[col] = ""          # 属性列：参考物无此数据，留空
        rows.append(row)

    for col in header:
        if col in id_cols:
            filled.append(col)
        elif col in PROVENANCE_FILL:
            filled.append(col)
        else:
            blanked.append(col)

    print(f"\n-- 列处理 --")
    print(f"  填充 {len(filled)} 列: {filled}")
    print(f"  留空 {len(blanked)} 列（参考物无该属性）")
    if len(blanked) > 0:
        print(f"    {blanked[:12]}{' …' if len(blanked) > 12 else ''}")

    missing_fill = [k for k in PROVENANCE_FILL if k not in header]
    if missing_fill:
        print(f"\n  ⚠ 以下预填列在真实 schema 中不存在，已忽略: {missing_fill}")
        print("     （若 pool_class 不在其中，说明 backfill 未写入该 CSV，需另行确认）")

    # 4) 落盘
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        w.writerows(rows)
    print(f"\n已写 {OUT_CSV.relative_to(ROOT)}  ({len(rows)} 行 × {len(header)} 列)")

    # 5) 身份冻结记录
    freeze = {
        "module": "reference_pool_input",
        "version": "0.1.0",
        "created_at_utc": datetime.datetime.utcnow().isoformat() + "Z",
        "panel": "Plan A — three non-covalent references",
        "schema_source": str(SCHEMA_SRC.relative_to(ROOT)),
        "schema_n_columns": len(header),
        "columns_filled": filled,
        "columns_blank": blanked,
        "resolution_limit_kcal_mol": 2.63,
        "resolution_limit_note": (
            "取自 V3 候选队列三 replica 的 σ_PBSA (2 SE)。这是流程【内部】的精度估计，"
            "故据此评判参考对照是内部一致性陈述，不是外部验证。"
            "参考化合物自身的 σ_PBSA 须在其三重轨迹跑完后另行核验（见 expected_checks）。"
        ),
        "compounds": [
            {k: v for k, v in r.items()} | {"smiles_sha256": sha256_str(r["smiles"])}
            for r in REFERENCES
        ],
        "expected_checks": {
            "node01_stereo_enumeration": {
                "ref_AVI_1084": "不应枚举出异构体（无立体中心）",
                "ref_X77": "不应枚举出异构体（1 个中心，已定义为 R）",
                "ref_ensitrelvir": "C=N 应按 E 处理；若判为未定义，说明反斜杠丢失",
            },
            "node02_identity_invariant": "三者皆应 built，qc_fail = 0",
            "node03_qc": "三者皆应 built，qc_fail = 0；共晶配体可能含本池未遇过的官能团",
            "node05_invariants": "六道不变量皆过，尤其 I2b（图身份）",
            "sigma_ref_consistency": (
                "phase3 跑完后，用三个参考化合物自身的三重轨迹计算 σ_PBSA，"
                "确认与候选队列的 2.63 kcal/mol 阈值一致，再断言可分辨性。零额外成本。"
            ),
        },
        "downstream_policy": {
            "entry_point": "Node 01（不经 Node 00 —— 参考物无历史运行来源）",
            "compass_treatment": (
                "投影而非并入：以 31 个候选计算各轴 median/MAD，"
                "再把参考化合物投影到该标度；不参与定义标度"
            ),
            "output_arm": "md_solv_v3_ref",
        },
    }
    OUT_JSON.write_text(json.dumps(freeze, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"已写 {OUT_JSON.relative_to(ROOT)}")
    print(f"  sha256: {sha256_str(OUT_JSON.read_text(encoding='utf-8'))}")

    print("\n" + "=" * 74)
    print("下一步：确认 Node 01 的 finalize(df, ...) 需要哪些列，")
    print("      以及本表留空的属性列是否会被容忍。")
    print("=" * 74)


if __name__ == "__main__":
    main()
