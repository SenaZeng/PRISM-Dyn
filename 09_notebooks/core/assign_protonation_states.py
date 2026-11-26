#!/usr/bin/env python3
"""
core/assign_protonation_states.py  (v4 · 吡啶/吡咯拆分 + pyridine_rich 软标记 + 完整候选层级)

决策 3B · 微观态表 A 的 pH-7.4 定态候选生成器与分流工具（非最终电荷裁决器）。

修订历史：
  [v2] 显式 SMARTS 排除酰胺/内酰胺/磺酰胺/脲；普通酚中性、活化酚复核；难例多候选；单测门控。
  [v3 · checkagain6] 源电荷不静默忽略；conf/复核判据统一；芳基胺分类；酚盐识别；
       列名改电荷层级；互斥三层=28 断言。
  [v4 · checkagain7 + checkagain8]
   1. candidate_charge_levels 包含触发复核的态（003:−1;0；008:0;1；015:0;1）。
   2. 芳香氮拆分：pyridine_like_N（无H，低优先碱）vs pyrrolic_N（[nH]，非碱）。
   3. pyridine_rich 软标记（≥4 吡啶型 N）：不改分层（选项 C），记录可审计的中性假设，
      设 secondary_pKa_review=CONDITIONAL；仅在进 Top-N / Tier1 / 邻强给电子 / 状态敏感
      翻转时升级 Epik。【不】并入 needs_pKa_review（否则变成选项 A）。
   4. 报告措辞：'自动高置信可冻结' -> '规则建议中性主态（可审计假设）'。
   5. 新增 3 单测（酚盐层级[-1,0]、芳基胺层级[0,1]、铵盐层级[0,1]）+ 软标记一致性断言。

输出 -> results/V3_microstate_tableA_draft.csv  (候选表，非最终真值)
        results/V3_protonation_report.md

用法:
    conda activate md_env
    python core/assign_protonation_states.py       # 单元测试门控 -> 处理 28 个
"""
import sys, csv, json
from pathlib import Path

ROOT    = Path("/mnt/hpc/home/wanghongmei/lang2")
MAP_CSV = ROOT / "results/TableS5_Dyn28_to_Gen36_Tanimoto_Overlap.csv"
OUT_CSV = ROOT / "results/V3_microstate_tableA_draft.csv"
OUT_MD  = ROOT / "results/V3_protonation_report.md"
PH = 7.4

SENSITIVITY = {"pool_mol_017", "pool_mol_022", "pool_mol_027"}
TIER1_OLD   = "pool_mol_027"
PYRIDINE_RICH_THRESHOLD = 4

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")
except ImportError:
    sys.exit("需要 rdkit（conda activate md_env）")

# ---- 显式 SMARTS 规则 ----
AMIDE_N          = Chem.MolFromSmarts("[N;X3][C,S](=[O,S])")
LACTAM_N         = Chem.MolFromSmarts("[N;r][C;r](=O)")
SULFONAMIDE_N    = Chem.MolFromSmarts("[N][S](=O)(=O)")
UREA_CARBAMATE_N = Chem.MolFromSmarts("[N][C](=O)[N,O]")
GUANIDINE_AMIDINE = Chem.MolFromSmarts("[N]=[C]([N])[N,C]")
PHENOL           = Chem.MolFromSmarts("[O;H1]-[c]")
PHENOXIDE        = Chem.MolFromSmarts("[O-]-[c]")
CARBOXYL         = Chem.MolFromSmarts("[C](=O)[O;H1]")
EWG_NITRO        = Chem.MolFromSmarts("[N+](=O)[O-]")


def match_atoms(mol, patt):
    return set(m[0] for m in mol.GetSubstructMatches(patt)) if patt else set()


def source_charge_info(mol):
    q = Chem.GetFormalCharge(mol)
    charged = [dict(atom_idx=a.GetIdx(), element=a.GetSymbol(),
                    formal_charge=a.GetFormalCharge())
               for a in mol.GetAtoms() if a.GetFormalCharge() != 0]
    return q, charged


def classify_nitrogens(mol):
    """排除酰胺/内酰胺/磺酰胺/脲；芳香氮拆吡啶型/吡咯型；再分芳基胺、伯/仲/叔胺。"""
    nonbasic = (match_atoms(mol, AMIDE_N) | match_atoms(mol, LACTAM_N)
                | match_atoms(mol, SULFONAMIDE_N) | match_atoms(mol, UREA_CARBAMATE_N))
    guan = match_atoms(mol, GUANIDINE_AMIDINE)
    result = dict(nonbasic_N=[], pyridine_like_N=[], pyrrolic_N=[],
                  aryl_substituted_amine=[], primary_amine=[], secondary_amine=[],
                  tertiary_amine=[], guanidine_amidine=[])
    for a in mol.GetAtoms():
        if a.GetSymbol() != "N":
            continue
        idx = a.GetIdx()
        if idx in nonbasic:
            result["nonbasic_N"].append(idx)
        elif idx in guan:
            result["guanidine_amidine"].append(idx)
        elif a.GetIsAromatic():
            # checkagain7: 吡啶型(无H,低优先碱) vs 吡咯型([nH],非碱)
            if a.GetTotalNumHs() > 0:
                result["pyrrolic_N"].append(idx)
            else:
                result["pyridine_like_N"].append(idx)
        elif any(nb.GetIsAromatic() for nb in a.GetNeighbors()):
            result["aryl_substituted_amine"].append(idx)
        else:
            nH = a.GetTotalNumHs()
            if nH >= 2:   result["primary_amine"].append(idx)
            elif nH == 1: result["secondary_amine"].append(idx)
            else:         result["tertiary_amine"].append(idx)
    return result


def phenol_status(mol):
    phenols = match_atoms(mol, PHENOL)
    phenoxides = match_atoms(mol, PHENOXIDE)
    n_phenol, n_phenoxide = len(phenols), len(phenoxides)
    if not phenols and not phenoxides:
        return 0, 0, 0
    has_nitro = bool(match_atoms(mol, EWG_NITRO))
    n_halogen = sum(1 for a in mol.GetAtoms() if a.GetSymbol() in ("F", "Cl", "Br", "I"))
    activated = n_phenol if (has_nitro or n_halogen >= 3) else 0
    return n_phenol, activated, n_phenoxide


def analyze(mol_id, smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return dict(mol_id=mol_id, source_smiles=smi, status="SMILES解析失败",
                    source_formal_charge="", source_charged_atoms="",
                    suggested_formal_charge="", candidate_charge_levels="",
                    candidate_protonation_count="", n_primary_amine=0, n_secondary_amine=0,
                    n_tertiary_amine=0, n_aryl_substituted_amine=0, n_pyridine_like_N=0,
                    n_pyrrolic_N=0, n_nonbasic_N=0, n_guanidine_amidine=0, n_phenol=0,
                    n_phenol_activated=0, n_phenoxide=0, excluded_sites="",
                    needs_pKa_review="YES", is_hard_case="", is_old_tier1="",
                    pyridine_rich="", pyridine_neutral_assumption="", secondary_pKa_review="",
                    confidence="parse_fail", confirmed_by_human="")

    src_q, src_charged = source_charge_info(mol)
    cls = classify_nitrogens(mol)
    phe_n, phe_act, phe_ox = phenol_status(mol)

    n_basic = (len(cls["primary_amine"]) + len(cls["secondary_amine"])
               + len(cls["tertiary_amine"]))
    n_guan = len(cls["guanidine_amidine"])
    n_aryl_amine = len(cls["aryl_substituted_amine"])
    n_pyridine = len(cls["pyridine_like_N"])
    has_cooh = bool(match_atoms(mol, CARBOXYL))

    # ---- 候选电荷层级：含触发复核的态（checkagain7 问题1）----
    level_values = {0}
    for k in range(1, n_basic + n_guan + 1):
        level_values.add(k)
    if has_cooh:
        level_values.add(-1)
    level_values.add(src_q)                       # 源电荷态
    if phe_ox >= 1:
        level_values.add(src_q + 1)               # 酚盐 -> 中性酚
    if n_aryl_amine >= 1 and src_q == 0:
        level_values.add(1)                       # 中性芳基胺 -> 可能 +1
    levels_sorted = sorted(level_values)

    # 位点标签（与层级对应，仅描述计数）
    def level_label(q):
        if q < 0:  return f"{q}_deprotonated"
        if q == 0: return "all_neutral"
        return f"{q}_basic_protonated"

    excluded = f"内酰胺/酰胺/磺酰胺/脲 N ({len(cls['nonbasic_N'])}个) 永久排除质子化"

    # ---- pyridine_rich 软标记（checkagain8 选项 C）----
    pyridine_rich = (n_pyridine >= PYRIDINE_RICH_THRESHOLD)

    # ---- 统一复核判据（注意：pyridine_rich 不并入，否则变成选项 A）----
    hard = (n_basic >= 2) or (n_guan >= 1) or (phe_act >= 1)
    needs_review = (
        (mol_id in SENSITIVITY) or (n_basic >= 2) or (n_guan >= 1)
        or (phe_act >= 1) or (phe_ox >= 1) or (n_aryl_amine >= 1)
        or (n_basic == 1) or (src_q != 0)
    )

    # ---- 建议主态 ----
    if src_q != 0:
        suggested = ""
    elif n_basic == 0 and n_guan == 0 and n_aryl_amine == 0:
        suggested = -1 if has_cooh else 0
    else:
        suggested = ""

    confidence = "review" if needs_review else "high"

    return dict(
        mol_id=mol_id, source_smiles=smi,
        source_formal_charge=src_q,
        source_charged_atoms=json.dumps(src_charged, ensure_ascii=False) if src_charged else "",
        suggested_formal_charge=suggested,
        candidate_charge_levels=";".join(str(q) for q in levels_sorted),
        candidate_protonation_count=";".join(level_label(q) for q in levels_sorted),
        n_primary_amine=len(cls["primary_amine"]),
        n_secondary_amine=len(cls["secondary_amine"]),
        n_tertiary_amine=len(cls["tertiary_amine"]),
        n_aryl_substituted_amine=n_aryl_amine,
        n_pyridine_like_N=n_pyridine,
        n_pyrrolic_N=len(cls["pyrrolic_N"]),
        n_nonbasic_N=len(cls["nonbasic_N"]),
        n_guanidine_amidine=n_guan,
        n_phenol=phe_n, n_phenol_activated=phe_act, n_phenoxide=phe_ox,
        excluded_sites=excluded,
        needs_pKa_review=("YES" if needs_review else "no"),
        is_hard_case=("YES" if hard else ""),
        is_old_tier1=("YES" if mol_id == TIER1_OLD else ""),
        pyridine_rich=("YES" if pyridine_rich else ""),
        pyridine_neutral_assumption=("provisional_neutral_at_pH7.4" if pyridine_rich else ""),
        secondary_pKa_review=("CONDITIONAL" if pyridine_rich else ""),
        confidence=confidence,
        confirmed_by_human="",
    )


# ================= 单元测试（门控）=================
def run_unit_tests():
    tests = [
        ("乙酰胺", "CC(=O)N", lambda r: r["n_nonbasic_N"] == 1 and r["n_primary_amine"] == 0,
         "酰胺 N 不得当碱性胺"),
        ("β-内酰胺", "O=C1CCN1", lambda r: r["n_nonbasic_N"] == 1 and r["n_primary_amine"] == 0,
         "环内酰胺 N 不得质子化"),
        ("乙胺", "CCN", lambda r: r["n_primary_amine"] == 1, "伯胺应识别为质子化候选"),
        ("苯酚", "c1ccccc1O", lambda r: r["n_phenol"] == 1 and r["n_phenol_activated"] == 0
                                        and r["suggested_formal_charge"] == 0, "普通酚默认中性"),
        ("乙酸", "CC(=O)O", lambda r: "-1" in str(r["candidate_charge_levels"]),
         "羧酸应生成去质子化候选"),
        ("吡啶", "c1ccncc1", lambda r: r["n_pyridine_like_N"] == 1 and r["n_pyrrolic_N"] == 0,
         "吡啶归 pyridine_like（低优先候选）"),
        ("吡咯", "c1cc[nH]c1", lambda r: r["n_pyrrolic_N"] == 1 and r["n_pyridine_like_N"] == 0,
         "吡咯归 pyrrolic（非碱）"),
        ("对硝基苯酚", "Oc1ccc(cc1)[N+](=O)[O-]",
         lambda r: r["n_phenol"] == 1 and r["n_phenol_activated"] >= 1
                   and r["needs_pKa_review"] == "YES", "活化酚进复核"),
        # checkagain7 新增 3 个：候选层级须含触发复核的态
        ("显式酚盐", "[O-]c1ccccc1",
         lambda r: r["source_formal_charge"] == -1 and r["needs_pKa_review"] == "YES"
                   and set(str(r["candidate_charge_levels"]).split(";")) == {"-1", "0"},
         "源酚盐：候选层级=[-1,0]"),
        ("显式铵盐", "C[NH2+]C",
         lambda r: r["source_formal_charge"] == 1 and r["needs_pKa_review"] == "YES"
                   and set(str(r["candidate_charge_levels"]).split(";")) == {"0", "1"},
         "源铵盐：候选层级=[0,1]"),
        ("芳基胺", "Nc1ccncc1",
         lambda r: r["n_aryl_substituted_amine"] == 1 and r["needs_pKa_review"] == "YES"
                   and set(str(r["candidate_charge_levels"]).split(";")) == {"0", "1"},
         "芳基胺：候选层级=[0,1]，进复核"),
    ]
    print("=== 分类器单元测试（checkagain5+6+7，共 %d 个）===" % len(tests))
    all_pass = True
    for name, smi, check, desc in tests:
        r = analyze(name, smi)
        try:
            ok = bool(check(r))
        except Exception:
            ok = False
        if not ok:
            all_pass = False
            print(f"  [FAIL] {name}: {desc}")
            print(f"         src_q={r.get('source_formal_charge')} levels={r.get('candidate_charge_levels')!r} "
                  f"pyridine={r.get('n_pyridine_like_N')} pyrrolic={r.get('n_pyrrolic_N')} "
                  f"aryl_amine={r.get('n_aryl_substituted_amine')} review={r.get('needs_pKa_review')}")
        else:
            print(f"  [PASS] {name}: {desc}")
    print(f"\n单元测试{'全部通过 ✓' if all_pass else '有失败 ✗ —— 修复后再处理 28 个'}\n")
    return all_pass


def load_smiles():
    m = {}
    with open(MAP_CSV, newline="") as f:
        rd = csv.reader(f)
        header = [h.strip().lstrip("\ufeff") for h in next(rd)]
        i_id, i_smi = header.index("mol_id"), header.index("smiles")
        for row in rd:
            if len(row) > max(i_id, i_smi):
                m[row[i_id].strip()] = row[i_smi].strip()
    return m


def main():
    if not run_unit_tests():
        sys.exit("⚠️ 单元测试未全过，中止。请修复分类器后重跑。")

    SMILES = load_smiles()
    mols = sorted(SMILES.keys(), key=lambda x: (len(x), x))
    rows = [analyze(m, SMILES[m]) for m in mols]

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    keys = ["mol_id", "source_smiles", "source_formal_charge", "source_charged_atoms",
            "suggested_formal_charge", "candidate_charge_levels", "candidate_protonation_count",
            "n_primary_amine", "n_secondary_amine", "n_tertiary_amine",
            "n_aryl_substituted_amine", "n_pyridine_like_N", "n_pyrrolic_N", "n_nonbasic_N",
            "n_guanidine_amidine", "n_phenol", "n_phenol_activated", "n_phenoxide",
            "excluded_sites", "needs_pKa_review", "is_hard_case", "is_old_tier1",
            "pyridine_rich", "pyridine_neutral_assumption", "secondary_pKa_review",
            "confidence", "confirmed_by_human"]
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    # ---- 互斥三层 + 软标记一致性断言 ----
    hard = [r for r in rows if r.get("is_hard_case") == "YES"]
    review_only = [r for r in rows if r.get("needs_pKa_review") == "YES"
                   and r.get("is_hard_case") != "YES"]
    auto = [r for r in rows if r.get("needs_pKa_review") != "YES"]
    assert len(auto) + len(review_only) + len(hard) == len(rows), "分层计数不互斥"
    assert all(not (r["confidence"] == "review" and r["needs_pKa_review"] != "YES")
               for r in rows), "conf=review 但 needs_pKa_review!=YES 不一致"

    pyr_rich = [r for r in rows if r.get("pyridine_rich") == "YES"]
    pyr_rich_in_review = [r for r in pyr_rich if r.get("needs_pKa_review") == "YES"]
    pyr_rich_in_auto = [r for r in pyr_rich if r.get("needs_pKa_review") != "YES"]
    # 软标记不应【单独因 pyridine_rich】把分子推入复核（选项 C，非 A）

    # ---- checkagain9: 所有数字从实际列表动态生成，杜绝硬编码与 CSV 脱节 ----
    mandatory_epik = review_only + hard
    mandatory_epik_ids = [r["mol_id"] for r in mandatory_epik]
    conditional_pyridine = pyr_rich_in_auto
    conditional_pyridine_ids = [r["mol_id"] for r in conditional_pyridine]
    assert len(pyr_rich) == len(pyr_rich_in_review) + len(pyr_rich_in_auto), \
        "pyridine_rich 计数不自洽"
    assert len(mandatory_epik) == len(review_only) + len(hard), "强制 Epik 计数不自洽"

    src_charged = [r for r in rows if r.get("source_formal_charge") not in (0, "", None)]

    L = ["# 微观态表 A · pH-7.4 候选（草稿第4轮，非最终真值）\n",
         "⚠️ **候选生成器 + 分流工具，非最终电荷裁决器。**\n",
         "本轮修（checkagain7+8）：候选层级含触发复核的态；芳香氮拆吡啶型/吡咯型；"
         "加 pyridine_rich 软标记（选项 C，不改分层）；措辞改'规则建议中性主态'。\n",
         "## 互斥三层（计数严格 = 28）",
         f"- 规则建议中性主态（无碱/无酸/普通酚，源电荷=0）：{len(auto)} 个",
         f"  （注：这是'规则未发现须复核的强电离位点'，**非** pKa 已证明；可作暂定微观态）",
         f"- 普通复核（单一真碱 / 芳基胺 / 源带电 / 酚盐）：{len(review_only)} 个",
         f"- 难例（≥2真碱 / 胍脒 / 活化酚）：{len(hard)} 个",
         f"- 合计：{len(auto)+len(review_only)+len(hard)} = {len(rows)} ✓\n"]

    L.append(f"## 多吡啶软标记（pyridine_rich，选项 C）")
    L.append(f"- 共 {len(pyr_rich)} 个：{[r['mol_id'] for r in pyr_rich]}")
    L.append(f"- 其中强制复核组内 {len(pyr_rich_in_review)} 个，规则中性组内 {len(pyr_rich_in_auto)} 个")
    L.append(f"- 采用可审计假设：普通未活化吡啶型氮在 pH 7.4 下暂定以中性态为主"
             f"（如需入 SI/审计材料，应补文献 DOI 或保留此定性表述，不写单一数值）")
    L.append(f"- secondary_pKa_review=CONDITIONAL：仅在进 Top-N / Tier1 / 邻强给电子 / "
             f"状态敏感翻转时升级 Epik\n")

    if src_charged:
        L.append("## ⚠️ 源 SMILES 已带电荷（不得静默中和/保留，必须复核）")
        for r in src_charged:
            L.append(f"- {r['mol_id']}: source_formal_charge={r['source_formal_charge']}, "
                     f"候选层级={r['candidate_charge_levels']}, 带电原子={r['source_charged_atoms']}")
        L.append("")

    L.append("## 难例明细（电荷层级，非具体微观态）")
    for r in hard:
        L.append(f"### {r['mol_id']}" + (" ⚠️ 旧Tier1#1" if r.get("is_old_tier1") == "YES" else ""))
        L.append(f"- 氮：伯胺 {r['n_primary_amine']} / 仲胺 {r['n_secondary_amine']} / "
                 f"叔胺 {r['n_tertiary_amine']} / 芳基胺 {r['n_aryl_substituted_amine']} / "
                 f"吡啶型 {r['n_pyridine_like_N']} / 吡咯型 {r['n_pyrrolic_N']} / "
                 f"非碱 {r['n_nonbasic_N']} / 胍脒 {r['n_guanidine_amidine']}")
        L.append(f"- **候选电荷层级**：{r['candidate_charge_levels']}（电荷计数，非具体微观态）")
        L.append(f"- 排除：{r['excluded_sites']}")
        if r.get("is_old_tier1") == "YES":
            L.append("- ⚠️ **027：真伯胺上 0/+1/+2 层级，内酰胺永久排除；具体微观态由 Epik 定**")
        L.append("")

    L.append("## 送 Epik 前的顺序")
    L.append("1. 源电荷已记录+候选层级已含 ✓  2. conf/复核已统一 ✓  3. 酚盐/芳基胺/吡啶拆分 ✓")
    L.append("4. 本轮已重跑 ✓  5. 互斥分层=28 已断言 ✓  6. pyridine_rich 软标记 ✓")
    L.append(f"7. → 强制 Epik {len(mandatory_epik_ids)} 个：{mandatory_epik_ids}；"
             f"pyridine_rich 共 {len(pyr_rich)} 个，其中 {len(pyr_rich_in_review)} 个已在强制组，"
             f"另有 {len(conditional_pyridine_ids)} 个按条件升级：{conditional_pyridine_ids}")
    L.append("\n未过 pKa 复核者：禁止进入正式长程 MD；可进入低成本微观态敏感性分析。")
    OUT_MD.write_text("\n".join(L))

    print(f"完成（草稿第4轮）。互斥三层：规则中性 {len(auto)} / 普通复核 {len(review_only)} / 难例 {len(hard)}")
    print(f"  合计 {len(auto)+len(review_only)+len(hard)} = {len(rows)} ✓")
    print(f"pyridine_rich 软标记：{len(pyr_rich)} 个（复核组内 {len(pyr_rich_in_review)}，中性组内 {len(pyr_rich_in_auto)}）")
    if src_charged:
        print(f"⚠️ 源已带电荷需复核：{[r['mol_id'] for r in src_charged]}")
    print(f"⚠️ 规则中性组是'暂定中性态'，非 pKa 已证明。强制 Epik {len(mandatory_epik_ids)} 个；"
          f"pyridine_rich 共 {len(pyr_rich)} 个，其中新增条件性复核 {len(conditional_pyridine_ids)} 个。")
    print(f"\n输出:\n  {OUT_CSV}\n  {OUT_MD}")


if __name__ == "__main__":
    main()
