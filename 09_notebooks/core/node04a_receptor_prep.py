# -*- coding: utf-8 -*-
"""
V3_clean/core/node04a_receptor_prep.py  ·  v0.2.0
=================================================
受体制备（Node 04 对接前置）—— 把活性位点定义与对接盒做成【可审计、代码断言】的产物。

三受体（催化二联体随病毒不同）：
  6W63  SARS-CoV-2 Mpro : His41 / Cys145   共晶配体 X77 → 盒中心可自动
  3V3M  SARS-CoV-1 Mpro : His41 / Cys145
  4YLU  MERS-CoV  Mpro : His41 / Cys148   ← 更正手稿 Cys145

【科学保险 —— 断言，非注释】：
  1) 盒中心到催化二联体中点 ≤ MAX_CENTER_TO_DYAD_ANGSTROM（拒绝 ~29.82 Å 坏回退）。
  2) 盒中心【自动优先】：有共晶配体 → 用其质心；无 → 允许固定，但必须可追溯
     （fixed_from_verified_reference + 来源），禁止无源硬编码。
  3) v0.2.0 新增：清洗后受体【必须】不含非水杂原子（共晶配体/离子）——
     否则口袋被占，docking 静默出错（checkagain35 指出的最重要漏洞）。
  4) v0.2.0 新增：盒尺寸【必须】容纳候选池最大配体（数据证明，非"沿用 22.5"）。

【本模块不做（须在你 HPC 以真实结构完成）】：PDB→PDBQT 转换（工具相关）。
【质子化事实（v0.2.4 更正）】：清洗阶段【删除】全部 H；极性 H 由 prep 阶段 Meeko
  按残基模板以【标准态】加回 —— 无 pKa 预测、无滴定，故【不得】标注为 "pH 7.4"。
（工具相关，做成可配置命令并记录）。离线仅验证几何/清洁/尺寸断言逻辑。
"""

import json
import math
import hashlib
import subprocess
from datetime import datetime, timezone
from pathlib import Path

RECEPTOR_PREP_VERSION = "0.2.6"

# 盒中心到催化二联体中点最大允许距离（Å）。4YLU 自动盒中心经核验 ~4.11 Å；8 Å 给足余量，
# 远小于 29.82 Å 坏回退，故断言必然拒绝回退。
MAX_CENTER_TO_DYAD_ANGSTROM = 8.0
# 自动盒中心与"记录参考"固定值的交叉核验容差（Å）：超过即告警（可能一方有误）。
CROSS_CHECK_TOL_ANGSTROM = 6.0
# 盒边到最大配体长轴的每侧最小余量（Å）。低于此告警；盒边 < 最大长轴则硬失败。
MIN_BOX_MARGIN_PER_SIDE = 1.5

_CATALYTIC_ATOM = {"HIS": ("NE2", "CA"), "CYS": ("SG", "CA")}
_WATER = {"HOH", "WAT", "H2O", "DOD", "TIP", "TIP3", "SOL"}

# ─────────── 每个直系同源物的权威配置 ───────────
RECEPTORS = {
    "6W63": {
        "virus": "SARS-CoV-2", "chain": "A",
        "catalytic_dyad": [("HIS", 41), ("CYS", 145)],
        "cocrystal_ligand": "X77",                 # 有共晶配体 → 自动质心优先
        "box_center_reference": (-17.369, 18.129, -31.220),   # 交叉核验用
        "box_center_reference_source": "prior PRISM-Gen docking box (to verify)",
        "box_size": (25.0, 25.0, 25.0),
    },
    "3V3M": {
        "virus": "SARS-CoV-1", "chain": "A",
        "catalytic_dyad": [("HIS", 41), ("CYS", 145)],
        "cocrystal_ligand": None,                  # 若你的 3V3M 有共晶配体，填其 HET 以启用自动
        "box_center_reference": (21.845, -31.086, -4.129),
        "box_center_reference_source": "prior PRISM-Gen docking box (to verify)",
        "box_size": (25.0, 25.0, 25.0),
    },
    "4YLU": {
        "virus": "MERS-CoV", "chain": "A",
        "catalytic_dyad": [("HIS", 41), ("CYS", 148)],   # 更正 Cys148
        "cocrystal_ligand": None,                  # ← 必填：4YLU 共晶配体 HET code
        "box_center_reference": None,
        "box_center_reference_source": None,
        "box_size": (25.0, 25.0, 25.0),
    },
}


def utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _dist(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def parse_pdb(path):
    """返回 [(record, atom_name, res_name, chain, res_seq, (x,y,z), element)]。"""
    atoms = []
    for ln in Path(path).read_text().splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            rec = ln[0:6].strip(); name = ln[12:16].strip(); res = ln[17:20].strip()
            ch = ln[21:22].strip()
            try:
                seq = int(ln[22:26])
            except ValueError:
                continue
            x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
            elem = ln[76:78].strip() or name[:1]
            atoms.append((rec, name, res, ch, seq, (x, y, z), elem))
    return atoms


def dyad_midpoint(atoms, dyad, chain):
    """催化二联体中点（His 催化原子 与 Cys 催化原子中点；缺则回退 CA）。"""
    picked = []
    for resname, resseq in dyad:
        chosen = None
        for want in _CATALYTIC_ATOM.get(resname, ("CA", "CA")):
            for rec, name, res, ch, seq, xyz, _ in atoms:
                if rec == "ATOM" and ch == chain and res == resname \
                        and seq == resseq and name == want:
                    chosen = (want, xyz); break
            if chosen:
                break
        if chosen is None:
            raise ValueError(f"未找到 {resname}{resseq}（链 {chain}）催化原子或 CA "
                             f"—— 检查链号/残基号/结构完整性")
        picked.append((resname, resseq, chosen[0], chosen[1]))
    mid = tuple((picked[0][3][i] + picked[1][3][i]) / 2.0 for i in range(3))
    detail = [{"res": f"{p[0]}{p[1]}", "atom": p[2], "xyz": [round(c, 3) for c in p[3]]}
              for p in picked]
    return mid, detail


def ligand_centroid(atoms, het_code, chain=None):
    """共晶配体重原子质心（HETATM 且 resname==het_code）。"""
    pts = [xyz for rec, name, res, ch, seq, xyz, elem in atoms
           if rec == "HETATM" and res == het_code and elem != "H"
           and (chain is None or ch == chain)]
    if not pts:
        raise ValueError(f"未找到共晶配体 HET='{het_code}'"
                         f"{'' if chain is None else '（链 '+chain+'）'} —— 无法自动定盒中心，"
                         f"且禁止硬编码回退")
    c = tuple(sum(p[i] for p in pts) / len(pts) for i in range(3))
    return c, len(pts)


def resolve_box_center(cfg, atoms):
    """盒中心：自动优先（共晶配体质心），否则固定但必须可追溯。
    返回 (center, source, extra)。"""
    het = cfg.get("cocrystal_ligand")
    ref = cfg.get("box_center_reference")
    if het:
        c, n = ligand_centroid(atoms, het, cfg.get("chain"))
        extra = {"cocrystal_ligand": het, "n_ligand_heavy_atoms": n}
        if ref is not None:                        # 交叉核验：自动 vs 记录参考
            d = _dist(c, tuple(ref))
            extra["cross_check_vs_reference_angstrom"] = round(d, 3)
            extra["cross_check_reference"] = [round(x, 3) for x in ref]
            if d > CROSS_CHECK_TOL_ANGSTROM:
                extra["cross_check_warning"] = (
                    f"自动质心与记录参考相差 {d:.1f} Å (> {CROSS_CHECK_TOL_ANGSTROM}) —— 请核实")
        return c, "auto_from_cocrystal_ligand", extra
    # 无共晶配体 → 固定，但须可追溯
    if ref is None:
        raise ValueError("无共晶配体且无 box_center_reference —— 无法定盒中心（禁止硬编码）")
    if not cfg.get("box_center_reference_source"):
        raise ValueError("固定盒中心必须提供 box_center_reference_source（可追溯来源）")
    return tuple(ref), "fixed_from_verified_reference", \
        {"reference_source": cfg["box_center_reference_source"]}


def assert_box_sane(center, atoms, cfg):
    """盒中心必须落在催化二联体附近（拒绝 ~29.82 Å 回退）。"""
    mid, detail = dyad_midpoint(atoms, cfg["catalytic_dyad"], cfg["chain"])
    d = _dist(center, mid)
    if d > MAX_CENTER_TO_DYAD_ANGSTROM:
        raise AssertionError(
            f"盒中心距催化二联体中点 {d:.2f} Å > {MAX_CENTER_TO_DYAD_ANGSTROM} Å "
            f"—— 疑似硬编码回退/错误位点。center={tuple(round(c,3) for c in center)} "
            f"dyad_mid={tuple(round(c,3) for c in mid)}")
    return {"center_to_dyad_angstrom": round(d, 3),
            "dyad_midpoint": [round(c, 3) for c in mid],
            "dyad_atoms": detail,
            "assertion": f"<= {MAX_CENTER_TO_DYAD_ANGSTROM} A OK"}


def check_clean_receptor(clean_pdb, allowed_het=frozenset()):
    """清洗后受体的非水杂原子清单（应为空）。allowed_het 白名单（如结构性离子，需理由）。"""
    atoms = parse_pdb(clean_pdb)
    het = sorted({res for rec, name, res, ch, seq, xyz, elem in atoms
                  if rec == "HETATM" and res not in _WATER and res not in allowed_het})
    n_water = sum(1 for rec, name, res, *_ in atoms if rec == "HETATM" and res in _WATER)
    return {"remaining_non_water_het": het, "n_remaining_non_water_het": len(het),
            "n_water_atoms_present": n_water, "allowed_het": sorted(allowed_het)}


def assert_clean_receptor(clean_pdb, allowed_het=frozenset()):
    """★ 断言清洗后受体不含非水杂原子（共晶配体/离子未删 = 口袋被占 = 静默错误）。"""
    r = check_clean_receptor(clean_pdb, allowed_het)
    if r["n_remaining_non_water_het"] > 0:
        raise AssertionError(
            f"清洗后受体仍含非水 HETATM {r['remaining_non_water_het']} —— "
            f"共晶配体/杂原子未删，口袋被占，docking 会静默出错。"
            f"（如为必要辅因子，请显式加入 allowed_het 并说明理由）")
    return r


def assert_box_fits(box_size, max_ligand_axis, min_margin=MIN_BOX_MARGIN_PER_SIDE,
                    n_tight=None, n_over=None):
    """盒尺寸必须容纳候选池最大配体。盒边 < 最大长轴 → 硬失败；余量不足 → 记录告警。"""
    edge = min(box_size)
    margin = (edge - max_ligand_axis) / 2.0
    rec = {"box_min_edge_angstrom": edge,
           "max_ligand_longest_axis_angstrom": round(max_ligand_axis, 2),
           "margin_per_side_angstrom": round(margin, 2),
           "min_required_margin_angstrom": min_margin,
           "n_ligands_over_box": n_over, "n_ligands_tight": n_tight,
           "adequate": margin >= min_margin}
    if edge < max_ligand_axis:
        raise AssertionError(
            f"盒边 {edge} Å < 最大配体长轴 {max_ligand_axis:.1f} Å —— 最大分子放不进盒，"
            f"必须增大盒尺寸")
    if margin < min_margin:
        rec["warning"] = (f"最大配体每侧余量仅 {margin:.1f} Å (< {min_margin}) —— "
                          f"极端拉长分子采样受限；核实这些分子或增大盒")
    return rec


def read_pdbqt_dyad_state(pdbqt_path, chain, dyad):
    """v0.2.4 · 直接从产出的 PDBQT 读催化二联体状态（正向断言的依据）。
    不依赖任何 stdout 文本；按 PDB 列格式解析残基与原子名。
    His 互变异构体由 H 原子名判定：HD1→δ(HID)，HE2→ε(HIE)，两者皆有→HIP(带电)。
    Cys：有 HG→中性硫醇(SH)，无 HG→硫醇盐(S-)。"""
    seen = {}
    for ln in Path(pdbqt_path).read_text().splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        if ln[21] != chain:
            continue
        try:
            seq = int(ln[22:26])
        except ValueError:
            continue
        seen.setdefault(seq, {"resname": ln[17:20].strip(), "atoms": set()})
        seen[seq]["atoms"].add(ln[12:16].strip())
    out = {}
    for resname, seq in dyad:
        key = f"{resname}{seq}"
        rec = seen.get(seq)
        if rec is None:
            out[key] = {"present": False, "state": None, "atoms": []}
            continue
        A = rec["atoms"]
        if resname.upper().startswith("HIS"):
            d, e = ("HD1" in A), ("HE2" in A)
            state = ("HIP_charged" if (d and e) else
                     "HID_delta" if d else "HIE_epsilon" if e else "no_polar_H")
            heavy_ok = "NE2" in A and "ND1" in A
        elif resname.upper().startswith("CYS"):
            state = "SH_neutral_thiol" if "HG" in A else "S-_thiolate"
            heavy_ok = "SG" in A
        else:
            state, heavy_ok = "unknown_residue_type", True
        out[key] = {"present": True, "resname_in_pdbqt": rec["resname"],
                    "state": state, "key_heavy_atoms_ok": heavy_ok,
                    "atoms": sorted(A)}
    return out


def heavy_atom_census(clean_pdb, chain):
    """v0.2.6 · 清洗后 PDB 的重原子普查（清洗阶段已删 H，故 ATOM 行即重原子）。
    返回 (总重原子数, {"A:306": 该残基重原子数, ...})，用于转换前后的守恒断言。"""
    total, per_res = 0, {}
    for ln in Path(clean_pdb).read_text().splitlines():
        if not ln.startswith("ATOM"):
            continue
        if ln[21] != chain:
            continue
        total += 1
        key = f"{ln[21]}:{int(ln[22:26])}"
        per_res[key] = per_res.get(key, 0) + 1
    return total, per_res


def resolve_tool(tool_cmd_template):
    """v0.2.5 · 把命令模板的第一个 token（可执行文件）解析为【绝对路径】。
    动机：subprocess(shell=True) 走 /bin/sh，不会 source conda 激活脚本；
    "在登录 shell 里验证过" 不等于 "Jupyter kernel 的 PATH 里有"。
    解析顺序：PATH → 当前解释器的 bin 目录（同一 conda env）。
    失败时给出可操作的诊断，而不是让 /bin/sh 抛 command not found。"""
    import shutil, sys, os
    exe = tool_cmd_template.split()[0]
    if os.path.sep in exe:                       # 已经是路径
        cand = Path(exe)
        found = str(cand) if cand.is_file() else None
    else:
        found = shutil.which(exe)
        if not found:                            # 同 env 的 bin/ 兜底
            sibling = Path(sys.executable).parent / exe
            found = str(sibling) if sibling.is_file() else None
    if not found:
        raise AssertionError(
            f"受体制备工具 '{exe}' 不可执行。\n"
            f"  当前解释器: {sys.executable}\n"
            f"  该 env 的 bin: {Path(sys.executable).parent}\n"
            f"  CONDA_DEFAULT_ENV: {os.environ.get('CONDA_DEFAULT_ENV')}\n"
            f"  PATH 前 3 段: {os.environ.get('PATH','').split(os.pathsep)[:3]}\n"
            "  处理：确认 Meeko 装在【运行本 notebook 的这个 env】里（conda list meeko），\n"
            "        或把 TOOL_CMD_TEMPLATE 的第一个 token 换成绝对路径。\n"
            "  说明：conda activate 只影响登录 shell，不影响已启动的 Jupyter kernel。")
    ver = None
    try:
        import meeko as _mk
        ver = getattr(_mk, "__version__", None)
    except Exception:
        pass
    resolved = tool_cmd_template.replace(exe, found, 1)
    return resolved, {"tool_exe": exe, "tool_resolved_path": found,
                      "meeko_version": ver, "python_executable": sys.executable}


def receptor_to_pdbqt(clean_pdb, out_pdbqt, tool_cmd_template, dry_run=True, chain="A"):
    """PDB→受体 PDBQT。tool_cmd_template 可配置，如
      ADFR:  'prepare_receptor -r {in} -o {out}'
      Meeko: 'mk_prepare_receptor.py --read_pdb {in} -o {out} -p'
    绝不在此重新决定配体质子化。"""
    out_pdbqt = Path(out_pdbqt)
    out_base = out_pdbqt.with_suffix("")      # 去 .pdbqt —— Meeko 的 -o 需 basename
    resolved_template, tool_info = resolve_tool(tool_cmd_template)   # ★ v0.2.5 先解析再跑
    cmd = resolved_template.format(**{"in": str(clean_pdb), "out": str(out_pdbqt),
                                      "out_base": str(out_base)})
    if dry_run:
        return {"cmd": cmd, "executed": False, **tool_info}
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0 or not out_pdbqt.is_file():
        raise RuntimeError(f"受体 PDBQT 生成失败: {(r.stderr or r.stdout)[:300]}")
    # 转换后 QC：确认受体已质子化（加了极性 H）且无残留配体
    lines = [l for l in out_pdbqt.read_text().splitlines() if l.startswith("ATOM")]
    n_hd = sum(1 for l in lines if l.split()[-1] in ("HD", "H"))
    # 捕获 Meeko -a 丢弃的残基（缺原子残基）——记录以便审计
    import re as _re
    dropped = []
    mtc = _re.search(r"Template matching failed for: \[([^\]]*)\]",
                     (r.stdout or "") + (r.stderr or ""))
    if mtc:
        dropped = [x.strip().strip("'\"") for x in mtc.group(1).split(",") if x.strip()]
    # ★ v0.2.6 fail-closed：重原子守恒。转换只应【加氢】与【丢弃已声明的缺原子残基】，
    #   不应静默丢失任何其它重原子。此断言使 stdout 抓取的 fail-open 失去危害：
    #   若某残基被丢弃而正则没抓到，守恒式立刻不成立。
    n_heavy_pdbqt = len(lines) - n_hd
    total_clean, per_res = heavy_atom_census(clean_pdb, chain)
    n_heavy_dropped = sum(per_res.get(k, 0) for k in dropped)
    unknown_drop = [k for k in dropped if k not in per_res]
    expected = total_clean - n_heavy_dropped
    if n_heavy_pdbqt != expected:
        raise AssertionError(
            f"重原子不守恒：clean_pdb {total_clean} - 已声明丢弃 {n_heavy_dropped} "
            f"({dropped or '无'}) = 期望 {expected}，但 PDBQT 实得 {n_heavy_pdbqt} "
            f"（差 {n_heavy_pdbqt - expected}）。存在【未被记录】的原子丢失/新增，禁止进入对接。")

    return {"cmd": cmd, "executed": True, "pdbqt_sha256": file_sha256(out_pdbqt),
            "n_pdbqt_atoms": len(lines), "n_polar_h": n_hd, "protonated": n_hd > 0,
            "n_heavy_atoms_clean_pdb": total_clean,
            "n_heavy_atoms_pdbqt": n_heavy_pdbqt,
            "n_heavy_atoms_dropped": n_heavy_dropped,
            "heavy_atom_conservation": "OK",
            "dropped_keys_not_in_clean_pdb": unknown_drop,
            "dropped_incomplete_residues": dropped,
            "dropped_source": "stdout_scrape_AUDIT_ONLY",
            "protonation_method": "meeko_residue_template_standard_state_no_pKa",
            **tool_info}


def prepare_receptor(pdb_id, pdb_path, out_dir, receptors=RECEPTORS,
                     tool_cmd_template=None, clean_pdb_path=None,
                     allowed_het=frozenset(), pool_max_axis=None, dry_run=True):
    """单受体：解析 → 定盒中心 → 断言活性位点 → (有 clean_pdb)断言口袋已清空
    → (有 pool_max_axis)断言盒容纳最大配体 → (有工具)转换 → 记录。"""
    cfg = receptors[pdb_id]
    atoms = parse_pdb(pdb_path)                     # 原始结构（含共晶配体，用于自动定中心）
    center, src, extra = resolve_box_center(cfg, atoms)
    sane = assert_box_sane(center, atoms, cfg)

    clean_check = None
    if clean_pdb_path:                              # ★ 口袋清空断言（最重要新增）
        clean_check = assert_clean_receptor(clean_pdb_path, allowed_het)

    box_fit = None
    if pool_max_axis is not None:                   # ★ 盒尺寸数据核验
        box_fit = assert_box_fits(cfg["box_size"], pool_max_axis)

    out_dir = Path(out_dir)
    out_pdbqt = out_dir / f"{pdb_id}_receptor.pdbqt"
    conv = None
    if tool_cmd_template and clean_pdb_path:
        conv = receptor_to_pdbqt(clean_pdb_path, out_pdbqt, tool_cmd_template,
                                 dry_run, chain=cfg["chain"])
        # ★ 安全断言：Meeko -a 丢弃的缺原子残基【绝不能】包含催化二联体
        if conv and conv.get("dropped_incomplete_residues"):
            dyad_keys = {f"{cfg['chain']}:{s}" for _, s in cfg["catalytic_dyad"]}
            hit = dyad_keys & set(conv["dropped_incomplete_residues"])
            if hit:
                raise AssertionError(
                    f"催化二联体残基 {hit} 因缺原子被丢弃 —— 活性位点被破坏，必须先修复结构")
        # ★ v0.2.4 fail-closed：上面的 stdout 抓取会 fail-OPEN（Meeko 措辞一变即静默放行）。
        #   故改为对【产物本身】正向断言：二联体必须在 PDBQT 中、关键重原子齐全、且已加极性 H。
        if conv and conv.get("executed"):
            if not conv.get("protonated"):
                raise AssertionError(
                    f"{pdb_id}: PDBQT 中极性 H 计数为 0 —— 受体未质子化，禁止进入对接")
            st = read_pdbqt_dyad_state(out_pdbqt, cfg["chain"], cfg["catalytic_dyad"])
            conv["dyad_state"] = st
            for key, rec in st.items():
                if not rec["present"]:
                    raise AssertionError(
                        f"{pdb_id}: 催化二联体残基 {key} 不在产出的 PDBQT 中 —— 活性位点被破坏")
                if not rec.get("key_heavy_atoms_ok", True):
                    raise AssertionError(
                        f"{pdb_id}: {key} 关键重原子缺失（His 需 ND1/NE2，Cys 需 SG）")

    return {
        "pdb_id": pdb_id, "virus": cfg["virus"], "chain": cfg["chain"],
        "catalytic_dyad": [f"{n}{s}" for n, s in cfg["catalytic_dyad"]],
        "box_center": [round(c, 3) for c in center],
        "box_center_source": src,
        "box_center_extra": extra,
        "box_size": list(cfg["box_size"]),
        "box_sanity": sane,
        "box_fit": box_fit,
        "clean_receptor_check": clean_check,
        "input_pdb": str(pdb_path),
        "input_pdb_sha256": file_sha256(pdb_path) if Path(pdb_path).exists() else None,
        "clean_pdb": str(clean_pdb_path) if clean_pdb_path else None,
        "receptor_pdbqt": str(out_pdbqt),
        "conversion": conv,
        "notes": "清洗阶段删除全部 H；极性 H 由 Meeko residue template 以标准态在 prep 阶段加回"
                 "（无 pKa 预测/无滴定）；催化二联体=中性 His + 中性 Cys 硫醇；"
                 "本模块不重新决定配体质子化",
    }


def run_receptor_prep(pdb_paths, out_dir, receptors=RECEPTORS,
                      tool_cmd_template=None, clean_pdbs=None,
                      allowed_het=frozenset(), pool_max_axis=None, dry_run=True):
    """三受体批量制备 + 写 manifest。"""
    out_dir = Path(out_dir)
    clean_pdbs = clean_pdbs or {}
    mans = {}
    for pdb_id, pth in pdb_paths.items():
        mans[pdb_id] = prepare_receptor(
            pdb_id, pth, out_dir, receptors,
            tool_cmd_template=tool_cmd_template,
            clean_pdb_path=clean_pdbs.get(pdb_id),
            allowed_het=allowed_het, pool_max_axis=pool_max_axis, dry_run=dry_run)
    # ★ v0.2.4 fail-closed：三受体的二联体质子化态【必须一致】。
    #   逐受体不同的 His41 互变异构体会直接偏置 BBS（三靶取最大）分布。
    states = {}
    for pid, m in mans.items():
        c = (m or {}).get("conversion") or {}
        ds = c.get("dyad_state")
        if not ds:
            continue
        states[pid] = {("HIS" if k.upper().startswith("HIS") else "CYS"): v["state"]
                       for k, v in ds.items()}
    if len(states) > 1:
        uniq = {tuple(sorted(v.items())) for v in states.values()}
        if len(uniq) > 1:
            raise AssertionError(
                f"催化二联体质子化态跨受体不一致 —— BBS 跨靶点不可比: {states}")

    manifest = {
        "module": "node04a_receptor_prep", "version": RECEPTOR_PREP_VERSION,
        "created_at_utc": utc_now(), "dry_run": dry_run,
        "max_center_to_dyad_angstrom": MAX_CENTER_TO_DYAD_ANGSTROM,
        "pool_max_ligand_axis_angstrom": pool_max_axis,
        "dyad_protonation_uniform": (states if states else None),
        "policy": {
            "catalytic_dyad": "per-ortholog; MERS(4YLU)=His41/Cys148 (corrected from Cys145)",
            "box_center": "auto from co-crystal ligand centroid where present; "
                          "else fixed_from_verified_reference (traceable); no bare hardcode",
            "box_center_sanity": f"asserted <= {MAX_CENTER_TO_DYAD_ANGSTROM} A from dyad midpoint "
                                 "(rejects ~29.82 A fallback)",
            "clean_receptor": "asserted NO non-water HETATM in cleaned receptor "
                              "(co-crystal ligand must be removed; else occupied-pocket silent error)",
            "box_size": "asserted to accommodate the largest candidate (data-driven, not assumed)",
            "protonation": "cleaning stage STRIPS all H; polar H re-added by Meeko "
                           "residue templates at STANDARD STATE during prep "
                           "(no pKa prediction, no titration -- NOT a pH 7.4 calculation); "
                           "catalytic dyad = neutral His + neutral Cys thiol; "
                           "ligand protonation NEVER re-decided",
            "dyad_integrity": "asserted POSITIVELY on the output PDBQT (dyad residues + key "
                              "atoms present); stdout scraping kept only as audit record",
            "dyad_tautomer": "recorded from PDBQT H-atom names and asserted IDENTICAL across "
                             "all receptors (cross-target comparability of BBS)",
            "tool_resolution": "converter resolved to an ABSOLUTE path before execution and "
                               "recorded in manifest (shell=True does not honour conda activate)",
            "heavy_atom_conservation": "asserted heavy(PDBQT) == heavy(clean PDB) - heavy(declared "
                                       "dropped residues); makes the stdout scrape non-load-bearing",
        },
        "receptors": mans,
    }
    # ★ v0.2.4 fail-closed：正式运行（dry_run=False）时【禁止】写出一份转换未完成的 manifest。
    #   历史事故：TOOL_CMD_TEMPLATE 为空 → conversion=null，但 manifest 照写照冻结，
    #   随后被当作 Node 04 的锚点。manifest 必须只在完整时才存在。
    if not dry_run:
        incomplete = [pid for pid, m in mans.items()
                      if not ((m.get("conversion") or {}).get("executed"))]
        if incomplete:
            raise AssertionError(
                f"拒绝写出不完整 manifest：{incomplete} 未完成 PDB→PDBQT 转换 "
                f"(tool_cmd_template={'空' if not tool_cmd_template else '已设'})。"
                "请设置 TOOL_CMD_TEMPLATE 后重跑；不完整的 manifest 绝不可作为 Node 04 锚点。")

    if not dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "Node04a_receptor_manifest.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return manifest


# ─────────── v0.2.1: 统一受体清洗（三受体同一流程，保证可比） ───────────
def clean_receptor_pdb(raw_pdb, out_pdb, chain="A", keep_altloc=(" ", "A"),
                       keep_het=frozenset(), strip_hydrogens=True):
    """统一清洗：只留指定链的蛋白 ATOM；去共晶配体/溶剂/离子/全部水；
    多构象只留 altloc ' ' 或 'A'（并清空 altloc 指示符）。结晶结构通常无 H，
    全部 H 在此删除；极性 H 由 prep 阶段 Meeko 以标准态模板加回（非 pH 7.4 计算）。
    keep_het 白名单默认空（全去）。
    返回清洗审计（保留原子数、丢弃分类）。"""
    import collections
    kept, dropped = [], collections.Counter()
    last_atom_line = None
    for ln in Path(raw_pdb).read_text().splitlines():
        if ln.startswith("ATOM"):
            ch = ln[21]; alt = ln[16]
            if ch != chain:
                dropped[f"chain_{ch}"] += 1; continue
            if alt not in keep_altloc:
                dropped["altloc_alt"] += 1; continue
            elem = (ln[76:78].strip() or ln[12:16].strip()[:1])
            if strip_hydrogens and elem == "H":   # 统一去 H；极性 H 由 prep 阶段 Meeko 标准态模板加回
                dropped["hydrogen"] += 1; continue
            ln = ln[:16] + " " + ln[17:]          # 清空 altloc 指示符
            kept.append(ln); last_atom_line = ln
        elif ln.startswith("HETATM"):
            res = ln[17:20].strip(); ch = ln[21]
            if res in keep_het and ch == chain:
                ln = ln[:16] + " " + ln[17:]
                kept.append(ln); continue
            dropped[("water" if res in _WATER else f"HET_{res}")] += 1
    lines = list(kept)
    if last_atom_line is not None:                # 收尾 TER
        lines.append("TER")
    lines.append("END")
    Path(out_pdb).write_text("\n".join(lines) + "\n")
    return {"raw_pdb": str(raw_pdb), "clean_pdb": str(out_pdb),
            "chain": chain, "kept_altloc": list(keep_altloc),
            "keep_het": sorted(keep_het), "stripped_hydrogens": strip_hydrogens,
            "n_atoms_kept": len(kept), "dropped": dict(dropped)}
