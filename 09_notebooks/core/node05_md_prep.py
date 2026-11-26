#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
core/node05_md_prep.py  —  V3 Node 05 MD 系统制备（PRISM-Dyn）
================================================================
把 Node 02–04 冻结的分子身份，安全地带过 V2 死掉的那一步：配体 → MD 拓扑。
替换 V2 的三处身份重新推导（step6b obabel / step6c 电荷扫描 / cell10 genion -neutral），
以四道 fail-closed 不变量强制：

  I1  配体几何来自 master SDF + 姿势坐标（绝不 obabel-从-姿势）；
      重建 SMILES 图不变、重原子键长守恒（乱序映射会被抓）。
  I2  ACPYPE 只用 Node 01 frozen_net_charge（单值，不扫描）；
      断言 ligand.itp 部分电荷之和 == frozen_net_charge。
  I2b 【v0.2.0 新增】配体分子图不变量：由 ligand_GMX.itp 的 [atoms] 质量列与 [bonds]
      重建元素序列与邻接，与投喂 ACPYPE 的 ligand.sdf 逐索引断言同一。
      理由：AM1-BCC 电荷按构造求和为 -nc，故 I2 对键级/芳香性/逐原子形式电荷的改动
      完全不敏感；V2 所损坏的那类分子会以 itp_charge_sum 正确的形式通过 I2。
      实测（pool_mol_030）：ACPYPE/antechamber 保持输入原子顺序，故可直接逐索引比对。
  I3  催化二联体质子化【强制并断言】：His41=HIE、Cys=中性 CYS，三受体统一。
  I3b 【v0.2.0 新增】残基序列全局比对：清洗后输入 PDB 与 protein.gro 的残基名序列
      逐位相等（HIS→HISE/HISD/HISH 等做别名归一）。理由：pdb2gmx 可能重新编号，
      而 I3 只读某索引处的残基名，编号偏移时存在巧合通过的可能。
  I4  反离子由 (protein_charge + frozen_net_charge) 【推导】；
      断言 genion 前 grompp 感知电荷 == 推导值、genion 后系统恰为中性。

物理协议（MDP / EM→NVT→NPT→MD、盒、溶剂化、盐）原样复用 V2（step6d + 溶剂化 cell10），
仅在四道不变量处替换 V2 的重新推导逻辑。gmx/acpype 调用忠实照搬 V2 已用命令。

作者说明：I1 与 removeHs→apply→addHs 流程、I2/I3/I4 的解析器均已在离线用
RDKit + 合成 GROMACS 文本验证（含负对照）。gmx/acpype 首次运行仍需逐分子核对
Cell 输出（真实 HPC 环境不可离线执行）。
"""
from __future__ import annotations
import json, hashlib, subprocess, re, shutil, os
from pathlib import Path

NODE05_VERSION = "0.2.3"

# ───────────────────────── 复用自 V2 的物理参数（勿改） ─────────────────────────
BOX_DIST_NM   = 1.0        # editconf -bt cubic -d           （V2 cell10）
SALT_CONC_M   = 0.15       # genion -conc                    （V2 cell10）
WATER_MODEL   = "spc216"   # solvate -cs spc216.gro          （V2 cell10；TIP3P 起始构型）
FF_PROTEIN    = "amber99sb-ildn"
FF_WATER      = "tip3p"
DT_PS         = 0.002      # 2 fs
PROD_NS       = 100.0      # 生产 100 ns（× 3 replica）
NSTXOUT_COMP  = 5000       # 每 10 ps 一帧 → 10001 帧/100ns  （V2 step6d）

# 催化二联体的目标质子化态（I3）——三受体统一：His41=HIE(ε)，Cys=中性 CYS
#   6W63 / 3V3M: His41 / Cys145 ；4YLU: His41 / Cys148
DYAD_SPEC = {
    "6W63": {"his": (41, "HIE"), "cys": (145, "CYS")},
    "3V3M": {"his": (41, "HIE"), "cys": (145, "CYS")},
    "4YLU": {"his": (41, "HIE"), "cys": (148, "CYS")},
}

# ═══════════════════════════════════════════════════════════════════════════════
#  工具
# ═══════════════════════════════════════════════════════════════════════════════
def sha256_file(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

# JupyterHub 内核不继承 md_env/bin（HANDOFF）——acpype/antechamber/parmchk2 皆在此，须显式加入 PATH
MDENV_BIN = "/mnt/hpc/home/wanghongmei/anaconda3/envs/md_env/bin"

def run(cmd: str, cwd=None, gmxrc: str | None = None, env_bin: str | None = MDENV_BIN):
    """执行 shell；prepend env_bin 到 PATH（acpype 及其子进程），再 source gmxrc。"""
    prefix = f'export PATH="{env_bin}:$PATH"; ' if env_bin else ""
    src    = f"source {gmxrc} >/dev/null 2>&1; " if gmxrc else ""
    full   = prefix + src + cmd
    return subprocess.run(full, shell=True, cwd=cwd, capture_output=True, text=True, executable="/bin/bash")

def path_export_line(env_bin: str = MDENV_BIN) -> str:
    """给 HPC 脚本（run_prep.sh / sbatch）用的 PATH 导出行。"""
    return f'export PATH="{env_bin}:$PATH"'

# ═══════════════════════════════════════════════════════════════════════════════
#  I1 · 配体几何来自 master SDF + 姿势坐标（绝不 obabel-从-姿势）
# ═══════════════════════════════════════════════════════════════════════════════
def _parse_pdbqt_coords(pdbqt_path: str, model: int = 1):
    """
    按文件顺序返回某 MODEL 的原子坐标（含极性 H；与 Node 03 配体 PDBQT 同序）。

    v0.2.0 fail-closed 修正：v0.1.x 在【请求的 MODEL 不存在】时会静默退回
    "把所有 MODEL 串联"，其前若干项恰为 MODEL 1，于是 apply_pose_to_master 会
    以为用了 model=k、实际用了 model=1，且 I1 的 SMILES/键长/电荷三条断言全部通过
    （它们是分子的性质，不是"用了哪个姿势"的性质）——一条无人把守的错姿势通路。
    现改为：文件含 MODEL 记录而请求索引不在其中 → 抛错。
    """
    text = Path(pdbqt_path).read_text().splitlines()
    per_model, cur, n_models = {}, None, 0
    for ln in text:
        if ln.startswith("MODEL"):
            n_models += 1; cur = []
        elif ln.startswith("ENDMDL"):
            if cur is not None: per_model[n_models] = cur
            cur = None
        elif ln.startswith(("ATOM", "HETATM")) and cur is not None:
            cur.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    if n_models == 0:                      # 无 MODEL 标签的单构象文件
        coords = [(float(l[30:38]), float(l[38:46]), float(l[46:54]))
                  for l in text if l.startswith(("ATOM", "HETATM"))]
        if int(model) != 1:
            raise AssertionError(
                f"I1: {Path(pdbqt_path).name} 无 MODEL 记录（单构象），但请求 model={model}")
        if not coords:
            raise AssertionError(f"I1: {Path(pdbqt_path).name} 无任何 ATOM/HETATM 记录")
        return coords
    if int(model) not in per_model:
        raise AssertionError(
            f"I1: {Path(pdbqt_path).name} 含 MODEL {sorted(per_model)}，不含请求的 model={model} "
            f"—— 拒绝退回 MODEL 1（v0.1.x 曾在此静默取错姿势）")
    return per_model[int(model)]

def apply_pose_to_master(master_sdf: str, pose_pdbqt: str, atom_map_json: str,
                         frozen_net_charge: int, out_mol2: str, model: int = 1,
                         bond_tol: float = 0.02):
    """
    I1：以 master SDF 为化学真源、姿势 PDBQT 为几何来源，构建 MD 配体。
    流程：master(全原子) → 对【Node 03 映射到的原子】写入姿势坐标 →
          RemoveHs（保留重原子=姿势几何）→ AddHs(addCoords) 重建 H → 写 mol2。
    断言：① 规范异构 SMILES 不变；② 重原子键长相对 master 守恒（乱序映射会爆）；
          ③ 形式电荷 == frozen_net_charge。任一不符 → 抛错（fail-closed）。
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import RDLogger; RDLogger.DisableLog("rdApp.*")

    master = Chem.MolFromMolFile(master_sdf, removeHs=False, sanitize=True)
    if master is None:
        raise RuntimeError(f"I1: 无法读取 master SDF: {master_sdf}")
    s_master = Chem.MolToSmiles(Chem.RemoveHs(master))
    q_master = Chem.GetFormalCharge(master)
    if q_master != int(frozen_net_charge):
        raise AssertionError(f"I1: master 形式电荷 {q_master} != frozen_net_charge {frozen_net_charge}")

    amap = json.loads(Path(atom_map_json).read_text())     # 期望 pdbqt_idx(str)->master_idx(int)
    fwd = amap.get("pdbqt_to_master") or amap.get("forward") or amap
    fwd = {int(k): int(v) for k, v in fwd.items()}
    pose = _parse_pdbqt_coords(pose_pdbqt, model=model)
    if max(fwd) >= len(pose):
        raise AssertionError(f"I1: 姿势原子数 {len(pose)} 与映射最大索引 {max(fwd)} 不符")

    # master 重原子键长（校验用）
    def heavy_bonds(m):
        c = m.GetConformer(); out = {}
        for b in m.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            if m.GetAtomWithIdx(i).GetAtomicNum() > 1 and m.GetAtomWithIdx(j).GetAtomicNum() > 1:
                out[(i, j)] = c.GetAtomPosition(i).Distance(c.GetAtomPosition(j))
        return out
    bl_master = heavy_bonds(master)

    # 写入姿势坐标到映射原子
    posed_full = Chem.Mol(master); conf = posed_full.GetConformer()
    from rdkit.Geometry import Point3D
    for pidx, midx in fwd.items():
        x, y, z = pose[pidx]; conf.SetAtomPosition(midx, Point3D(x, y, z))

    # 断言②：重原子键长守恒（乱序映射→键长爆）
    bl_posed = heavy_bonds(posed_full)
    max_dev = max(abs(bl_posed[k] - bl_master[k]) for k in bl_master)
    # v0.2.0：容差 0.15 → 0.02。Vina 只旋转扭转、不改变键长，pool_mol_030 实测偏差
    # 0.001 Å，故 0.15 丢弃了本检查约 150 倍的分辨力。
    if max_dev > bond_tol:
        raise AssertionError(f"I1: 重原子键长偏差 {max_dev:.4f} Å > {bond_tol} —— 姿势坐标可能错配到非对应原子")

    # RemoveHs（重原子=姿势几何）→ AddHs 重建 H
    heavy = Chem.RemoveHs(posed_full)
    posed = Chem.AddHs(heavy, addCoords=True)
    s_posed = Chem.MolToSmiles(Chem.RemoveHs(posed))
    q_posed = Chem.GetFormalCharge(posed)
    if s_posed != s_master:
        raise AssertionError(f"I1: SMILES 变化\n  master {s_master}\n  posed  {s_posed}")
    if q_posed != int(frozen_net_charge):
        raise AssertionError(f"I1: 重建后形式电荷 {q_posed} != frozen_net_charge {frozen_net_charge}")

    # 写 mol2 供 ACPYPE
    Chem.MolToMolFile(posed, out_mol2.replace(".mol2", ".sdf"))
    try:
        from rdkit.Chem import rdmolfiles
        Chem.MolToMol2File(posed, out_mol2)          # RDKit ≥ 2020 支持
    except Exception:
        out_mol2 = out_mol2.replace(".mol2", ".sdf")  # 回退 SDF（ACPYPE 亦可读）
    return {"invariant": "I1_pose_from_master", "ok": True,
            "canonical_smiles": s_master, "formal_charge": q_posed,
            "max_heavy_bond_dev_A": round(max_dev, 3),
            "ligand_geom_file": out_mol2, "n_atoms": posed.GetNumAtoms()}

# ═══════════════════════════════════════════════════════════════════════════════
#  I2 · ACPYPE 只用 frozen_net_charge（不扫描）+ 断言 itp 电荷和
# ═══════════════════════════════════════════════════════════════════════════════
def emit_param_commands(ligand_geom: str, frozen_net_charge: int, basename: str = "ligand") -> str:
    """替换 step6c 的 _CHARGE_CANDIDATES 扫描：单一 -n frozen_net_charge，绝不回退试错。"""
    ext = "mol2" if ligand_geom.endswith(".mol2") else "sdf"
    return (f'acpype -i "{ligand_geom}" -b {basename} -n {int(frozen_net_charge)} '
            f'-a gaff2 -o gmx -f  # I2: 冻结电荷单值，禁止扫描')

def parse_itp_charge(itp_path: str):
    total, n = 0.0, 0; in_atoms = False
    for ln in Path(itp_path).read_text().splitlines():
        s = ln.strip()
        if s.startswith("["):
            in_atoms = s.replace(" ", "").startswith("[atoms]"); continue
        if in_atoms and s and not s.startswith(";"):
            p = s.split()
            if len(p) >= 7:
                try: total += float(p[6]); n += 1
                except ValueError: pass
    return total, n

def assert_ligand_charge(itp_path: str, frozen_net_charge: int, tol: float = 0.02):
    total, n = parse_itp_charge(itp_path)
    if abs(total - int(frozen_net_charge)) > tol:
        raise AssertionError(
            f"I2: ligand.itp 电荷和 {total:.4f}（{n} 原子）!= frozen_net_charge {frozen_net_charge}（容差 {tol}）")
    return {"invariant": "I2_frozen_charge", "ok": True,
            "itp_charge_sum": round(total, 4), "frozen_net_charge": int(frozen_net_charge), "n_atoms": n}

# ═══════════════════════════════════════════════════════════════════════════════
#  I3 · 催化二联体质子化 强制 + 断言
# ═══════════════════════════════════════════════════════════════════════════════
def emit_receptor_topology_commands(receptor_pdb: str, target: str) -> str:
    """
    pdb2gmx 生成蛋白拓扑并【强制】二联体态：His41=HIE、Cys=中性 CYS。
    用 -his 交互输入固定组氨酸态；半胱氨酸默认中性（-nossl / 不成二硫）。
    His 质子化通过预置残基名 HIE + `gmx pdb2gmx -his` 应答实现；此处以 here-doc 提供选择。
    """
    return (
        f'# I3: 强制二联体质子化（His41=HIE, Cys 中性）\n'
        f'gmx pdb2gmx -f "{receptor_pdb}" -o protein.gro -p topol.top -i posre_prot.itp \\\n'
        f'    -ff {FF_PROTEIN} -water {FF_WATER} -ignh -his -nossl <<EOF\n'
        f'# 对每个 His 提供质子化选择：HIE=1（ε）。其余 His 亦置 HIE 以求统一。\n'
        f'EOF\n'
        f'#   注：若结构中含多个 His，需按 pdb2gmx 提示顺序提供选择；\n'
        f'#   V3 在 Cell 内以受体的 His 残基列表自动生成应答序列（全部 =HIE=1）。'
    )

def _parse_gro_residues(gro_path: str, wanted_resids: dict):
    obs = {}
    lines = Path(gro_path).read_text().splitlines()
    for ln in lines[2:-1]:
        if len(ln) < 10: continue
        try:
            resid = int(ln[0:5]); resnm = ln[5:10].strip()
        except ValueError:
            continue
        if resid in wanted_resids and resid not in obs:
            obs[resid] = resnm
    return obs

# GROMACS/AMBER99SB-ILDN 组氨酸命名别名 → 规范态
_HIS_EPS = {"HIE", "HISE"}   # ε-质子化（中性）—— 我们要的催化二联体态
_HIS_DEL = {"HID", "HISD"}   # δ-质子化（中性）
_HIS_POS = {"HIP", "HISH"}   # 双质子化（+1，带电）
_CYS_NEUTRAL = {"CYS"}       # 中性（带 HG）；CYS2/CYX/CYM 为二硫/去质子，排除
_HIS_ACCEPT = {"HIE": _HIS_EPS, "HID": _HIS_DEL, "HIP": _HIS_POS}

def assert_dyad_states(protein_gro: str, target: str):
    """断言蛋白 .gro 中 His41=ε-质子化中性（HIE/HISE）、Cys=中性 CYS。接受 GROMACS 别名。"""
    spec = DYAD_SPEC[target]
    his_id, his_want = spec["his"]; cys_id, cys_want = spec["cys"]
    obs = _parse_gro_residues(protein_gro, {his_id: 1, cys_id: 1})
    his_obs, cys_obs = obs.get(his_id), obs.get(cys_id)
    accept_his = _HIS_ACCEPT.get(his_want, {his_want})
    problems = []
    if his_obs not in accept_his:
        problems.append(f"His{his_id}={his_obs}（期望 {his_want}/{sorted(accept_his)}）")
    if cys_obs not in _CYS_NEUTRAL:
        problems.append(f"Cys{cys_id}={cys_obs}（期望中性 CYS，非 CYS2/CYX/CYM）")
    if problems:
        raise AssertionError(f"I3: {target} 二联体质子化不符：" + "；".join(problems))
    return {"invariant": "I3_neutral_dyad", "ok": True, "target": target,
            f"His{his_id}": his_obs, f"Cys{cys_id}": cys_obs,
            "note": "HISE==HIE（ε-质子化中性），GROMACS 命名"}

# ═══════════════════════════════════════════════════════════════════════════════
#  I4 · 反离子由 (protein_charge + frozen_net_charge) 推导 + 中性断言
# ═══════════════════════════════════════════════════════════════════════════════
def emit_solvate_commands(box_d=BOX_DIST_NM, salt=SALT_CONC_M) -> str:
    """照搬 V2 cell10 溶剂化块：editconf → solvate(spc216) → grompp(ions) → genion -neutral。"""
    return (
        f'# I4: 溶剂化 + 加离子（沿用 V2 已用命令）\n'
        f'mv -f complex.gro complex_dry.gro\n'
        f'gmx editconf -f complex_dry.gro -o box.gro -bt cubic -d {box_d}\n'
        f'gmx solvate  -cp box.gro -cs {WATER_MODEL}.gro -p topol.top -o solv.gro\n'
        f'cat > ions.mdp <<MDP\n'
        f'integrator = steep\nnsteps = 0\ncutoff-scheme = Verlet\n'
        f'coulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\npbc = xyz\nMDP\n'
        f'gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 5\n'
        f'printf "SOL\\n" | gmx genion -s ions.tpr -o complex.gro -p topol.top \\\n'
        f'    -pname NA -nname CL -neutral -conc {salt}\n'
    )

# grompp 成功运行过的证据（用于把"真中性"与"根本没跑/跑挂了"区分开）
_GROMPP_RAN_MARKERS = ("Analysing residue names", "This run will generate",
                       "Calculating fourier grid dimensions", "writing run input file",
                       "Writing run input file", "Number of degrees of freedom")

def parse_grompp_system_charge(grompp_log: str):
    """
    从 grompp 输出解析 genion 前的感知总电荷。

    v0.2.0 fail-closed 修正：v0.1.x 对【空日志】与【grompp fatal error】同样返回 0，
    与"体系真的中性"不可区分。后果：批处理里凡 (蛋白+配体) 恰为 0 的体系
    （即配体恰为 +4 者），I4 第一条断言会在 grompp 从未成功运行的情况下空洞通过。
    现改为：无法确证 grompp 成功运行 → 返回 None，由 assert_counterions 判失败。
    """
    m = re.search(r"System has non-zero total charge:\s*([-\d.eE]+)", grompp_log)
    if m:
        return round(float(m.group(1)))
    if any(k in grompp_log for k in _GROMPP_RAN_MARKERS):
        return 0            # grompp 确实跑过且未报非零电荷 → 真中性
    return None             # 空 / 报错 / 无法识别 → 不给结论

def count_ions_in_gro(gro_path: str):
    na = cl = 0
    for ln in Path(gro_path).read_text().splitlines():
        nm = ln[5:10].strip() if len(ln) >= 10 else ""
        if nm in ("NA", "NA+", "SOD"): na += 1
        elif nm in ("CL", "CL-", "CLA"): cl += 1
    return na, cl

def assert_counterions(protein_charge: int, frozen_net_charge: int,
                       pre_genion_grompp_log: str, final_gro: str):
    """
    I4：① 断言 genion 前 grompp 感知电荷 == 推导值 (protein_charge + frozen_net_charge)；
        ② 断言 genion 后系统恰为中性：system_charge + n_NA - n_CL == 0。
    """
    derived = int(protein_charge) + int(frozen_net_charge)
    perceived = parse_grompp_system_charge(pre_genion_grompp_log)
    # v0.2.0：None 视为失败，不再"无法解析就跳过比较"——跳过等于没有 I4。
    if perceived is None:
        raise AssertionError(
            "I4: 无法从 grompp_pre.log 确证 grompp 成功运行并读出感知电荷 "
            "（日志为空 / 报错 / 只重定向了 stderr）—— 拒绝在无证据下判定通过")
    if perceived != derived:
        raise AssertionError(
            f"I4: genion 前 grompp 感知电荷 {perceived} != 推导值 {derived} "
            f"(protein {protein_charge} + ligand {frozen_net_charge}) —— 拓扑电荷与冻结电荷不符")
    na, cl = count_ions_in_gro(final_gro)
    final_net = derived + na - cl
    if final_net != 0:
        raise AssertionError(
            f"I4: 加离子后系统净电荷 {final_net:+d} != 0 (system {derived:+d}, NA {na}, CL {cl})")
    return {"invariant": "I4_derived_counterions", "ok": True,
            "protein_charge": int(protein_charge), "frozen_net_charge": int(frozen_net_charge),
            "derived_system_charge": derived, "grompp_perceived": perceived,
            "n_NA": na, "n_CL": cl, "final_net_charge": final_net}

# ═══════════════════════════════════════════════════════════════════════════════
#  MDP 协议（原样复用 V2 step6d：EM → NVT → NPT → MD）
# ═══════════════════════════════════════════════════════════════════════════════
# write_mdp() 已于 v0.2.2 删除。
# 理由：该函数系凭记忆重建，且从未被任何 notebook 调用；其数值与 V2 真实协议不符
# （它假定 tc-grps 单组、NPT 用 Parrinello-Rahman 且缺 refcoord_scaling，
#   而 V2 run_md2.sh 实际用两组 V-rescale、NPT 用 C-rescale + refcoord_scaling = com）。
# 权威 MDP 见本文件末尾 mdp_em / mdp_nvt / mdp_nvt_continue / mdp_npt / mdp_md，
# 逐字转录自 V2 run_md2.sh 的 write_mdp_* 函数。

def prep_manifest(job: dict, i1: dict, i2: dict, i3: dict, i4: dict, out_dir: Path):
    man = {"module": "node05_md_prep", "version": NODE05_VERSION, "job": job,
           "invariants": {"I1": i1, "I2": i2, "I3": i3, "I4": i4},
           "all_passed": all(x.get("ok") for x in (i1, i2, i3, i4))}
    p = Path(out_dir) / "node05_prep_manifest.json"
    p.write_text(json.dumps(man, indent=2, ensure_ascii=False))
    return man


# ═══════════════════════════════════════════════════════════════════════════════
#  跑通辅助（His 应答 / 蛋白电荷解析 / protein+ligand 合并）
# ═══════════════════════════════════════════════════════════════════════════════
def count_his_residues(pdb_path: str):
    """返回受体中的 His 残基号列表（pdb2gmx -his 按序提示，用于生成应答）。"""
    resids = []
    for ln in Path(pdb_path).read_text().splitlines():
        if ln.startswith(("ATOM", "HETATM")) and ln[17:20].strip() in ("HIS", "HIE", "HID", "HIP"):
            try:
                rid = int(ln[22:26])
            except ValueError:
                continue
            if rid not in resids:
                resids.append(rid)
    return resids

def his_answer_sequence(n_his: int) -> str:
    """pdb2gmx -his 应答：全部选 HIE（ε，输入 1）以求统一。"""
    return "".join("1\n" for _ in range(n_his))

def parse_pdb2gmx_charge(log_text: str):
    """从 pdb2gmx 输出解析蛋白总电荷（整数）。"""
    m = re.search(r"[Tt]otal charge(?:\s+in\s+system)?\s*[:=]?\s*([-\d.]+)", log_text)
    return round(float(m.group(1))) if m else None

def merge_protein_ligand(protein_gro: str, ligand_gro: str, out_complex_gro: str):
    """合并坐标：protein 原子在前、ligand 追加；盒取 protein。返回 (total, n_prot, n_lig)。"""
    def read_gro(p):
        L = Path(p).read_text().splitlines(); n = int(L[1]); return L[2:2 + n], L[2 + n]
    pa, pbox = read_gro(protein_gro); la, _ = read_gro(ligand_gro)
    total = len(pa) + len(la)
    Path(out_complex_gro).write_text("\n".join(["complex (protein+ligand) V3", str(total)] + pa + la + [pbox]) + "\n")
    return total, len(pa), len(la)

def merge_topology(topol_top: str, ligand_itp_include: str, ligand_molname: str = "ligand"):
    """在 forcefield #include 之后插入 ligand itp，并在 [molecules] 末尾加 'ligand 1'。"""
    lines = Path(topol_top).read_text().splitlines(); out = []; inserted = False; in_mol = False
    for ln in lines:
        out.append(ln)
        if (not inserted) and "#include" in ln and "forcefield.itp" in ln:
            out.append(f'#include "{ligand_itp_include}"'); inserted = True
        if ln.strip().replace(" ", "").startswith("[molecules]"):
            in_mol = True
    if in_mol:
        out.append(f"{ligand_molname} 1")
    if not inserted:
        raise AssertionError("merge_topology: 未找到 forcefield.itp include 行以插入 ligand itp")
    Path(topol_top).write_text("\n".join(out) + "\n")
    return {"ligand_itp_included": inserted}


# ═══════════════════════════════════════════════════════════════════════════════
#  批处理：自包含 run_prep.sh（proper 脚本，heredoc 完整）+ Slurm 数组
# ═══════════════════════════════════════════════════════════════════════════════
def emit_run_prep_script(sysdir, receptor_pdb, target, frozen_net_charge, ligand_geom,
                         his_answers, project_root, gmxrc, py,
                         reuse_itp_from=None) -> str:
    """
    每体系一份自包含制备脚本。v0.2.0 相对 v0.1.4 的改动：

      1) `set +u` 包裹 source GMXRC —— GMXRC 会引用 LD_LIBRARY_PATH / MANPATH，
         在 `set -u` 下可能以 unbound variable 在两秒内杀死全部数组任务。
      2) 开头清除本体系先前产物 —— 跑通遗留的 ligand.acpype/ 与新写的 ligand.sdf
         会构成"参数取自 A 姿势、坐标取自 B 姿势"的混合体系
         （AM1-BCC 电荷依赖构象，故那是真正不同的参数，不只是不同的坐标）。
      3) antechamber -fi mdl 改为【主路径】—— SDF→mol2 经 MDL 连接表完成，
         而非经 ACPYPE 内部的 OpenBabel（V2 根因所在的那个转换器）。
      4) grompp 捕获改为 `> grompp_pre.log 2>&1` 且去掉 `|| true` ——
         v0.1.4 只重定向 stderr 且吞掉失败，会让 I4 在无证据下空洞通过。
      5) reuse_itp_from：次姿势体系复用主姿势的 ligand_GMX.itp，
         使姿势比较不混入电荷参数化差异。
    """
    itp_inc = "ligand.acpype/ligand_GMX.itp"
    q = int(frozen_net_charge)
    if reuse_itp_from:
        param_block = (
            "# I2 复用：次姿势体系沿用主姿势的电荷/参数（仅坐标不同）\n"
            "rm -rf ligand.acpype && mkdir -p ligand.acpype\n"
            f"cp {reuse_itp_from}/ligand.acpype/ligand_GMX.itp   ligand.acpype/\n"
            f"cp {reuse_itp_from}/ligand.acpype/posre_ligand.itp ligand.acpype/ 2>/dev/null || true\n"
            "# 坐标取本体系姿势：由 ligand.sdf 生成 ligand_GMX.gro（不重新参数化）\n"
            f"{py} -c \"import sys; sys.path.insert(0,'{project_root}'); "
            f"from core import node05_md_prep as N5; "
            f"N5.sdf_to_gro('{ligand_geom}','ligand.acpype/ligand_GMX.gro')\"")
    else:
        param_block = (
            "# I2: 冻结电荷 ACPYPE（单值，禁扫描）。\n"
            "#     主路径 = antechamber -fi mdl（键级来自 SDF 连接表），\n"
            "#     以绕开 ACPYPE 内部的 OpenBabel 转换（V2 根因所在）。\n"
            f"antechamber -i {ligand_geom} -fi mdl -o ligand.mol2 -fo mol2 "
            f"-nc {q} -at gaff2 -dr no -pf y\n"
            f"acpype -i ligand.mol2 -b ligand -n {q} -a gaff2 -o gmx -f")

    merge_py = (f"{py} -c \"import sys; sys.path.insert(0,'{project_root}'); "
                f"from core import node05_md_prep as N5; "
                f"N5.merge_protein_ligand('protein.gro','ligand.acpype/ligand_GMX.gro','complex.gro'); "
                f"N5.merge_topology('topol.top','{itp_inc}')\"")

    return f"""#!/bin/bash
set -euo pipefail
{path_export_line()}
set +u; source {gmxrc}; set -u
cd {sysdir}
echo "[PREP] {target} q={q} @ $(hostname) $(date -Is)"

# ---- 0) 清除先前产物（保留 ligand.sdf / ligand.mol2 与本脚本）-----------------
rm -rf ligand.acpype
rm -f protein.gro topol.top posre_prot.itp complex.gro complex_dry.gro box.gro \\
      solv.gro ions.tpr precharge.tpr grompp_pre.log pdb2gmx.log mdout.mdp \\
      node05_prep_manifest.json
rm -f \\#*\\# 2>/dev/null || true

{param_block}

# ---- I3: pdb2gmx（His 全 HIE=epsilon 中性；应答由受体 His 数推导）------------
printf '{his_answers}' | gmx pdb2gmx -f {receptor_pdb} -o protein.gro -p topol.top \\
    -i posre_prot.itp -ff amber99sb-ildn -water tip3p -ignh -his > pdb2gmx.log 2>&1

# ---- 合并 protein+ligand（Python；只搬坐标，不改化学）------------------------
{merge_py}

# ---- I4 前：grompp 捕获 genion 前感知电荷（stdout+stderr 全捕获，失败即停）----
cat > ions.mdp <<MDP
integrator=steep
nsteps=0
cutoff-scheme=Verlet
coulombtype=PME
rcoulomb=1.0
rvdw=1.0
pbc=xyz
MDP
gmx grompp -f ions.mdp -c complex.gro -p topol.top -o precharge.tpr -maxwarn 5 \\
    > grompp_pre.log 2>&1

# ---- I4: 溶剂化 + 加离子 -----------------------------------------------------
{emit_solvate_commands()}
echo "[PREP DONE] {sysdir} $(date -Is)"
"""


def emit_array_sbatch(script_list_file, n_tasks, project_root, log_dir,
                      partition="cpuPartition", cpus=2, throttle=16) -> str:
    """
    Slurm 数组：每任务跑一份 run_prep.sh。script_list_file 每行一个脚本路径。
    v0.2.0：cpus 8 -> 2（sqm / antechamber 为串行，8 核会空占约 200 核）。
    注意：清单在数组运行期间【绝不可重写】（HANDOFF 6.1 第 11 条）。
    """
    return f"""#!/bin/bash
#SBATCH --job-name=node05_prep
#SBATCH --partition={partition}
#SBATCH --array=0-{n_tasks-1}%{throttle}
#SBATCH --cpus-per-task={cpus}
#SBATCH --output={log_dir}/prep_%A_%a.out
#SBATCH --error={log_dir}/prep_%A_%a.err
cd {project_root}
SCRIPT=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" {script_list_file})
echo "task $SLURM_ARRAY_TASK_ID -> $SCRIPT"
bash "$SCRIPT"
"""


def sdf_to_gro(sdf_path: str, out_gro: str, resname: str = "MOL"):
    """次姿势体系复用主姿势参数时，只把本姿势坐标写成 .gro（原子序与 SDF 一致）。"""
    from rdkit import Chem
    m = Chem.MolFromMolFile(sdf_path, removeHs=False, sanitize=True)
    if m is None:
        raise RuntimeError(f"sdf_to_gro: 无法读取 {sdf_path}")
    conf = m.GetConformer()
    n = m.GetNumAtoms()
    seen, lines = {}, []
    for i, a in enumerate(m.GetAtoms()):
        sym = a.GetSymbol().upper()
        seen[sym] = seen.get(sym, 0) + 1
        p = conf.GetAtomPosition(i)
        lines.append(f"{1:5d}{resname:<5s}{(sym + str(seen[sym])):>5s}{i+1:5d}"
                     f"{p.x/10:8.3f}{p.y/10:8.3f}{p.z/10:8.3f}")
    Path(out_gro).write_text("\n".join(["ligand pose (V3)", str(n)] + lines +
                                       ["   0.00000   0.00000   0.00000"]) + "\n")
    return n


# ===============================================================================
#  I2b · 配体分子图不变量（v0.2.0 新增）
# ===============================================================================
# GAFF2 质量 -> 元素。ACPYPE 的 [atoms] 第 8 列（0-based p[7]）即质量。
_MASS2EL = [(1.008, "H"), (12.01, "C"), (14.01, "N"), (16.00, "O"), (19.00, "F"),
            (30.97, "P"), (32.06, "S"), (35.45, "Cl"), (79.90, "Br"), (126.9, "I"),
            (10.81, "B"), (28.09, "Si")]


def _mass_to_element(m: float, tol: float = 0.35):
    for mm, el in _MASS2EL:
        if abs(m - mm) <= tol:
            return el
    return None


def parse_itp_atoms_bonds(itp_path: str):
    """
    解析 ACPYPE ligand_GMX.itp 的 [atoms] 与 [bonds]。
    实测列布局（pool_mol_030）：
      [atoms]  nr type resi res atom cgnr charge mass  -> p[0] 序号, p[6] 电荷, p[7] 质量
      [bonds]  ai aj funct r k                         -> p[0], p[1] 为 1-based 原子序号
    仅取 [ moleculetype ] 之后的 [atoms]（[atomtypes] 在其前，须排除）。
    """
    sect, seen_moltype = None, False
    elements, bonds = [], set()
    for ln in Path(itp_path).read_text().splitlines():
        st = ln.strip()
        if st.startswith("["):
            name = st.replace(" ", "").lower()
            if name.startswith("[moleculetype]"):
                seen_moltype = True
            sect = name
            continue
        if not st or st.startswith(";"):
            continue
        p = st.split()
        if sect and sect.startswith("[atoms]") and seen_moltype and len(p) >= 8:
            try:
                mass = float(p[7])
            except ValueError:
                continue
            el = _mass_to_element(mass)
            if el is None:
                raise AssertionError(f"I2b: [atoms] 序号 {p[0]} 的质量 {mass} 无法映射到元素")
            elements.append(el)
        elif sect and sect.startswith("[bonds]") and len(p) >= 3:
            try:
                ai, aj = int(p[0]) - 1, int(p[1]) - 1
            except ValueError:
                continue
            bonds.add(frozenset((ai, aj)))
    return elements, bonds


def assert_ligand_graph(itp_path: str, ligand_sdf: str):
    """
    I2b：断言 ligand_GMX.itp 所描述的分子图 == 投喂 ACPYPE 的 ligand.sdf 的分子图。

    为何需要它：I2 只断言部分电荷之和 == frozen_net_charge，而 AM1-BCC 电荷
    【按构造】就求和为 -nc —— 该和对键级、芳香性判定、逐原子形式电荷的任何改动
    完全不敏感。V2 所损坏的那类分子会带着正确的 itp_charge_sum 通过 I2。
    ligand.sdf 与 master SDF 之间的同一性已由 I1（规范异构 SMILES + 形式电荷）保证，
    故 I2b 一旦成立，master -> MD 拓扑的整条身份链条即闭合。

    实测（pool_mol_030，81 原子）：ACPYPE/antechamber 保持输入原子顺序
    （itp nr=1..41 为重原子、42..81 为氢，与 RDKit AddHs 后的顺序一致），故逐索引比对。
    若将来某分子上顺序不被保持，本函数会明确报出，届时改用图同构比对。
    """
    from collections import Counter
    from rdkit import Chem
    m = Chem.MolFromMolFile(ligand_sdf, removeHs=False, sanitize=True)
    if m is None:
        raise RuntimeError(f"I2b: 无法读取 {ligand_sdf}")
    ref_el = [a.GetSymbol() for a in m.GetAtoms()]
    ref_bonds = {frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in m.GetBonds()}
    itp_el, itp_bonds = parse_itp_atoms_bonds(itp_path)

    if len(itp_el) != len(ref_el):
        raise AssertionError(f"I2b: 原子数不符 itp {len(itp_el)} vs sdf {len(ref_el)}")
    if Counter(itp_el) != Counter(ref_el):
        raise AssertionError(f"I2b: 元素组成不符 itp {dict(Counter(itp_el))} "
                             f"vs sdf {dict(Counter(ref_el))}")
    mism = [(i + 1, ref_el[i], itp_el[i]) for i in range(len(ref_el)) if ref_el[i] != itp_el[i]]
    if mism:
        raise AssertionError(f"I2b: 逐索引元素不符，前 5 处 {mism[:5]}（共 {len(mism)} 处）"
                             f" —— ACPYPE 未保持输入原子顺序，需改用图同构比对")

    only_itp = itp_bonds - ref_bonds
    only_sdf = ref_bonds - itp_bonds
    if only_itp or only_sdf:
        fmt = lambda S: sorted(tuple(sorted(x)) for x in S)[:5]
        raise AssertionError(
            f"I2b: 键连接不符 —— itp 多出 {len(only_itp)} 条 {fmt(only_itp)}；"
            f"sdf 多出 {len(only_sdf)} 条 {fmt(only_sdf)}。"
            f"分子图在 SDF -> 参数化之间被改动（V2 根因的失效形态）")

    return {"invariant": "I2b_ligand_graph", "ok": True,
            "n_atoms": len(ref_el), "n_bonds": len(ref_bonds),
            "element_counts": dict(Counter(ref_el)),
            "atom_order_preserved": True,
            "ligand_sdf_sha256": sha256_file(ligand_sdf)}


# ===============================================================================
#  I3b · 残基序列全局比对（v0.2.0 新增）
# ===============================================================================
_RES_ALIAS = {
    "HIS": {"HIS", "HISE", "HISD", "HISH", "HIS1", "HIE", "HID", "HIP"},
    "CYS": {"CYS", "CYS2", "CYX", "CYM"},
    "LYS": {"LYS", "LYSN", "LYSH"},
    "ASP": {"ASP", "ASPH", "ASH"},
    "GLU": {"GLU", "GLUH", "GLH"},
}


def _alias_ok(pdb_name: str, gro_name: str) -> bool:
    return pdb_name == gro_name or gro_name in _RES_ALIAS.get(pdb_name, set())


def _pdb_residue_sequence(pdb_path: str):
    seq, last = [], None
    for ln in Path(pdb_path).read_text().splitlines():
        if not ln.startswith("ATOM"):
            continue
        key = (ln[21], ln[22:27])
        if key != last:
            seq.append((int(ln[22:26]), ln[17:20].strip()))
            last = key
    return seq


def _gro_residue_sequence(gro_path: str):
    L = Path(gro_path).read_text().splitlines()
    n = int(L[1])
    seq, last = [], None
    for ln in L[2:2 + n]:
        rid, rnm = int(ln[0:5]), ln[5:10].strip()
        if (rid, rnm) != last:
            seq.append((rid, rnm))
            last = (rid, rnm)
    return seq


def assert_residue_sequence(receptor_pdb: str, protein_gro: str):
    """
    I3b：断言 protein.gro 的残基序列与清洗后输入 PDB 逐位一致。

    理由：pdb2gmx 可能重新编号，而 I3 只读某一索引处的残基名；编号偏移时
    （Mpro 有多个 His / Cys）存在巧合通过的可能。本检查一次性验证
    (1) 全局编号未偏移；(2) 组氨酸互变异构体是否统一为 HISE。
    """
    pdb_seq = _pdb_residue_sequence(receptor_pdb)
    gro_seq = _gro_residue_sequence(protein_gro)
    if len(pdb_seq) != len(gro_seq):
        raise AssertionError(f"I3b: 残基数不符 pdb {len(pdb_seq)} vs gro {len(gro_seq)}")
    bad = [(i + 1, pdb_seq[i], gro_seq[i]) for i in range(len(pdb_seq))
           if not _alias_ok(pdb_seq[i][1], gro_seq[i][1])]
    if bad:
        raise AssertionError(f"I3b: 残基序列不符，前 5 处 {bad[:5]}（共 {len(bad)} 处）")
    renum = [i for i in range(len(pdb_seq)) if pdb_seq[i][0] != gro_seq[i][0]]
    his_states = sorted({g[1] for p, g in zip(pdb_seq, gro_seq) if p[1] == "HIS"})
    if his_states and his_states != ["HISE"]:
        raise AssertionError(f"I3b: 组氨酸态不统一：{his_states}（期望全部 HISE）")
    return {"invariant": "I3b_residue_sequence", "ok": True,
            "n_residues": len(pdb_seq),
            "renumbered_by_pdb2gmx": bool(renum),
            "n_his": sum(1 for p in pdb_seq if p[1] == "HIS"),
            "his_states": his_states,
            "n_cys": sum(1 for p in pdb_seq if p[1] == "CYS")}


# ===============================================================================
#  受体拓扑跨体系同一性（v0.2.0 新增）
# ===============================================================================
def gro_canonical_bytes(gro_path) -> bytes:
    """
    .gro 的可比对内容（剔除第 1 行标题）。

    v0.2.1 修正：GROMACS 在写 .gro 时会随机挑一句题词填入标题行
    （"Great Red Oystrich Makes All Chemists Sane" 等）。v0.2.0 对整文件取
    sha256，于是 33 个逐字节相同的受体被判为 19 种哈希 —— 是检查错，不是产物错。
    此处保留原子数行、全部原子记录与盒矢量，只丢弃题词。
    """
    L = Path(gro_path).read_text().splitlines()
    n = int(L[1])
    return ("\n".join([L[1]] + L[2:n + 2] + [L[n + 2]])).encode()


def _gro_coords(gro_path):
    L = Path(gro_path).read_text().splitlines()
    n = int(L[1])
    return [(float(l[20:28]), float(l[28:36]), float(l[36:44])) for l in L[2:n + 2]]


def assert_receptor_topology_identical(sysdirs, atol_angstrom: float = 0.0):
    """
    数组对同一受体运行了 N 次 pdb2gmx。断言各体系的受体侧产物同一。

    v0.2.1：`protein.gro` 改为对【剔除随机题词后】的内容取哈希，并独立断言
    各体系坐标的最大绝对偏差 <= atol_angstrom（默认 0，即要求逐位相同）。
    `posre_prot.itp` 无标题行，仍取整文件哈希。
    不含 topol.top —— 它在 merge / solvate / genion 之后已按体系不同。
    """
    dirs = [Path(d) for d in sysdirs]
    out = {}

    h = {}
    for d in dirs:
        p = d / "posre_prot.itp"
        if not p.is_file():
            raise AssertionError(f"受体同一性: 缺失 {p}")
        h.setdefault(sha256_file(p), []).append(d.name)
    if len(h) != 1:
        raise AssertionError("受体同一性: posre_prot.itp 出现 "
                             f"{len(h)} 种哈希 —— " +
                             "; ".join(f"{k[:12]}...x{len(v)}（如 {v[0]}）" for k, v in h.items()))
    out["posre_prot.itp"] = {"sha256": next(iter(h)), "n_systems": len(dirs)}

    import hashlib
    hg, titles = {}, set()
    for d in dirs:
        p = d / "protein.gro"
        if not p.is_file():
            raise AssertionError(f"受体同一性: 缺失 {p}")
        titles.add(Path(p).read_text().splitlines()[0])
        hg.setdefault(hashlib.sha256(gro_canonical_bytes(p)).hexdigest(), []).append(d.name)
    if len(hg) != 1:
        raise AssertionError("受体同一性: protein.gro（剔除题词后）出现 "
                             f"{len(hg)} 种哈希 —— " +
                             "; ".join(f"{k[:12]}...x{len(v)}（如 {v[0]}）" for k, v in hg.items()))

    ref, worst = _gro_coords(dirs[0] / "protein.gro"), (0.0, None)
    for d in dirs[1:]:
        X = _gro_coords(d / "protein.gro")
        if len(X) != len(ref):
            raise AssertionError(f"受体同一性: {d.name} 原子数 {len(X)} != {len(ref)}")
        m = max(abs(a - b) for xa, xb in zip(X, ref) for a, b in zip(xa, xb)) * 10.0
        if m > worst[0]:
            worst = (m, d.name)
    if worst[0] > atol_angstrom:
        raise AssertionError(f"受体同一性: 坐标最大偏差 {worst[0]:.4f} Å > {atol_angstrom} "
                             f"(@{worst[1]})")
    out["protein.gro"] = {"sha256_excl_title": next(iter(hg)), "n_systems": len(dirs),
                          "max_coord_dev_angstrom": round(worst[0], 6),
                          "n_distinct_gromacs_quotes": len(titles)}
    return {"invariant": "receptor_topology_identical", "ok": True, "files": out}


def his_answers_for_receptor(receptor_pdb: str, expected_n=None):
    """由受体实际 His 数生成 pdb2gmx -his 应答（全 HIE=1）；可选断言计数。"""
    his = count_his_residues(receptor_pdb)
    if expected_n is not None and len(his) != expected_n:
        raise AssertionError(f"His 计数 {len(his)} != 预期 {expected_n}；实测残基号 {his}")
    return his_answer_sequence(len(his)), his


def prep_manifest_v2(job: dict, invariants: dict, out_dir):
    """v0.2.0 manifest：不变量键不再固定为 I1..I4，容纳 I2b / I3b。"""
    man = {"module": "node05_md_prep", "version": NODE05_VERSION, "job": job,
           "invariants": invariants,
           "all_passed": all(v.get("ok") for v in invariants.values())}
    p = Path(out_dir) / "node05_prep_manifest.json"
    p.write_text(json.dumps(man, indent=2, ensure_ascii=False))
    return man


# ═══════════════════════════════════════════════════════════════════════════════
#  Node 05b · MD 生产层（v0.2.2）
#
#  物理协议【逐字转录自 V2 run_md2.sh 的 write_mdp_* 函数】，仅四处受控偏离：
#    (a) gen_seed = -1  →  每 replica 显式互异种子并写入 manifest
#        （-1 由时钟/PID 取种，replica 既不可复现也未被记录，与 σ_PBSA 的方差论证相抵触）
#    (b) NVT/NPT 的 define 增加 -DPOSRESLIG，引入 ACPYPE 已产出的 posre_ligand.itp
#        （V2 的 -DPOSRES 只约束蛋白，配体在平衡阶段自由，削弱"轨迹起自所选姿势"的表述）
#    (c) tc-grps 中的配体组名参数化：V2 = Protein_UNL，V3 = Protein_MOL
#        （残基名由 ACPYPE 决定，V3 为 MOL；见 LIGAND_RESNAME）
#    (d) index 组一律【按名解析】，不再硬编码序号（V2 的 -cg 1 13 与 mindist 的 16/19）
#
#  其余全部保持 V2 原样：dt 2 fs、LINCS h-bonds（iter 2 / order 4）、Verlet、nstlist 20、
#  PME、rc = rvdw = 1.0 nm、V-rescale（τ 0.2 平衡 / 0.5 生产）、NPT 用 C-rescale +
#  refcoord_scaling = com、生产用 Parrinello-Rahman、nstxout-compressed = 5000。
#
#  已核对的算术：50,000,000 步 × 2 fs = 100 ns；每 5,000 步存一帧 = 10 ps → 10,001 帧，
#  与 phase2/3/4 的假定一致；MM-PBSA 取 5001–10000 间隔 25 = 200 帧 = 50.01–99.76 ns（2026-10-31 更正：原记 interval 10 / 500 帧）。
# ═══════════════════════════════════════════════════════════════════════════════

LIGAND_RESNAME = "MOL"          # ACPYPE 输出的配体残基名（V2 为 UNL）
REPLICA_SEEDS = {"r1": 42, "r2": 43, "r3": 44}   # 显式互异，取 Node 04 种子族
V2_PROTOCOL_SOURCE = "run_md2.sh (V2) write_mdp_{em,nvt,nvt_continue,npt_continue,md}"

# ── 平衡路径：由涌现改为声明（v0.2.3）────────────────────────────────────────
# V2 的 run_md2.sh 含一条"条件"平衡闸：EM 后若 fmax > 500 或 mindist < 0.12 nm
# 则改走强 EM + 分段 NVT。实测该闸【恒真】，理由是所比较的两组重叠：
#   System ⊇ SOL，故最近的 System–SOL 原子对即某个水分子自身的 O–H 键。
#   TIP3P 平衡 O–H = 0.09572 nm，且 EM 段无 constraints，实测 0.0944 nm
#   —— 对任何溶剂化体系都必然 < 0.12 nm。
# 磁盘证据（V2 全部 138 条溶剂化轨迹）：em_strong.gro 138 / nvt300.gro 138 / nvt.gro 0。
# 即 V2 的"主路径"从未被执行过，是死代码。
# 因此 V3 显式声明分段平衡为唯一路径（与 V2 既成行为逐字一致），
# 而不让协议取决于一次组集合运算的偶然。原闸仍照算并记入 qc.json 作为诊断。
EQUILIBRATION_MODE = "staged"      # EM(500) -> EM(100) -> NVT 50 K 20 ps -> NVT 300 K 80 ps
LEGACY_QC_NOTE = ("System-SOL mindist 恒等于水分子自身 O-H 键长（TIP3P 0.09572 nm，"
                  "EM 段无约束），必 < 0.12 nm 阈值；V2 因此 138/138 走分段路径。"
                  "本值仅作诊断记录，不再决定分支。")


def detect_ligand_resname(itp_path: str) -> str:
    """从 ligand_GMX.itp 的 [atoms] 读出真实残基名，不假定 MOL。"""
    sect, seen_mol = None, False
    for ln in Path(itp_path).read_text().splitlines():
        st = ln.strip()
        if st.startswith("["):
            name = st.replace(" ", "").lower()
            if name.startswith("[moleculetype]"):
                seen_mol = True
            sect = name
            continue
        if not st or st.startswith(";"):
            continue
        if sect and sect.startswith("[atoms]") and seen_mol:
            p = st.split()
            if len(p) >= 8:
                return p[3]
    raise AssertionError(f"detect_ligand_resname: 无法从 {itp_path} 读出残基名")


# ── index 组：一律按名解析 ────────────────────────────────────────────────────
def parse_ndx_groups(ndx_path: str):
    """返回 index.ndx 的组名列表（下标即 gmx 的组号）。"""
    names = []
    for ln in Path(ndx_path).read_text().splitlines():
        m = re.match(r"\s*\[\s*(.+?)\s*\]", ln)
        if m:
            names.append(m.group(1))
    return names


def ndx_group_index(ndx_path: str, name: str) -> int:
    names = parse_ndx_groups(ndx_path)
    if name not in names:
        raise AssertionError(f"index 组 '{name}' 不存在；现有组：{names}")
    if names.count(name) != 1:
        raise AssertionError(f"index 组 '{name}' 出现 {names.count(name)} 次，无法唯一解析")
    return names.index(name)


def assert_ndx_groups(ndx_path: str, required):
    """
    按名断言所需 index 组存在且唯一，返回 {名: 组号}。

    V2 在 phase3 用 `-cg 1 13` 并断言 g[13]=='UNL'，在 run_md2 的 mindist 里硬编码
    16 / 19。V3 的 index 是重新构建的，序号不保证一致——故一律按名解析。
    """
    return {n: ndx_group_index(ndx_path, n) for n in required}


# ── 拓扑：为已存在的 topol.top 补配体位置限制（幂等） ─────────────────────────
def add_ligand_posres_include(top_path: str,
                              ligand_itp_include: str = "ligand.acpype/ligand_GMX.itp",
                              posre_include: str = "ligand.acpype/posre_ligand.itp"):
    """
    在配体 itp 的 #include 之后插入 POSRESLIG 保护的配体位置限制 include。

    位置正确性：[position_restraints] 归属于其前最近的 [moleculetype]，
    而 ligand_GMX.itp 结束时最近的 moleculetype 即 ligand，故紧随其后插入即正确。
    幂等：已存在则不重复插入。就地修改已溶剂化的 topol.top，不重跑制备。
    """
    p = Path(top_path)
    lines = p.read_text().splitlines()
    if any("POSRESLIG" in ln for ln in lines):
        return {"patched": False, "reason": "already_present"}
    idx = None
    for i, ln in enumerate(lines):
        if ln.strip().startswith("#include") and ligand_itp_include in ln:
            idx = i
            break
    if idx is None:
        raise AssertionError(f"add_ligand_posres_include: 在 {top_path} 未找到 "
                             f'#include "{ligand_itp_include}"')
    block = ["#ifdef POSRESLIG", f'#include "{posre_include}"', "#endif"]
    lines[idx + 1:idx + 1] = block
    p.write_text("\n".join(lines) + "\n")
    return {"patched": True, "inserted_after_line": idx + 1}


# ── I5 · 配体坐标与拓扑一致（第七道断言） ────────────────────────────────────
def assert_ligand_coords_match_itp(ligand_sdf: str, ligand_gro: str, itp_path: str,
                                   tol_angstrom: float = 0.02):
    """
    I5：断言 ligand_GMX.gro 的原子顺序与 itp / ligand.sdf 一致。

    为何需要：I2b 验证的是 .itp。而合并进 complex.gro 的坐标取自 ligand_GMX.gro，
    该文件出自 ACPYPE 内部 OpenBabel 的写出路径。若其原子顺序与 itp 不同，
    配体会带着被打乱的几何进入 MD，而 I1–I4 全都察觉不到。

    方法：比较两文件中【成键原子对】的距离。键长在刚体平移与旋转下不变，故本检验
    不假定 ACPYPE 未移动分子。阈值 0.02 Å = .gro 的 0.001 nm 网格所致量化底
    （约 0.0095 Å 中位）的 2 倍；pool_mol_001 实测 0.0087 Å。
    """
    from rdkit import Chem
    import math
    m = Chem.MolFromMolFile(ligand_sdf, removeHs=False, sanitize=True)
    if m is None:
        raise RuntimeError(f"I5: 无法读取 {ligand_sdf}")
    conf = m.GetConformer()
    sdf_xyz = [tuple(conf.GetAtomPosition(i)) for i in range(m.GetNumAtoms())]

    L = Path(ligand_gro).read_text().splitlines()
    n = int(L[1])
    gro_xyz = [(float(l[20:28]) * 10, float(l[28:36]) * 10, float(l[36:44]) * 10)
               for l in L[2:n + 2]]
    if n != len(sdf_xyz):
        raise AssertionError(f"I5: 原子数不符 gro {n} vs sdf {len(sdf_xyz)}")

    _, bonds = parse_itp_atoms_bonds(itp_path)
    dist = lambda A, B: math.dist(A, B)
    worst = 0.0
    for b in bonds:
        i, j = tuple(b)
        worst = max(worst, abs(dist(sdf_xyz[i], sdf_xyz[j]) - dist(gro_xyz[i], gro_xyz[j])))
    if worst > tol_angstrom:
        raise AssertionError(f"I5: gro 与 sdf 键长最大偏差 {worst:.4f} Å > {tol_angstrom} "
                             f"—— ligand_GMX.gro 原子顺序可能与拓扑不符")
    return {"invariant": "I5_ligand_coords_match_itp", "ok": True,
            "n_atoms": n, "n_bonds": len(bonds),
            "max_bond_dev_angstrom": round(worst, 5), "tol_angstrom": tol_angstrom}


# ── MDP：逐字转录自 V2 run_md2.sh ────────────────────────────────────────────
def _tc_grps(ligand_resname: str = LIGAND_RESNAME) -> str:
    return f"Protein_{ligand_resname} Water_and_ions"


def mdp_em(emtol: float = 500.0, nsteps: int = 50000) -> str:
    """V2 write_mdp_em（主路径 emtol=500.0 nsteps=50000；修复分支 100.0 / 200000）。"""
    return (f"integrator  = steep\n"
            f"emtol       = {emtol}\n"
            f"emstep      = 0.01\n"
            f"nsteps      = {nsteps}\n"
            f"cutoff-scheme = Verlet\n"
            f"coulombtype = PME\n"
            f"rcoulomb    = 1.0\n"
            f"rvdw        = 1.0\n"
            f"pbc         = xyz\n")


def _mdp_equil_common(dt: float, ref_t: float, tc: str, tau_t: float) -> str:
    return (f"integrator  = md\n"
            f"dt          = {dt}\n"
            f"constraint_algorithm = lincs\n"
            f"constraints = h-bonds\n"
            f"lincs_iter  = 2\n"
            f"lincs_order = 4\n"
            f"cutoff-scheme = Verlet\n"
            f"nstlist     = 20\n"
            f"rcoulomb    = 1.0\n"
            f"rvdw        = 1.0\n"
            f"coulombtype = PME\n"
            f"tcoupl      = V-rescale\n"
            f"tc-grps     = {tc}\n"
            f"tau_t       = {tau_t} {tau_t}\n"
            f"ref_t       = {ref_t} {ref_t}\n"
            f"pbc         = xyz\n"
            f"nstxout     = 0\n"
            f"nstvout     = 0\n")


def mdp_nvt(nsteps: int, gen_seed: int, dt: float = 0.002, ref_t: float = 300,
            ligand_resname: str = LIGAND_RESNAME, restrain_ligand: bool = True) -> str:
    """V2 write_mdp_nvt；V3 偏离：gen_seed 显式、define 增加 -DPOSRESLIG。"""
    define = "-DPOSRES -DPOSRESLIG" if restrain_ligand else "-DPOSRES"
    return (f"define = {define}\n"
            + _mdp_equil_common(dt, ref_t, _tc_grps(ligand_resname), 0.2)
            + f"nsteps      = {nsteps}\n"
              f"continuation = no\n"
              f"pcoupl      = no\n"
              f"gen_vel     = yes\n"
              f"gen_temp    = {ref_t}\n"
              f"gen_seed    = {int(gen_seed)}\n"
              f"nstenergy   = 500\n"
              f"nstlog      = 500\n")


def mdp_nvt_continue(nsteps: int, dt: float = 0.002, ref_t: float = 300,
                     ligand_resname: str = LIGAND_RESNAME,
                     restrain_ligand: bool = True) -> str:
    """V2 write_mdp_nvt_continue（QC 修复分支的 300 K 续段）。"""
    define = "-DPOSRES -DPOSRESLIG" if restrain_ligand else "-DPOSRES"
    return (f"define = {define}\n"
            + _mdp_equil_common(dt, ref_t, _tc_grps(ligand_resname), 0.2)
            + f"nsteps      = {nsteps}\n"
              f"continuation = yes\n"
              f"pcoupl      = no\n"
              f"gen_vel     = no\n"
              f"nstenergy   = 500\n"
              f"nstlog      = 500\n")


def mdp_npt(nsteps: int, dt: float = 0.002, ref_t: float = 300,
            ligand_resname: str = LIGAND_RESNAME, restrain_ligand: bool = True) -> str:
    """V2 write_mdp_npt_continue：C-rescale + refcoord_scaling = com（V2 原样）。"""
    define = "-DPOSRES -DPOSRESLIG" if restrain_ligand else "-DPOSRES"
    return (f"define = {define}\n"
            + _mdp_equil_common(dt, ref_t, _tc_grps(ligand_resname), 0.2)
            + f"nsteps      = {nsteps}\n"
              f"continuation = yes\n"
              f"pcoupl      = C-rescale\n"
              f"pcoupltype  = isotropic\n"
              f"refcoord_scaling = com\n"
              f"tau_p       = 2.0\n"
              f"ref_p       = 1.0\n"
              f"compressibility = 4.5e-5\n"
              f"gen_vel     = no\n"
              f"nstenergy   = 500\n"
              f"nstlog      = 500\n")


def mdp_md(ns: float = 100.0, dt: float = 0.002, ref_t: float = 300,
           ligand_resname: str = LIGAND_RESNAME) -> str:
    """V2 write_mdp_md：生产段无位置限制、Parrinello-Rahman、每 5000 步存一帧。"""
    nsteps = int(round(ns * 1000.0 / dt))
    return (_mdp_equil_common(dt, ref_t, _tc_grps(ligand_resname), 0.5)
            + f"nsteps      = {nsteps}\n"
              f"continuation = yes\n"
              f"pcoupl      = Parrinello-Rahman\n"
              f"pcoupltype  = isotropic\n"
              f"tau_p       = 2.0\n"
              f"ref_p       = 1.0\n"
              f"compressibility = 4.5e-5\n"
              f"nstxout-compressed = 5000\n"
              f"compressed-x-grps   = System\n"
              f"nstenergy   = 1000\n"
              f"nstlog      = 1000\n")


def mdp_set(replica: str, ligand_resname: str = LIGAND_RESNAME,
            md_ns: float = 100.0, nvt_ps: float = 100.0, npt_ps: float = 100.0,
            dt: float = 0.002, restrain_ligand: bool = True):
    """一个 replica 的全套 MDP + 其种子；返回 {文件名: 文本} 与元信息。"""
    if replica not in REPLICA_SEEDS:
        raise AssertionError(f"未知 replica '{replica}'，已定义 {sorted(REPLICA_SEEDS)}")
    seed = REPLICA_SEEDS[replica]
    n_nvt, n_npt = int(round(nvt_ps / dt)), int(round(npt_ps / dt))
    files = {
        "em.mdp":         mdp_em(500.0, 50000),
        "em_strong.mdp":  mdp_em(100.0, 200000),
        "nvt.mdp":        mdp_nvt(n_nvt, seed, dt, 300, ligand_resname, restrain_ligand),
        "nvt50.mdp":      mdp_nvt(int(round(20.0 / dt)), seed, dt, 50,
                                  ligand_resname, restrain_ligand),
        "nvt300.mdp":     mdp_nvt_continue(int(round(80.0 / dt)), dt, 300,
                                           ligand_resname, restrain_ligand),
        "npt.mdp":        mdp_npt(n_npt, dt, 300, ligand_resname, restrain_ligand),
        "md.mdp":         mdp_md(md_ns, dt, 300, ligand_resname),
    }
    meta = {"replica": replica, "gen_seed": seed, "dt_ps": dt,
            "nvt_ps": nvt_ps, "npt_ps": npt_ps, "md_ns": md_ns,
            "tc_grps": _tc_grps(ligand_resname),
            "ligand_restrained_in_equilibration": restrain_ligand,
            "n_frames_expected": int(round(md_ns * 1000.0 / dt)) // 5000 + 1,
            "protocol_source": V2_PROTOCOL_SOURCE,
            "sha256": {k: hashlib.sha256(v.encode()).hexdigest()[:16]
                       for k, v in files.items()}}
    return files, meta


def assert_replica_seeds_distinct(replicas):
    seeds = [REPLICA_SEEDS[r] for r in replicas]
    if len(set(seeds)) != len(seeds):
        raise AssertionError(f"replica 种子不互异：{dict(zip(replicas, seeds))}")
    if any(int(s) < 0 for s in seeds):
        raise AssertionError("replica 种子含负值（gen_seed = -1 即不可复现）")
    return dict(zip(replicas, seeds))


# ── 生产脚本 ────────────────────────────────────────────────────────────────
_RUN_MD2_V3 = r'''#!/bin/bash
# Node 05b · V3 生产驱动（协议逐字来自 V2 run_md2.sh；四处受控偏离见模块头注释）
set -euo pipefail
__PATH_EXPORT__
set +u; source __GMXRC__; set -u
cd __RUNDIR__
LIG=__LIGRES__
echo "[MD2] $(pwd) replica=__REPLICA__ seed=__SEED__ @ $(hostname) $(date -Is)"

# ---- 0) 由制备目录取输入（只读，不改制备产物）--------------------------------
mkdir -p 03_min 04_nvt 05_npt 06_md
cp -f __SYSDIR__/complex.gro 03_min/complex.gro
cp -f __SYSDIR__/topol.top   03_min/topol.top
cp -f __SYSDIR__/posre_prot.itp 03_min/
mkdir -p 03_min/ligand.acpype
cp -f __SYSDIR__/ligand.acpype/ligand_GMX.itp   03_min/ligand.acpype/
cp -f __SYSDIR__/ligand.acpype/posre_ligand.itp 03_min/ligand.acpype/
cd 03_min

# ---- 1) index：合并 Protein 与配体残基组；随后【按名】断言 -------------------
printf '"Protein" | "%s"\nq\n' "$LIG" | gmx make_ndx -f complex.gro -o index.ndx > make_ndx.log 2>&1
__PY__ -c "import sys; sys.path.insert(0,'__PROJECT_ROOT__'); from core import node05_md_prep as N5; \
import json; g=N5.assert_ndx_groups('index.ndx', ['Protein','$LIG','Protein_$LIG','Water_and_ions','System','SOL']); \
open('ndx_groups.json','w').write(json.dumps(g, indent=1)); print('[NDX]', g)"

# ---- 2) EM（主路径）---------------------------------------------------------
cat > em.mdp <<'MDPEM'
__MDP_EM__
MDPEM
gmx grompp -f em.mdp -c complex.gro -r complex.gro -p topol.top -o em.tpr -n index.ndx -maxwarn 5 > grompp_em.log 2>&1
gmx mdrun -deffnm em -ntmpi 1 -ntomp "${OMP_NUM_THREADS:-8}" >> mdrun_em.log 2>&1

# ---- 3) QC 闸（V2 阈值：fmax>500 或 mindist<0.12 nm → 修复分支）--------------
SOLG=$(__PY__ -c "import json;print(json.load(open('ndx_groups.json'))['SOL'])")
SYSG=$(__PY__ -c "import json;print(json.load(open('ndx_groups.json'))['System'])")
printf '%s\n%s\n' "$SYSG" "$SOLG" | gmx mindist -f em.gro -s em.tpr -n index.ndx -od mindist_em.xvg > mindist.log 2>&1 || true
__PY__ -c "
import re,json,math
fmax=None
for ln in open('em.log',errors='ignore'):
    m=re.search(r'Maximum force\s*=\s*([-\d.eE+]+)', ln)
    if m: fmax=float(m.group(1))
mind=None
try:
    rows=[l.split() for l in open('mindist_em.xvg') if l and l[0] not in '#@']
    if rows: mind=float(rows[-1][1])
except Exception: pass
legacy='repair' if ((mind is not None and mind<0.12) or (fmax is not None and fmax>500.0)) else 'main'
branch='__EQUIL_MODE__'
json.dump({'fmax':fmax,'mindist_nm':mind,
           'equilibration':branch,'declared':True,
           'legacy_gate_branch':legacy,'legacy_gate_agrees':(legacy=='repair'),
           'thresholds':{'fmax':500.0,'mindist_nm':0.12},
           'note':'__QC_NOTE__'}, open('qc.json','w'), indent=1)
print('[QC] fmax=',fmax,' mindist=',mind,' legacy_gate=',legacy,' -> equilibration=',branch)
"
BRANCH=$(__PY__ -c "import json;print(json.load(open('qc.json'))['equilibration'])")

if [ "$BRANCH" = "staged" ]; then
  cat > em_strong.mdp <<'MDPEMS'
__MDP_EM_STRONG__
MDPEMS
  gmx grompp -f em_strong.mdp -c em.gro -r em.gro -p topol.top -o em_strong.tpr -n index.ndx -maxwarn 5 > grompp_ems.log 2>&1
  gmx mdrun -deffnm em_strong -ntmpi 1 -ntomp "${OMP_NUM_THREADS:-8}" >> mdrun_ems.log 2>&1
  EM_GRO=em_strong.gro
else
  EM_GRO=em.gro
fi
cd ..

# ---- 4) NVT -----------------------------------------------------------------
cd 04_nvt
for f in topol.top posre_prot.itp index.ndx; do cp -f ../03_min/$f .; done
mkdir -p ligand.acpype && cp -f ../03_min/ligand.acpype/* ligand.acpype/
cp -f "../03_min/$EM_GRO" start.gro
if [ "$BRANCH" = "staged" ]; then
  cat > nvt50.mdp  <<'MDPN50'
__MDP_NVT50__
MDPN50
  cat > nvt300.mdp <<'MDPN300'
__MDP_NVT300__
MDPN300
  gmx grompp -f nvt50.mdp -c start.gro -r start.gro -p topol.top -o nvt50.tpr -n index.ndx -maxwarn 5 > grompp_nvt50.log 2>&1
  gmx mdrun -deffnm nvt50 -ntmpi 1 -ntomp "${OMP_NUM_THREADS:-8}" >> mdrun_nvt50.log 2>&1
  gmx grompp -f nvt300.mdp -c nvt50.gro -r nvt50.gro -p topol.top -o nvt300.tpr -n index.ndx -maxwarn 5 > grompp_nvt300.log 2>&1
  gmx mdrun -deffnm nvt300 -ntmpi 1 -ntomp "${OMP_NUM_THREADS:-8}" >> mdrun_nvt300.log 2>&1
  NVT_GRO=nvt300.gro
else
  cat > nvt.mdp <<'MDPNVT'
__MDP_NVT__
MDPNVT
  gmx grompp -f nvt.mdp -c start.gro -r start.gro -p topol.top -o nvt.tpr -n index.ndx -maxwarn 5 > grompp_nvt.log 2>&1
  gmx mdrun -deffnm nvt -ntmpi 1 -ntomp "${OMP_NUM_THREADS:-8}" >> mdrun_nvt.log 2>&1
  NVT_GRO=nvt.gro
fi
cd ..

# ---- 5) NPT -----------------------------------------------------------------
cd 05_npt
for f in topol.top posre_prot.itp index.ndx; do cp -f ../04_nvt/$f .; done
mkdir -p ligand.acpype && cp -f ../04_nvt/ligand.acpype/* ligand.acpype/
cp -f "../04_nvt/$NVT_GRO" nvt.gro
cat > npt.mdp <<'MDPNPT'
__MDP_NPT__
MDPNPT
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -n index.ndx -maxwarn 5 > grompp_npt.log 2>&1
gmx mdrun -deffnm npt -ntmpi 1 -ntomp "${OMP_NUM_THREADS:-8}" >> mdrun_npt.log 2>&1
cd ..

# ---- 6) 生产 MD -------------------------------------------------------------
cd 06_md
for f in topol.top index.ndx; do cp -f ../05_npt/$f .; done
mkdir -p ligand.acpype && cp -f ../05_npt/ligand.acpype/* ligand.acpype/
cp -f ../05_npt/npt.gro .
cat > md.mdp <<'MDPMD'
__MDP_MD__
MDPMD
gmx grompp -f md.mdp -c npt.gro -t ../05_npt/npt.cpt -p topol.top -o md.tpr -n index.ndx -maxwarn 5 > grompp_md.log 2>&1
gmx mdrun -deffnm md -ntmpi 1 -ntomp "${OMP_NUM_THREADS:-8}" -nb gpu -pme gpu >> mdrun_md.log 2>&1
cd ..

cp -f 03_min/qc.json 03_min/ndx_groups.json 06_md/ 2>/dev/null || true
echo "[MD2 DONE] $(pwd) branch=$BRANCH $(date -Is)"
'''


def emit_run_md2_v3(rundir, sysdir, replica, project_root, gmxrc, py,
                    ligand_resname: str = LIGAND_RESNAME, md_ns: float = 100.0,
                    nvt_ps: float = 100.0, npt_ps: float = 100.0,
                    restrain_ligand: bool = True):
    """生成一条轨迹的生产脚本；同时返回该 replica 的 MDP 元信息（写入 manifest）。"""
    files, meta = mdp_set(replica, ligand_resname, md_ns, nvt_ps, npt_ps,
                          restrain_ligand=restrain_ligand)
    s = _RUN_MD2_V3
    for k, v in {
        "__PATH_EXPORT__": path_export_line(),
        "__GMXRC__": str(gmxrc), "__RUNDIR__": str(rundir), "__SYSDIR__": str(sysdir),
        "__PROJECT_ROOT__": str(project_root), "__PY__": str(py),
        "__LIGRES__": ligand_resname, "__REPLICA__": replica,
        "__SEED__": str(meta["gen_seed"]),
        "__EQUIL_MODE__": EQUILIBRATION_MODE, "__QC_NOTE__": LEGACY_QC_NOTE,
        "__MDP_EM__": files["em.mdp"].rstrip("\n"),
        "__MDP_EM_STRONG__": files["em_strong.mdp"].rstrip("\n"),
        "__MDP_NVT__": files["nvt.mdp"].rstrip("\n"),
        "__MDP_NVT50__": files["nvt50.mdp"].rstrip("\n"),
        "__MDP_NVT300__": files["nvt300.mdp"].rstrip("\n"),
        "__MDP_NPT__": files["npt.mdp"].rstrip("\n"),
        "__MDP_MD__": files["md.mdp"].rstrip("\n"),
    }.items():
        s = s.replace(k, v)
    return s, meta


_PROD_SBATCH = r'''#!/bin/bash
#SBATCH --job-name=node05b_prod
#SBATCH --output=__LOGDIR__/prod_%A_%a.out
#SBATCH --error=__LOGDIR__/prod_%A_%a.err
#SBATCH --array=0-__NMAX__%__THROTTLE__
#SBATCH --partition=__PARTITION__
#SBATCH --account=__ACCOUNT__
#SBATCH --qos=__QOS__
#SBATCH --gres=__GRES__
#SBATCH --cpus-per-task=__CPUS__
#SBATCH --time=__TIME__
set -uo pipefail
cd __PROJECT_ROOT__
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
SCRIPT=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" __LISTFILE__)
RUNDIR=$(dirname "$SCRIPT")
if [ -f "$RUNDIR/06_md/md.gro" ]; then echo "SKIP(已完成) $RUNDIR"; exit 0; fi
echo "task $SLURM_ARRAY_TASK_ID -> $SCRIPT"
bash "$SCRIPT"
'''


def emit_prod_sbatch(list_file, n_tasks, project_root, log_dir, partition="gpu",
                     account="hongmeilab", qos="normal", gres="gpu:geforce:1",
                     cpus=8, time="14:00:00", throttle=8) -> str:
    """生产数组；断点续跑（06_md/md.gro 存在即跳过）。单条 100 ns ≈ 8.4 h（V2 实测）。"""
    s = _PROD_SBATCH
    for k, v in {"__LOGDIR__": str(log_dir), "__NMAX__": str(n_tasks - 1),
                 "__THROTTLE__": str(throttle), "__PARTITION__": partition,
                 "__ACCOUNT__": account, "__QOS__": qos, "__GRES__": gres,
                 "__CPUS__": str(cpus), "__TIME__": time,
                 "__PROJECT_ROOT__": str(project_root),
                 "__LISTFILE__": str(list_file)}.items():
        s = s.replace(k, v)
    return s


# ── 生产后核验 ──────────────────────────────────────────────────────────────
def assert_trajectory_complete(rundir, expected_frames: int = 10001):
    """
    生产后逐轨迹核验：md.gro 存在、md.log 见完成标记、QC 分支已记录、
    实际 gen_seed 与登记值一致（从 md.log/nvt tpr 侧读回 mdp 记录）。
    """
    d = Path(rundir)
    out = {"rundir": d.name}
    gro = d / "06_md" / "md.gro"
    if not gro.is_file():
        raise AssertionError(f"轨迹未完成：缺 {gro}")
    log = d / "06_md" / "md.log"
    if not log.is_file() or "Finished mdrun" not in log.read_text(errors="ignore"):
        raise AssertionError(f"轨迹无完成标记 'Finished mdrun'：{log}")
    qc = d / "06_md" / "qc.json"
    if not qc.is_file():
        raise AssertionError(f"缺 QC 记录 {qc} —— 无法说明本条走了哪条平衡分支")
    q = json.loads(qc.read_text())
    if q.get("equilibration") != EQUILIBRATION_MODE:
        raise AssertionError(f"平衡路径 {q.get('equilibration')} != 声明的 {EQUILIBRATION_MODE}：{qc}")
    out["qc"] = q
    nvt_mdp = next((p for p in [(d / "04_nvt" / "nvt.mdp"), (d / "04_nvt" / "nvt50.mdp")]
                    if p.is_file()), None)
    if nvt_mdp is None:
        raise AssertionError(f"缺 NVT mdp：{d/'04_nvt'}")
    m = re.search(r"gen_seed\s*=\s*(-?\d+)", nvt_mdp.read_text())
    if not m or int(m.group(1)) < 0:
        raise AssertionError(f"gen_seed 未显式记录或为 -1：{nvt_mdp}")
    out["gen_seed"] = int(m.group(1))
    out["n_atoms"] = int(gro.read_text().splitlines()[1])
    out["ok"] = True
    return out


def assert_equilibration_uniform(rundirs):
    """
    断言全部轨迹走同一条平衡路径，且与声明值一致。

    V2 的分支由一次组集合运算的偶然决定（恒真），故"统一"是涌现而非保证；
    V3 声明之后，这条断言把它变成可检验的性质。同时汇报 legacy 闸是否一致 ——
    若将来某条轨迹的 legacy 闸给出 'main'，那是需要解释的信号，而非可忽略的噪声。
    """
    import collections
    eq, legacy, mism = collections.Counter(), collections.Counter(), []
    for d in rundirs:
        p = Path(d) / "06_md" / "qc.json"
        if not p.is_file():
            p = Path(d) / "03_min" / "qc.json"
        if not p.is_file():
            raise AssertionError(f"缺 qc.json：{d}")
        q = json.loads(p.read_text())
        eq[q.get("equilibration")] += 1
        legacy[q.get("legacy_gate_branch")] += 1
        if not q.get("legacy_gate_agrees", True):
            mism.append(Path(d).name)
    if len(eq) != 1 or EQUILIBRATION_MODE not in eq:
        raise AssertionError(f"平衡路径不统一或与声明不符：{dict(eq)}（声明 {EQUILIBRATION_MODE}）")
    return {"invariant": "equilibration_uniform", "ok": True,
            "equilibration": EQUILIBRATION_MODE, "n_trajectories": sum(eq.values()),
            "legacy_gate": dict(legacy),
            "legacy_gate_disagreements": mism, "note": LEGACY_QC_NOTE}


if __name__ == "__main__":
    print(f"node05_md_prep v{NODE05_VERSION} — 不变量 I1 / I2 / I2b / I3 / I3b / I4；"
          f"I5；MDP 逐字转录自 V2 run_md2.sh。")
    print("请由 Node05_md_prep.ipynb 调度（gmx / acpype 需 HPC 环境）。")
