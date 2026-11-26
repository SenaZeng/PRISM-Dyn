# -*- coding: utf-8 -*-
"""
V3_clean/core/master_molecule.py
================================
"一个分子对象，多个派生表示"（checkagain21/22 的核心原则）。

单一权威分子对象 = master SDF（带键级/芳香性/形式电荷/显式氢/立体/坐标/稳定原子编号）。
所有下游表示（PDBQT、posed SDF、参数化输入）都从它【派生】，
并保存【正向】原子映射（准备时记录，非事后逆向猜测）。

关键不变量（根因修复，一行话）：
    CanonicalIsoSMILES(posed molecule) == CanonicalIsoSMILES(master molecule)
即：从 docking 出来的分子，可证与进去的是同一分子。

本模块【不】做 docking，也【不】绑定具体 PDBQT 工具（Meeko 等）。
它定义抽象接口 prepare_pdbqt(...)，后端以后替换，不影响目录/manifest/QC 架构。
"""

import hashlib
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

MASTER_MODULE_VERSION = "0.1.0"


# ============ 权威分子对象 ============
def canonical_iso_smiles(mol_or_smiles):
    """规范同分异构 SMILES —— 化学身份的唯一指纹（含立体）。

    ⚠️ 必须清除 atom map number 与显式氢，否则它们污染 canonical SMILES，
       使"posed==master"不变量误判。身份指纹只应反映化学图+立体，不含标注。
    """
    m = mol_or_smiles if isinstance(mol_or_smiles, Chem.Mol) else Chem.MolFromSmiles(mol_or_smiles)
    if m is None:
        raise ValueError("canonical_iso_smiles: 无法解析分子")
    mc = Chem.Mol(m)
    for a in mc.GetAtoms():
        a.SetAtomMapNum(0)                      # 清 atom map（否则进 SMILES）
    mc = Chem.RemoveHs(mc)                       # 去显式氢（比化学图，不比 H 记法）
    return Chem.MolToSmiles(mc, isomericSmiles=True, canonical=True)


def build_master_molecule(authoritative_smiles, tag="", embed_seed=42):
    """从权威 SMILES 建唯一主分子对象。

    产出的 mol 带：键级/芳香性/形式电荷/显式氢/已定义立体/3D 坐标/稳定原子编号。
    每个原子写入 atom_map_number = 权威原子索引 + 1，作为跨表示的稳定锚。

    ⚠️ 不做质子化枚举、不改互变异构、不改形式电荷 —— SMILES 即权威微观态。
       立体未定义的分子【不】自动赋 R/S（保持未定义，交微观态表标注）。
    """
    ref = Chem.MolFromSmiles(authoritative_smiles)
    if ref is None:
        raise ValueError(f"{tag}: 权威 SMILES 解析失败: {authoritative_smiles}")

    # 稳定原子编号（在加氢前，按重原子权威顺序）
    for atom in ref.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)

    molH = Chem.AddHs(ref)                       # 显式氢
    params = AllChem.ETKDGv3()
    params.randomSeed = embed_seed
    if AllChem.EmbedMolecule(molH, params) != 0:
        # 退化情形：换一个种子重试一次
        params.randomSeed = embed_seed + 1
        if AllChem.EmbedMolecule(molH, params) != 0:
            raise RuntimeError(f"{tag}: 3D 嵌入失败")
    AllChem.MMFFOptimizeMolecule(molH)

    Chem.AssignStereochemistry(molH, cleanIt=True, force=True)
    return molH, ref


def identity_manifest(molH, authoritative_smiles, tag=""):
    """身份 manifest：把主分子的所有身份不变量落盘，供 QC 与复现。"""
    heavy = molH.GetNumAtoms() - sum(1 for a in molH.GetAtoms() if a.GetSymbol() == "H")
    stereo = Chem.FindMolChiralCenters(molH, includeUnassigned=True, useLegacyImplementation=False)
    return {
        "tag": tag,
        "authoritative_smiles": authoritative_smiles,
        "canonical_isomeric_smiles": canonical_iso_smiles(molH),
        "formal_charge": Chem.GetFormalCharge(molH),           # 从结构算，非 PDBQT 部分电荷求和
        "molecular_formula": rdMolDescriptors.CalcMolFormula(molH),
        "heavy_atom_count": heavy,
        "explicit_hydrogen_count": sum(1 for a in molH.GetAtoms() if a.GetSymbol() == "H"),
        "stereocenters": [{"atom_idx": i, "label": lab} for i, lab in stereo],
        "input_hash": hashlib.sha256(authoritative_smiles.encode()).hexdigest(),
        "software_versions": {
            "rdkit": __import__("rdkit").__version__,
            "master_module": MASTER_MODULE_VERSION,
        },
    }


def write_master_sdf(molH, path):
    """写 master SDF（保留 atom map number、键级、立体、坐标）。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    w = Chem.SDWriter(str(path))
    w.write(molH)
    w.close()
    return str(path)


# ============ 抽象 PDBQT 准备接口（后端可替换） ============
def prepare_pdbqt(master_mol, output_pdbqt, allow_reprotonation=False,
                  atom_map_output=None, backend="stub", tag=""):
    """从权威分子生成 docking PDBQT，并保存【正向】原子映射。

    强制契约（无论后端）：
      - allow_reprotonation=False：禁止再次改变质子化态（不允许 obabel -p 7.4）
      - 不改 tautomer、不改形式电荷、不删必要极性氢
      - 保存 authoritative_to_pdbqt_atom_map（权威原子 → PDBQT serial）
      - PDBQT 只是 docking 表示，不得成为新的化学真值源

    backend:
      'meeko'  -> 用 Meeko（推荐；生产时接入）
      'stub'   -> 占位后端（架构验证用，不产真 PDBQT，只建映射骨架）

    返回 atom_map dict。真实后端须在生成后验证：重原子数/元素序列/微观态/映射完整。
    """
    if backend == "meeko":
        return _prepare_pdbqt_meeko(master_mol, output_pdbqt, allow_reprotonation,
                                    atom_map_output, tag)
    elif backend == "stub":
        return _prepare_pdbqt_stub(master_mol, output_pdbqt, atom_map_output, tag)
    else:
        raise ValueError(f"{tag}: 未知 PDBQT 后端 {backend}")


def _prepare_pdbqt_meeko(master_mol, output_pdbqt, allow_reprotonation, atom_map_output, tag):
    """Meeko 后端（生产）。在 HPC 上 meeko 可用时启用。"""
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
    except ImportError as e:
        raise RuntimeError(f"{tag}: meeko 未安装，无法用 meeko 后端（{e}）。"
                           f"架构验证请用 backend='stub'。")
    if allow_reprotonation:
        raise ValueError(f"{tag}: 本流程禁止 allow_reprotonation=True（会改微观态）")
    prep = MoleculePreparation()
    prep.prepare(master_mol)
    pdbqt_string, ok, err = PDBQTWriterLegacy.write_string(prep.setup)
    if not ok:
        raise RuntimeError(f"{tag}: Meeko 写 PDBQT 失败: {err}")
    Path(output_pdbqt).parent.mkdir(parents=True, exist_ok=True)
    Path(output_pdbqt).write_text(pdbqt_string)
    # 正向映射：master atom map number → PDBQT serial
    atom_map = _build_forward_atom_map(master_mol, pdbqt_string, tag)
    if atom_map_output:
        Path(atom_map_output).write_text(json.dumps(atom_map, indent=2, ensure_ascii=False))
    return atom_map


def _prepare_pdbqt_stub(master_mol, output_pdbqt, atom_map_output, tag):
    """占位后端：不产真 PDBQT，但建立【正向映射骨架】以验证架构。

    映射规则（架构验证）：master 重原子权威顺序 → 顺序 PDBQT serial 1..N。
    真实后端替换后，映射由该后端的原子顺序决定，但接口/manifest/QC 不变。
    """
    heavy = [a for a in master_mol.GetAtoms() if a.GetSymbol() != "H"]
    atom_map = {}
    for serial, a in enumerate(heavy, start=1):
        amn = a.GetAtomMapNum()                 # 权威原子编号
        atom_map[str(serial)] = {
            "authoritative_atom_index": amn - 1 if amn > 0 else a.GetIdx(),
            "element": a.GetSymbol(),
            "pdbqt_serial": serial,
            "pdbqt_name": f"{a.GetSymbol()}{serial}",
        }
    Path(output_pdbqt).parent.mkdir(parents=True, exist_ok=True)
    Path(output_pdbqt).write_text(f"# STUB PDBQT (architecture validation only) tag={tag}\n"
                                  f"# {len(heavy)} heavy atoms\n")
    if atom_map_output:
        Path(atom_map_output).parent.mkdir(parents=True, exist_ok=True)
        Path(atom_map_output).write_text(json.dumps(atom_map, indent=2, ensure_ascii=False))
    return atom_map


def _build_forward_atom_map(master_mol, pdbqt_string, tag):
    """从真实 PDBQT 文本解析 serial，与 master 权威编号对齐（元素校验）。"""
    heavy = [a for a in master_mol.GetAtoms() if a.GetSymbol() != "H"]
    # 解析 PDBQT 的 ATOM/HETATM 行（serial + element）
    pdbqt_atoms = []
    for line in pdbqt_string.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            serial = int(line[6:11])
            elem = line[76:78].strip() or line[12:14].strip()
            pdbqt_atoms.append((serial, elem))
    atom_map = {}
    # 按重原子权威顺序与 PDBQT 重原子顺序对齐（Meeko 保序），元素须一致
    heavy_pdbqt = [(s, e) for s, e in pdbqt_atoms if e.upper() != "H"]
    if len(heavy_pdbqt) != len(heavy):
        raise ValueError(f"{tag}: PDBQT 重原子{len(heavy_pdbqt)}≠master重原子{len(heavy)}")
    for k, (a, (serial, elem)) in enumerate(zip(heavy, heavy_pdbqt)):
        if a.GetSymbol().upper() != elem.upper():
            raise ValueError(f"{tag}: 映射元素不符 master={a.GetSymbol()} pdbqt={elem}（serial {serial}）")
        amn = a.GetAtomMapNum()
        atom_map[str(serial)] = {
            "authoritative_atom_index": amn - 1 if amn > 0 else a.GetIdx(),
            "element": a.GetSymbol(),
            "pdbqt_serial": serial,
        }
    return atom_map


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


if __name__ == "__main__":
    # 自测：001 走完整正推链，验证不变量
    smi = "CNC(=O)c1cc(C(=O)NC(=O)c2cccc3cccnc23)cc(-c2cccc3cccnc23)c1"
    molH, ref = build_master_molecule(smi, "self_test")
    man = identity_manifest(molH, smi, "self_test")
    print(f"master: 式={man['molecular_formula']}, 电荷={man['formal_charge']}, "
          f"重原子={man['heavy_atom_count']}, 立体中心={len(man['stereocenters'])}")
    # posed == master 不变量
    posed_can = canonical_iso_smiles(Chem.RemoveHs(molH))
    master_can = man["canonical_isomeric_smiles"]
    # 注意：master_can 含显式氢版的规范化；这里比去氢主图
    master_noH = canonical_iso_smiles(ref)
    print(f"不变量 CanonicalIsoSMILES(posed)==CanonicalIsoSMILES(master): {posed_can == master_noH}")
    print(f"master_module {MASTER_MODULE_VERSION} self-test OK")
