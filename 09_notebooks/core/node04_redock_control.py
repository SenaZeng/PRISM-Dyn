# -*- coding: utf-8 -*-
"""Node 04 阳性对照 · 共晶配体自对接（redocking control）
    cd /mnt/hpc/home/wanghongmei/lang2
    /mnt/hpc/.../md_env/bin/python core/node04_redock_control.py

把每个共晶配体重新对接回它【自己】的受体，用刚冻结的同一个盒与受体 PDBQT：
    X77 → 6W63 | 0EN → 3V3M | R30 → 4YLU

这是把"我们精心搭了对接"变成"我们的方案复现了三个已知晶体结构"的那一步，
同时给出每个靶点的天然抑制剂参考分数（候选优于/劣于天然物才是可解读的陈述）。

关键技术点（易错，均已显式处理）：
  * 键级：PDB 的 HETATM 不含键级。从 RCSB 的 ideal SDF 取模板，用
    AssignBondOrdersFromTemplate 赋回；模板文件的 sha256 记入 manifest。
  * 坐标：【保留晶体学坐标】，绝不重新 embed —— 参考姿势就是晶体姿势。
  * 质子化：沿用 ideal SDF 所载的形式电荷，不重新决定（与 Node 03 同一禁令）。
  * RMSD：必须【对称性感知】且【不叠合】。用 rdMolAlign.CalcRMS，
    不用 GetBestRMS（后者会先叠合，那样测的是形状相似而非位置复现）。
    仅比较重原子。判据：最优姿势 RMSD ≤ 2.0 Å。

⚠ 本脚本的化学部分未在本地测试过（无 RDKit/Meeko 环境）。
   每一步都带断言，首次运行请逐条看输出，不要只看最后的 ✅。
"""
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

R = Path("/mnt/hpc/home/wanghongmei/lang2")
sys.path.insert(0, str(R))
from core import node04_docking as N4          # noqa: E402

from rdkit import Chem                          # noqa: E402
from rdkit.Chem import AllChem, rdMolAlign      # noqa: E402

OUT = R / "results" / "04_docking" / "redock_control"
VINA = Path("/mnt/hpc/home/wanghongmei/apps/vina/vina_1.2.7_linux_x86_64")
MK_LIG = Path("/mnt/hpc/home/wanghongmei/anaconda3/envs/md_env/bin/mk_prepare_ligand.py")
RMSD_PASS = 2.0

PAIRS = {"6W63": "X77", "3V3M": "0EN", "4YLU": "R30"}


def extract_het(pdb_path, resname, chain="A"):
    """从原始 PDB 抽取共晶配体的 HETATM（丢 altloc B、丢 H），保留晶体坐标。"""
    keep = []
    for ln in Path(pdb_path).read_text().splitlines():
        if not ln.startswith("HETATM"):
            continue
        if ln[17:20].strip() != resname or ln[21] != chain:
            continue
        if ln[16] not in (" ", "A"):          # altloc
            continue
        if ln[76:78].strip().upper() == "H" or ln[12:16].strip().startswith("H"):
            continue
        keep.append(ln)
    if not keep:
        raise AssertionError(f"{pdb_path} 中找不到 {chain}:{resname} 的 HETATM")
    return "\n".join(keep) + "\nEND\n"


def fetch_template(het, dest):
    """RCSB 的理想构象 SDF —— 只用它的【键级与形式电荷】，不用它的坐标。

    HPC 无外网时：在本地下载后 scp 到 dest 指定位置即可，脚本会直接使用。
    【不】以任何方式猜测或重建键级 —— 猜错会让整个阳性对照失去意义且不可见。
    """
    url = f"https://files.rcsb.org/ligands/download/{het}_ideal.sdf"
    if not dest.is_file():
        try:
            dest.parent.mkdir(parents=True, exist_ok=True)
            urllib.request.urlretrieve(url, dest)
        except Exception as e:                                   # noqa: BLE001
            raise SystemExit(
                f"\n无法下载 {het} 的键级模板（{type(e).__name__}）—— gpu01 无外网。\n"
                f"请在本地执行：\n"
                f"    curl -O {url}\n"
                f"再上传到：\n"
                f"    {dest}\n"
                f"三个模板全部就位后重跑本脚本（已存在的文件不会重新下载）。\n"
                f"注：模板【仅】提供键级与形式电荷，坐标一律用晶体坐标；"
                f"其 sha256 会记入 manifest。\n") from e
    m = Chem.MolFromMolFile(str(dest), removeHs=True)
    if m is None:
        raise AssertionError(f"模板解析失败：{dest}（文件是否完整？）")
    return m


def build_reference(pdb_path, het, work):
    """晶体坐标 + 正确键级 → 参考分子（含 H，H 坐标由 RDKit 加）。"""
    block = extract_het(pdb_path, het)
    (work / f"{het}_crystal.pdb").write_text(block)
    raw = Chem.MolFromPDBBlock(block, removeHs=True, sanitize=False)
    if raw is None:
        raise AssertionError(f"{het}: PDB block 解析失败")
    tmpl = fetch_template(het, work / f"{het}_ideal.sdf")
    if raw.GetNumAtoms() != tmpl.GetNumAtoms():
        raise AssertionError(
            f"{het}: 晶体重原子 {raw.GetNumAtoms()} != 模板 {tmpl.GetNumAtoms()}"
            "（残缺配体或 altloc 处理不当）")
    tmpl_sha = N4.file_sha256(work / f"{het}_ideal.sdf")
    mol = AllChem.AssignBondOrdersFromTemplate(tmpl, raw)
    Chem.SanitizeMol(mol)
    molh = Chem.AddHs(mol, addCoords=True)     # 只补 H 坐标，重原子不动
    sdf = work / f"{het}_ref.sdf"
    w = Chem.SDWriter(str(sdf)); w.write(molh); w.close()
    return mol, sdf, tmpl_sha                  # mol = 重原子参考（晶体姿势）


def prep_ligand(sdf, out_pdbqt):
    """走与 Node 03 同一条 Meeko 路径；禁止任何改动质子化的开关。"""
    cmd = [str(MK_LIG), "-i", str(sdf), "-o", str(out_pdbqt)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not Path(out_pdbqt).is_file():
        raise AssertionError(f"配体 PDBQT 生成失败: {(r.stderr or r.stdout)[:300]}")
    return " ".join(cmd)


def pose_to_mol(pose_pdbqt, ref_mol):
    """把对接姿势读回 RDKit。用 Meeko 官方回转（保留化学身份，不从 PDBQT 读键级）。"""
    from meeko import PDBQTMolecule, RDKitMolCreate
    pm = PDBQTMolecule.from_file(str(pose_pdbqt), skip_typing=True)
    mols = RDKitMolCreate.from_pdbqt_mol(pm)
    mols = [m for m in mols if m is not None]
    if not mols:
        raise AssertionError(f"无法从 {pose_pdbqt} 还原分子")
    out = []
    for m in mols:
        for cid in range(m.GetNumConformers()):
            c = Chem.Mol(m, confId=cid) if m.GetNumConformers() > 1 else m
            out.append(Chem.RemoveHs(c))
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    man = json.loads((R / "results/04_docking/receptors/Node04a_receptor_manifest.json")
                     .read_text(encoding="utf-8"))
    report = {"module": "node04_redock_control", "version": "1.0",
              "created_at_utc": N4.utc_now(), "rmsd_pass_angstrom": RMSD_PASS,
              "rmsd_method": "rdMolAlign.CalcRMS (symmetry-aware, NO superposition), heavy atoms",
              "vina": str(VINA), "vina_sha256": N4.file_sha256(VINA),
              "exhaustiveness": N4.EXHAUSTIVENESS, "num_modes": N4.NUM_MODES,
              "seeds": list(N4.SEEDS), "targets": {}}

    for pid, het in PAIRS.items():
        rc = man["receptors"][pid]
        work = OUT / pid; work.mkdir(parents=True, exist_ok=True)
        print(f"\n── {pid} / {het} ──")

        ref_mol, ref_sdf, tmpl_sha = build_reference(rc["input_pdb"], het, work)
        n_heavy = ref_mol.GetNumAtoms()
        print(f"  参考姿势: {n_heavy} 重原子 (manifest 记 "
              f"{rc['box_center_extra'].get('n_ligand_heavy_atoms')})")
        assert n_heavy == rc["box_center_extra"].get("n_ligand_heavy_atoms"), \
            "重原子数与 Node 04a 记录不符 —— 抽取或 altloc 处理有误"

        lig_pdbqt = work / f"{het}.pdbqt"
        cmd_lig = prep_ligand(ref_sdf, lig_pdbqt)

        per_seed = []
        for seed in N4.SEEDS:
            pose = work / f"{het}__{pid}__s{seed}.pdbqt"
            row = {"ligand": het, "target": pid, "seed": seed,
                   "ligand_pdbqt": str(lig_pdbqt), "ligand_pdbqt_sha256": N4.file_sha256(lig_pdbqt),
                   "receptor_pdbqt": rc["receptor_pdbqt"],
                   "receptor_pdbqt_sha256": rc["conversion"]["pdbqt_sha256"],
                   "box_center": rc["box_center"], "box_size": rc["box_size"],
                   "out_pdbqt": str(pose), "out_json": str(work / f"run_s{seed}.json"),
                   "row_index": seed}
            rec = N4.run_one(row, VINA)
            if rec.get("status") != "ok":
                raise AssertionError(f"{pid}/{het} seed {seed} 对接失败: {rec.get('error')}")
            poses = pose_to_mol(pose, ref_mol)
            rmsds = [rdMolAlign.CalcRMS(p, ref_mol) for p in poses]
            per_seed.append({"seed": seed, "affinity_best": rec["affinity_best"],
                             "rmsd_best_pose": round(rmsds[0], 3),
                             "rmsd_min_over_modes": round(min(rmsds), 3),
                             "n_modes": len(rmsds)})
            print(f"    seed {seed}: ΔG {rec['affinity_best']:>7.2f}  "
                  f"RMSD(top) {rmsds[0]:.2f} Å  RMSD(min over {len(rmsds)} modes) {min(rmsds):.2f} Å")

        top = [x["rmsd_best_pose"] for x in per_seed]
        aff = [x["affinity_best"] for x in per_seed]
        ok = sum(1 for r in top if r <= RMSD_PASS)
        report["targets"][pid] = {
            "cocrystal_ligand": het, "n_heavy_atoms": n_heavy,
            "ligand_prep_cmd": cmd_lig,
            "bond_order_template": f"RCSB {het}_ideal.sdf", "template_sha256": tmpl_sha,
            "per_seed": per_seed,
            "rmsd_top_median": round(sorted(top)[len(top) // 2], 3),
            "n_seeds_within_2A": ok, "n_seeds": len(top),
            "reference_affinity_mean": round(sum(aff) / len(aff), 3),
            "verdict": "PASS" if ok >= 3 else "FAIL",
        }
        print(f"  → {ok}/{len(top)} 个种子的最优姿势在 {RMSD_PASS} Å 内 | "
              f"天然物参考分 {sum(aff)/len(aff):.2f} kcal/mol")

    p = OUT / "redock_control_manifest.json"
    p.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    verdicts = {k: v["verdict"] for k, v in report["targets"].items()}
    print(f"\n{'=' * 66}\n判定: {verdicts}\nmanifest: {p}")
    if all(v == "PASS" for v in verdicts.values()):
        print("✅ 三个已知晶体结构均被复现 —— 盒中心/受体制备/打分链路经阳性对照证实。")
        return 0
    print("❌ 有靶点未复现晶体姿势 —— 在全量前查明原因，不要继续。")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
