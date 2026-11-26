#!/usr/bin/env python3
"""Node 04a 独立核验器 v1.0
────────────────────────────────────────────────────────────────
目的：不复用 node04a_receptor_prep.py 的任何函数，仅从【磁盘上的文件】
重新推导 manifest 的每一项声明，再与 manifest 比对。
若两条独立代码路径得出同一结果，该结论才算被证实，而非被复述。

用法：
    python node04a_verify.py /path/to/Node04a_receptor_manifest.json
退出码：0 = 全部一致；1 = 存在不一致（不一致项逐条打印）。
"""
import hashlib
import json
import sys
from pathlib import Path


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def read_pdb_atoms(path, chain=None, records=("ATOM",)):
    """按 PDB 列格式解析，返回 [(resname, chain, seq, atomname, x, y, z, adtype)]。"""
    out = []
    for ln in Path(path).read_text().splitlines():
        if not ln.startswith(records):
            continue
        if chain and ln[21] != chain:
            continue
        try:
            seq = int(ln[22:26])
            xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        except ValueError:
            continue
        adtype = ln.split()[-1] if len(ln) > 70 else ""
        out.append((ln[17:20].strip(), ln[21], seq, ln[12:16].strip(), *xyz, adtype))
    return out


def verify(manifest_path):
    man = json.loads(Path(manifest_path).read_text())
    fails, notes = [], []

    def check(label, got, want):
        if got == want:
            notes.append(f"    ok   {label}: {got}")
        else:
            fails.append(f"    FAIL {label}: 重算={got} manifest={want}")

    print(f"manifest {manifest_path}")
    print(f"module {man['module']} v{man['version']}  dry_run={man['dry_run']}\n")

    tautomers = {}
    for pid, m in man["receptors"].items():
        print(f"── {pid} ({m['virus']}) ──")
        notes.clear()
        c = m.get("conversion") or {}
        if not c.get("executed"):
            fails.append(f"    FAIL {pid}: conversion 未执行")
            print("\n".join(fails)); continue

        # 1) 原始 PDB 与产出 PDBQT 的哈希
        for lbl, path, want in (("raw pdb sha256", m["input_pdb"], m["input_pdb_sha256"]),
                                ("pdbqt sha256", m["receptor_pdbqt"], c["pdbqt_sha256"])):
            if not Path(path).is_file():
                fails.append(f"    FAIL {lbl}: 文件不存在 {path}")
            else:
                check(lbl, sha256(path), want)

        # 2) PDBQT 原子普查（重新数，不看 manifest）
        pq = read_pdb_atoms(m["receptor_pdbqt"])
        n_polar = sum(1 for a in pq if a[-1] in ("HD", "H"))
        check("n_pdbqt_atoms", len(pq), c["n_pdbqt_atoms"])
        check("n_polar_h", n_polar, c["n_polar_h"])
        check("n_heavy_atoms_pdbqt", len(pq) - n_polar, c["n_heavy_atoms_pdbqt"])

        # 3) 清洗后 PDB 重原子（清洗阶段已删 H，故 ATOM 行即重原子）
        clean = read_pdb_atoms(m["clean_pdb"], chain=m["chain"])
        check("n_heavy_atoms_clean_pdb", len(clean), c["n_heavy_atoms_clean_pdb"])

        # 4) 守恒式独立重算：丢弃残基的重原子数从 clean PDB 现场统计
        per_res = {}
        for a in clean:
            per_res[f"{a[1]}:{a[2]}"] = per_res.get(f"{a[1]}:{a[2]}", 0) + 1
        n_drop = sum(per_res.get(k, 0) for k in c["dropped_incomplete_residues"])
        check("n_heavy_atoms_dropped", n_drop, c["n_heavy_atoms_dropped"])
        if len(clean) - n_drop != len(pq) - n_polar:
            fails.append(f"    FAIL 守恒: {len(clean)} - {n_drop} != {len(pq) - n_polar}")
        else:
            notes.append(f"    ok   守恒: {len(clean)} - {n_drop} = {len(pq) - n_polar}")

        # 5) 二联体：存在性 + 互变异构体（由 H 原子名判定）
        state = {}
        for tag in m["catalytic_dyad"]:
            resname, seq = tag[:3], int(tag[3:])
            names = {a[3] for a in pq if a[1] == m["chain"] and a[2] == seq}
            if not names:
                fails.append(f"    FAIL 二联体 {tag} 不在 PDBQT 中"); continue
            if resname == "HIS":
                d, e = "HD1" in names, "HE2" in names
                st = ("HIP_charged" if d and e else "HID_delta" if d
                      else "HIE_epsilon" if e else "no_polar_H")
                if not {"ND1", "NE2"} <= names:
                    fails.append(f"    FAIL {tag} 缺 ND1/NE2")
            else:
                st = "SH_neutral_thiol" if "HG" in names else "S-_thiolate"
                if "SG" not in names:
                    fails.append(f"    FAIL {tag} 缺 SG")
            state["HIS" if resname == "HIS" else "CYS"] = st
            check(f"二联体 {tag}", st, c["dyad_state"][tag]["state"])
        tautomers[pid] = state

        # 6) 盒中心：从原始 PDB 的共晶配体重新算质心
        het = m.get("box_center_extra", {}).get("cocrystal_ligand")
        if het and Path(m["input_pdb"]).is_file():
            lig = [a for a in read_pdb_atoms(m["input_pdb"], chain=m["chain"],
                                             records=("HETATM",))
                   if a[0] == het and not a[3].startswith("H")]
            if lig:
                ctr = [round(sum(a[4 + i] for a in lig) / len(lig), 3) for i in range(3)]
                check("box_center (质心重算)", ctr, m["box_center"])
                check("配体重原子数", len(lig),
                      m["box_center_extra"].get("n_ligand_heavy_atoms"))

        print("\n".join(notes))
        if fails:
            print("\n".join(fails))
        print()

    # 7) 跨受体二联体一致性
    uniq = {tuple(sorted(v.items())) for v in tautomers.values()}
    if len(uniq) == 1:
        print(f"跨受体二联体一致: {dict(tautomers[list(tautomers)[0]])}")
    else:
        fails.append(f"FAIL 跨受体二联体不一致: {tautomers}")

    print("\n" + ("=" * 60))
    if fails:
        print(f"❌ {len(fails)} 项不一致 —— manifest 未被独立证实，禁止进入 Node 04")
        return 1
    print("✅ 全部一致 —— manifest 已由独立代码路径证实，可作为 Node 04 锚点")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__); sys.exit(2)
    sys.exit(verify(sys.argv[1]))
