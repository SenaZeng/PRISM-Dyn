#!/usr/bin/env python3
"""Node 04 独立核验器 v1.0
────────────────────────────────────────────────────────────────
与 node04a_verify.py 同一形制：【不复用】node04_docking.py 的任何函数，
仅从磁盘上的 run JSON 与姿势文件重新推导 BBS 表的每一项，再与之比对。
两条独立代码路径吻合才算证实；一条路径自我复述不算。

    python core/node04_verify.py results/04_docking/Node04_manifest.json
退出码 0 = 全部一致；1 = 存在不一致。
"""
import csv
import hashlib
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def parse_affinities(pdbqt_path):
    """从姿势 PDBQT 重新解析亲和力 —— 不看 run JSON 里记的数。"""
    out = []
    for ln in Path(pdbqt_path).read_text().splitlines():
        if ln.startswith("REMARK VINA RESULT:"):
            out.append(float(ln.split(":", 1)[1].split()[0]))
    return out


def verify(manifest_path):
    man = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    d = man["docking"]
    fails, notes = [], []
    print(f"module {man['module']} v{man['version']}  {man['created_at_utc']}")

    # 1) 工具身份与工作清单
    v = Path(d["vina_path"])
    if not v.is_file() or sha256(v) != d["vina_sha256"]:
        fails.append(f"FAIL vina 二进制哈希不符或缺失: {v}")
    else:
        notes.append(f"ok   vina {d['vina_sha256'][:16]}…")

    root = Path(manifest_path).resolve().parent
    wm = root / "node04_work_manifest.jsonl"
    if not wm.is_file() or sha256(wm) != d["work_manifest_sha256"]:
        fails.append("FAIL 工作清单哈希不符或缺失")
        print("\n".join(notes + fails)); return 1
    rows = [json.loads(l) for l in wm.read_text(encoding="utf-8").splitlines() if l.strip()]
    if len(rows) != d["n_runs"]:
        fails.append(f"FAIL 清单行数 {len(rows)} != manifest n_runs {d['n_runs']}")
    else:
        notes.append(f"ok   工作清单 {len(rows)} 行")

    # 2) 逐次运行：姿势文件哈希 + 亲和力重新解析
    per = defaultdict(dict)
    n_missing = n_hash = n_aff = 0
    seeds = set(d["seeds"])
    for row in rows:
        jp, pp = Path(row["out_json"]), Path(row["out_pdbqt"])
        if not jp.is_file() or not pp.is_file():
            n_missing += 1; continue
        rec = json.loads(jp.read_text(encoding="utf-8"))
        if rec.get("pose_sha256") and sha256(pp) != rec["pose_sha256"]:
            n_hash += 1; continue
        aff = parse_affinities(pp)
        if not aff or abs(aff[0] - rec["affinity_best"]) > 1e-9:
            n_aff += 1; continue
        if len(aff) > d["num_modes"]:
            fails.append(f"FAIL {pp.name}: {len(aff)} 个构象 > num_modes {d['num_modes']}")
        per[(rec["ligand"], rec["target"])][rec["seed"]] = aff[0]
    for lbl, n in (("缺失文件", n_missing), ("姿势哈希漂移", n_hash), ("亲和力不符", n_aff)):
        (notes if n == 0 else fails).append(f"{'ok  ' if n == 0 else 'FAIL'} {lbl}: {n}")

    # 3) 独立重算 BBS，与 CSV 逐配体比对
    tbl = {r["ligand"]: r for r in csv.DictReader(
        open(man["outputs"]["bbs_table"], encoding="utf-8"))}
    if sha256(man["outputs"]["bbs_table"]) != man["outputs"]["bbs_table_sha256"]:
        fails.append("FAIL BBS 表哈希漂移")
    means = defaultdict(dict)
    for (lig, tgt), dd in per.items():
        if set(dd) != seeds:
            fails.append(f"FAIL {lig}/{tgt} 种子集合不完整: {sorted(dd)}")
            continue
        means[lig][tgt] = statistics.fmean(dd[s] for s in sorted(seeds))
    n_bbs_bad = 0
    for lig, mt in means.items():
        if lig not in tbl:
            fails.append(f"FAIL {lig} 不在 BBS 表中"); continue
        if abs(max(mt.values()) - float(tbl[lig]["BBS"])) > 1e-3:   # ★ 最不利靶点
            n_bbs_bad += 1
    (notes if n_bbs_bad == 0 else fails).append(
        f"{'ok  ' if n_bbs_bad == 0 else 'FAIL'} BBS 独立重算不符: {n_bbs_bad} / {len(means)}")
    if len(means) != len(tbl):
        fails.append(f"FAIL 配体数 重算 {len(means)} != 表 {len(tbl)}")
    else:
        notes.append(f"ok   配体数 {len(tbl)}")

    print("\n".join("    " + x for x in notes))
    if fails:
        print("\n".join("    " + x for x in fails))
    print("\n" + "=" * 62)
    if fails:
        print(f"❌ {len(fails)} 项不一致 —— Node 04 未被独立证实，禁止据此设 MD 门槛")
        return 1
    print("✅ 全部一致 —— BBS 表已由独立代码路径证实")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__); sys.exit(2)
    sys.exit(verify(sys.argv[1]))
