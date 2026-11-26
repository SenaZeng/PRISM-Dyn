# -*- coding: utf-8 -*-
"""
V3_clean/core/node04_docking.py
===============================
Node 04 — AutoDock Vina 对接（1765 配体 × 3 靶点 × 5 种子 = 26,475 次）。

设计原则（承接 Node 00–04a）：
  * 上游闸门 fail-closed：开跑前重算 Node 03 与 Node 04a 的全部哈希，任一不符即拒绝。
  * 幂等可续跑：输出已存在且哈希吻合的任务跳过；中断后是续跑而非重启。
  * 工具身份入链：vina 绝对路径 + sha256 + --version 记入 manifest。
  * `--cpu` 视作计算身份的一部分（已实测 cpu=1/4 输出逐位相同，但仍冻结并记录）。
  * BBS 方向 = 跨靶点【取最不利】：Vina 亲和力为负，故 BBS = max(三靶均值)。
    只有最弱靶点仍好的化合物才算广谱。以具名常量表达，不埋在比较符里。
  * pool_class 由 Node 02 冻结的 heavy_atom_count 重新推导（Node 00 CSV 曾被原地覆盖）。

Slurm 环境陷阱（已实测）：srun 进入的是 base 而非 md_env。
  因此 worker 一律使用【绝对路径】，不依赖 PATH 或 conda activate。
"""

import hashlib
import json
import os
import statistics
import subprocess
from datetime import datetime, timezone
from pathlib import Path

NODE04_VERSION = "0.1.0"

# ── 冻结的对接参数（改动即须升版本号）────────────────────────────
EXHAUSTIVENESS = 32
NUM_MODES = 9
VINA_CPU = 1                      # 每进程 1 核：核效率最高（4 核加速比仅 2.44×）
SEEDS = (42, 43, 44, 45, 46)      # 固定、记录、对每个配体-靶点对相同
BOX_SIZE = (27.0, 27.0, 27.0)

# ── BBS 方向：具名常量，使选择在 manifest 中可见 ──────────────────
BBS_RULE = "worst_case_across_targets"
BBS_NOTE = ("Vina affinity is negative; worst-case = max(per-target mean). "
            "A compound is broad-spectrum only if its WEAKEST protease is still good.")

# ── pool_class 阈值（与 Node 00 backfill 一致）────────────────────
HEAVY_ATOM_MIN = 12


def utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


# ══════════════════════════ 1 · 上游闸门 ══════════════════════════
def verify_upstream(node03_manifest, node03_jsonl, node04a_manifest,
                    check_every_ligand=True):
    """开跑前的 fail-closed 闸门。返回可入 manifest 的 gate 记录。

    重算：Node 03 manifest、Node 04a manifest、三个受体 PDBQT、
    以及（默认）全部 1765 个配体 PDBQT 的 sha256。任一不符即 raise。
    """
    n3m = Path(node03_manifest); n3j = Path(node03_jsonl); n4a = Path(node04a_manifest)
    for p in (n3m, n3j, n4a):
        if not p.is_file():
            raise AssertionError(f"上游产物缺失：{p}")

    m3 = json.loads(n3m.read_text(encoding="utf-8"))
    m4 = json.loads(n4a.read_text(encoding="utf-8"))
    recs = [json.loads(l) for l in n3j.read_text(encoding="utf-8").splitlines() if l.strip()]

    if len(recs) != m3["stats"]["n_built"]:
        raise AssertionError(
            f"Node03 行数 {len(recs)} != manifest n_built {m3['stats']['n_built']}")

    # 受体：三个 PDBQT 必须与 Node 04a 记录一致
    bad_rec = []
    for pid, m in m4["receptors"].items():
        c = m.get("conversion") or {}
        if not c.get("executed"):
            raise AssertionError(f"Node04a {pid}: conversion 未执行")
        p = Path(m["receptor_pdbqt"])
        if not p.is_file() or file_sha256(p) != c["pdbqt_sha256"]:
            bad_rec.append(pid)
    if bad_rec:
        raise AssertionError(f"受体 PDBQT 哈希漂移或缺失：{bad_rec}")

    # 配体：全量重算（~1765 次，数秒）
    bad_lig, missing = [], []
    if check_every_ligand:
        for r in recs:
            p = Path(r["pdbqt_path"])
            if not p.is_file():
                missing.append(r["name"])
            elif file_sha256(p) != r["pdbqt_sha256"]:
                bad_lig.append(r["name"])
    if missing or bad_lig:
        raise AssertionError(
            f"配体 PDBQT 缺失 {len(missing)} / 哈希漂移 {len(bad_lig)}："
            f"{(missing + bad_lig)[:10]}")

    return {
        "node03_manifest_sha256": file_sha256(n3m),
        "node03_jsonl_sha256": file_sha256(n3j),
        "node04a_manifest_sha256": file_sha256(n4a),
        "n_ligands_verified": len(recs) if check_every_ligand else 0,
        "n_receptors_verified": len(m4["receptors"]),
        "node03_version": m3["version"],
        "node04a_version": m4["version"],
        "assertion": "all upstream hashes re-derived and matched",
    }, recs, m4


# ══════════════════════════ 2 · pool_class ═══════════════════════
def derive_pool_class(node02_jsonl, heavy_min=HEAVY_ATOM_MIN, node00_csv=None,
                      pandas_mod=None):
    """由 Node 02 冻结记录重新推导 pool_class（首要来源）。
    若给出 Node 00 CSV 则逐分子交叉核对（次要来源），不符即 raise。"""
    recs = [json.loads(l) for l in
            Path(node02_jsonl).read_text(encoding="utf-8").splitlines() if l.strip()]
    derived = {r["name"]: ("fragment" if r["heavy_atom_count"] < heavy_min else "candidate")
               for r in recs}
    charge = {r["name"]: int(r["frozen_net_charge"]) for r in recs}   # 上游为字符串，此处强转

    xcheck = None
    if node00_csv and pandas_mod is not None and Path(node00_csv).is_file():
        df = pandas_mod.read_csv(node00_csv)
        if "pool_class" in df.columns and "name" in df.columns:
            ref = dict(zip(df["name"], df["pool_class"]))
            parent = {r["name"]: r.get("parent_id") for r in recs}
            mism = [n for n, c in derived.items()
                    if parent.get(n) in ref and ref[parent[n]] != c]
            if mism:
                raise AssertionError(
                    f"pool_class 与 Node 00 表不符（{len(mism)} 个）：{mism[:10]}")
            xcheck = {"source": str(node00_csv), "n_cross_checked": len(derived),
                      "n_mismatch": 0}
    n_frag = sum(1 for v in derived.values() if v == "fragment")
    return derived, charge, {"heavy_atom_min": heavy_min, "n_fragment": n_frag,
                             "n_candidate": len(derived) - n_frag,
                             "primary_source": "node02 frozen heavy_atom_count",
                             "cross_check_vs_node00": xcheck}


# ══════════════════════ 3 · 工作清单 + 分块 ═══════════════════════
def build_work_manifest(ligand_recs, node04a_manifest_obj, out_root, seeds=SEEDS):
    """生成 (配体 × 靶点 × 种子) 的全部行。每行自带预期输出路径。"""
    out_root = Path(out_root)
    rows = []
    for r in sorted(ligand_recs, key=lambda x: x["name"]):
        for pid, m in sorted(node04a_manifest_obj["receptors"].items()):
            for s in seeds:
                stem = f"{r['name']}__{pid}__s{s}"
                rows.append({
                    "ligand": r["name"],
                    "ligand_pdbqt": r["pdbqt_path"],
                    "ligand_pdbqt_sha256": r["pdbqt_sha256"],
                    "target": pid,
                    "receptor_pdbqt": m["receptor_pdbqt"],
                    "receptor_pdbqt_sha256": m["conversion"]["pdbqt_sha256"],
                    "box_center": m["box_center"],
                    "box_size": list(m.get("box_size", BOX_SIZE)),
                    "seed": s,
                    "out_pdbqt": str(out_root / "poses" / pid / f"{stem}.pdbqt"),
                    "out_json": str(out_root / "runs" / pid / f"{stem}.json"),
                })
    for i, row in enumerate(rows):
        row["row_index"] = i
    return rows


def chunk_bounds(n_rows, n_tasks):
    """把 n_rows 均分给 n_tasks 个 array 任务，返回 [(lo, hi), ...]（左闭右开）。"""
    if n_tasks < 1:
        raise ValueError("n_tasks must be >= 1")
    base, extra = divmod(n_rows, n_tasks)
    out, lo = [], 0
    for t in range(n_tasks):
        hi = lo + base + (1 if t < extra else 0)
        out.append((lo, hi))
        lo = hi
    assert lo == n_rows
    return out


# ══════════════════════ 4 · 单次运行（幂等）═══════════════════════
def parse_vina_poses(pdbqt_text):
    """从 Vina 输出 PDBQT 提取各 mode 亲和力（kcal/mol），按文件顺序。"""
    out = []
    for ln in pdbqt_text.splitlines():
        if ln.startswith("REMARK VINA RESULT:"):
            try:
                out.append(float(ln.split(":", 1)[1].split()[0]))
            except (IndexError, ValueError):
                pass
    return out


def run_one(row, vina_exe, exhaustiveness=EXHAUSTIVENESS, num_modes=NUM_MODES,
            cpu=VINA_CPU, dry_run=False, timeout_s=3600):
    """跑一次对接。幂等：若 out_json 已存在且其记录的输入哈希吻合，则跳过。"""
    out_json = Path(row["out_json"]); out_pdbqt = Path(row["out_pdbqt"])

    if out_json.is_file():
        try:
            prev = json.loads(out_json.read_text(encoding="utf-8"))
            same_inputs = (prev.get("ligand_pdbqt_sha256") == row["ligand_pdbqt_sha256"]
                           and prev.get("receptor_pdbqt_sha256") == row["receptor_pdbqt_sha256"]
                           and prev.get("seed") == row["seed"]
                           and prev.get("exhaustiveness") == exhaustiveness
                           and prev.get("cpu") == cpu)
            if same_inputs and prev.get("status") == "ok" and Path(prev["out_pdbqt"]).is_file() \
               and file_sha256(prev["out_pdbqt"]) == prev.get("pose_sha256"):
                prev["skipped_existing"] = True
                return prev
        except Exception:
            pass          # 记录损坏 → 重跑

    cx, cy, cz = row["box_center"]; sx, sy, sz = row["box_size"]
    cmd = [str(vina_exe),
           "--receptor", row["receptor_pdbqt"], "--ligand", row["ligand_pdbqt"],
           "--center_x", f"{cx}", "--center_y", f"{cy}", "--center_z", f"{cz}",
           "--size_x", f"{sx}", "--size_y", f"{sy}", "--size_z", f"{sz}",
           "--exhaustiveness", str(exhaustiveness), "--num_modes", str(num_modes),
           "--seed", str(row["seed"]), "--cpu", str(cpu),
           "--out", str(out_pdbqt)]
    if dry_run:
        return {**row, "status": "dry_run", "cmd": " ".join(cmd)}

    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    t0 = datetime.now(timezone.utc)
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s)
        rc, err = r.returncode, (r.stderr or "")[:400]
    except subprocess.TimeoutExpired:
        rc, err = -9, f"timeout after {timeout_s}s"
    dt = (datetime.now(timezone.utc) - t0).total_seconds()

    rec = {k: row[k] for k in ("ligand", "target", "seed", "row_index",
                               "ligand_pdbqt_sha256", "receptor_pdbqt_sha256")}
    rec.update({"cmd": " ".join(cmd), "exhaustiveness": exhaustiveness,
                "num_modes": num_modes, "cpu": cpu, "wall_seconds": round(dt, 2),
                "out_pdbqt": str(out_pdbqt), "created_at_utc": utc_now()})

    if rc != 0 or not out_pdbqt.is_file():
        rec.update({"status": "fail", "returncode": rc, "error": err,
                    "affinities": None, "affinity_best": None})
    else:
        aff = parse_vina_poses(out_pdbqt.read_text())
        if not aff:
            rec.update({"status": "fail", "returncode": rc,
                        "error": "no REMARK VINA RESULT in output",
                        "affinities": None, "affinity_best": None})
        else:
            rec.update({"status": "ok", "returncode": rc,
                        "affinities": aff, "affinity_best": aff[0],
                        "n_modes_returned": len(aff),
                        "pose_sha256": file_sha256(out_pdbqt)})
    out_json.write_text(json.dumps(rec, ensure_ascii=False, indent=2), encoding="utf-8")
    return rec


def run_slice(rows, lo, hi, vina_exe, **kw):
    """array 任务入口：处理 rows[lo:hi]。"""
    done = []
    for row in rows[lo:hi]:
        done.append(run_one(row, vina_exe, **kw))
    return done


# ══════════════════════════ 5 · 汇总 ═════════════════════════════
def collect(rows, pool_class=None, charge=None, seeds=SEEDS):
    """读全部 out_json，先断言完整，再算 per-(ligand,target) 均值/σ 与 BBS。"""
    per = {}
    missing, failed = [], []
    for row in rows:
        p = Path(row["out_json"])
        if not p.is_file():
            missing.append(row["out_json"]); continue
        rec = json.loads(p.read_text(encoding="utf-8"))
        if rec.get("status") != "ok" or rec.get("affinity_best") is None:
            failed.append((rec.get("ligand"), rec.get("target"), rec.get("seed")))
            continue
        per.setdefault((rec["ligand"], rec["target"]), []).append(
            (rec["seed"], rec["affinity_best"]))

    if missing or failed:
        raise AssertionError(
            f"汇总拒绝：缺失 {len(missing)} / 失败 {len(failed)}。"
            f"示例 {(missing[:3] + [str(x) for x in failed[:3]])}。"
            "先补跑（幂等，重提交即可），不得在不完整数据上计算 BBS。")

    lig_tgt = {}
    for (lig, tgt), vals in per.items():
        seen = sorted(s for s, _ in vals)
        if seen != sorted(seeds):
            raise AssertionError(f"{lig}/{tgt} 种子集合 {seen} != {sorted(seeds)}")
        a = [v for _, v in sorted(vals)]
        lig_tgt[(lig, tgt)] = {
            "mean": round(statistics.fmean(a), 4),
            "sigma_seed": round(statistics.stdev(a), 4) if len(a) > 1 else 0.0,
            "min": min(a), "max": max(a), "n_seeds": len(a),
        }

    ligands = sorted({l for l, _ in lig_tgt})
    targets = sorted({t for _, t in lig_tgt})
    out = []
    for lig in ligands:
        means = {t: lig_tgt[(lig, t)]["mean"] for t in targets if (lig, t) in lig_tgt}
        if len(means) != len(targets):
            raise AssertionError(f"{lig} 缺靶点：{set(targets) - set(means)}")
        bbs = max(means.values())                      # ★ 取最不利靶点
        row = {"ligand": lig, "BBS": round(bbs, 4), "BBS_rule": BBS_RULE,
               "BBS_limiting_target": max(means, key=means.get)}
        for t in targets:
            row[f"mean_{t}"] = lig_tgt[(lig, t)]["mean"]
            row[f"sigma_{t}"] = lig_tgt[(lig, t)]["sigma_seed"]
        row["sigma_seed_max"] = round(max(lig_tgt[(lig, t)]["sigma_seed"]
                                          for t in targets), 4)
        if pool_class:
            row["pool_class"] = pool_class.get(lig)
        if charge:
            row["formal_charge"] = charge.get(lig)
            row["charge_stratum"] = "neutral" if charge.get(lig) == 0 else "charged"
        out.append(row)
    return out


def stratify(bbs_rows):
    """按 pool_class × charge 分层报告 BBS —— 阈值设定前必须先看这张表。"""
    groups = {}
    for r in bbs_rows:
        k = (r.get("pool_class", "?"), r.get("charge_stratum", "?"))
        groups.setdefault(k, []).append(r["BBS"])
    out = []
    for (pc, cs), v in sorted(groups.items()):
        v = sorted(v)
        out.append({"pool_class": pc, "charge": cs, "n": len(v),
                    "BBS_mean": round(statistics.fmean(v), 3),
                    "BBS_median": round(statistics.median(v), 3),
                    "BBS_best": v[0], "BBS_worst": v[-1]})
    return out
