# -*- coding: utf-8 -*-
"""
V3_clean/core/node04_worker.py
==============================
Node 04 的 Slurm array 任务入口（薄封装）。

  python -m core.node04_worker --work-manifest <jsonl> --task-id N --n-tasks M \
         --vina <abs_path> --workers K

职责仅三条：
  1. 读工作清单，按 chunk_bounds 取出本任务负责的连续切片；
  2. 用 K 个进程并行执行该切片（每个 Vina 进程 --cpu 1）；
  3. 写一份任务级完成记录，供汇总阶段核对。

不做任何科学判断。所有对接语义都在 node04_docking.py 里。
Slurm 不继承 conda 环境 —— 故本模块只用标准库，vina 走绝对路径。
"""

import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from core import node04_docking as N4          # noqa: E402

_CFG = {}


def _init(cfg):
    global _CFG
    _CFG = cfg


def _do(row):
    """子进程执行体。异常不得杀死整个任务 —— 记录为 fail，由汇总阶段拦截。"""
    try:
        r = N4.run_one(row, _CFG["vina"],
                       exhaustiveness=_CFG["exhaustiveness"],
                       num_modes=_CFG["num_modes"],
                       cpu=_CFG["cpu"], timeout_s=_CFG["timeout_s"])
        return {"ligand": row["ligand"], "target": row["target"], "seed": row["seed"],
                "status": r.get("status"), "skipped": bool(r.get("skipped_existing")),
                "wall_seconds": r.get("wall_seconds")}
    except Exception as e:                       # noqa: BLE001
        return {"ligand": row["ligand"], "target": row["target"], "seed": row["seed"],
                "status": "error", "error": f"{type(e).__name__}: {e}"}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--work-manifest", required=True)
    ap.add_argument("--task-id", type=int, required=True)
    ap.add_argument("--n-tasks", type=int, required=True)
    ap.add_argument("--vina", required=True)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--exhaustiveness", type=int, default=N4.EXHAUSTIVENESS)
    ap.add_argument("--num-modes", type=int, default=N4.NUM_MODES)
    ap.add_argument("--cpu", type=int, default=N4.VINA_CPU)
    ap.add_argument("--timeout-s", type=int, default=3600)
    a = ap.parse_args(argv)

    vina = Path(a.vina)
    if not (vina.is_file() and os.access(vina, os.X_OK)):
        raise SystemExit(f"vina 不可执行（须绝对路径，Slurm 不继承 conda）：{vina}")

    wm = Path(a.work_manifest)
    rows = [json.loads(l) for l in wm.read_text(encoding="utf-8").splitlines() if l.strip()]
    lo, hi = N4.chunk_bounds(len(rows), a.n_tasks)[a.task_id]
    my = rows[lo:hi]
    print(f"[n04] task {a.task_id}/{a.n_tasks}  rows[{lo}:{hi}]  n={len(my)}  "
          f"workers={a.workers}  cpu={a.cpu}", flush=True)

    cfg = {"vina": str(vina), "exhaustiveness": a.exhaustiveness,
           "num_modes": a.num_modes, "cpu": a.cpu, "timeout_s": a.timeout_s}
    with ProcessPoolExecutor(max_workers=a.workers,
                             initializer=_init, initargs=(cfg,)) as ex:
        res = list(ex.map(_do, my, chunksize=1))

    n_ok = sum(1 for r in res if r["status"] == "ok")
    n_skip = sum(1 for r in res if r.get("skipped"))
    bad = [r for r in res if r["status"] not in ("ok",)]
    # v1.1 修复：任务记录文件名必须按【工作清单 + 分块数】命名空间隔离。
    #   否则 (task 0 / n=700 冒烟) 与 (task 0 / n=32 试点) 与 (task 0 / n=353 全量)
    #   会写到同一个 task_00000.json，互相覆盖 —— 对接结果不受影响（run JSON 按
    #   ligand__target__seed 命名），但任务级审计记录会丢失。
    done = (wm.parent / "tasks" / f"{wm.stem}__n{a.n_tasks}"
            / f"task_{a.task_id:05d}.json")
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text(json.dumps(
        {"task_id": a.task_id, "n_tasks": a.n_tasks, "slice": [lo, hi],
         "n_rows": len(my), "n_ok": n_ok, "n_skipped_existing": n_skip,
         "n_bad": len(bad), "bad": bad[:20],
         "work_manifest_sha256": N4.file_sha256(wm),
         "vina": str(vina), "vina_sha256": N4.file_sha256(vina),
         "exhaustiveness": a.exhaustiveness, "num_modes": a.num_modes, "cpu": a.cpu,
         "created_at_utc": N4.utc_now()},
        ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"[n04] task {a.task_id}: ok={n_ok} skipped={n_skip} bad={len(bad)}", flush=True)
    # 非零退出 → sacct 里立刻可见；幂等，重提交即补跑。
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
