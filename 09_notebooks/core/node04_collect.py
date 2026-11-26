# -*- coding: utf-8 -*-
"""Node 04 汇总 · 26,475 次运行 → BBS 表 + manifest
    cd /mnt/hpc/home/wanghongmei/lang2
    /mnt/hpc/.../md_env/bin/python core/node04_collect.py

顺序刻意如此：先断言完整，再计算；先按电荷分层，再决定门槛规则。
门槛规则【预先登记】—— 在看到哪些分子被选中【之前】就已写死，
以免"看到结果再定规则"把选择变成循环论证（与 pool_class 的反循环原则一致）。
"""
import csv
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

R = Path("/mnt/hpc/home/wanghongmei/lang2")
sys.path.insert(0, str(R))
from core import node04_docking as N4          # noqa: E402

OUT = R / "results" / "04_docking"
VINA = Path("/mnt/hpc/home/wanghongmei/apps/vina/vina_1.2.7_linux_x86_64")

# ── 预先登记的门槛规则（在看到任何 BBS 之前定死）──────────────────
THRESHOLD_RULE = {
    "id": "charge_stratified_if_shift_exceeds_noise",
    "statement": ("若 |median BBS(charged) − median BBS(neutral)| > median σ_seed，"
                  "则在各电荷分层内分别设 MD 门槛；否则在合并分布上设单一门槛。"),
    "rationale": ("Vina 无静电项、无去溶剂化惩罚 → 带电分子可能被系统性高估。"
                  "高估的带电分子进 MD 后会被 MM-PBSA 降级（可自我纠正）；"
                  "但被其挤到门槛下的优质中性候选永不进 MD（静默假阴性，不可挽回）。"
                  "故规则必须在看到结果前定死。"),
    "pre_registered": True,
}


def main():
    # ── 1 · 上游闸门（与开跑前同一套断言）────────────────────────
    gate, recs, m4a = N4.verify_upstream(
        R / "results/03_ligand_pdbqt/Node03_manifest.json",
        R / "results/03_ligand_pdbqt/identity_manifests_node03.jsonl",
        OUT / "receptors/Node04a_receptor_manifest.json")
    import pandas as pd
    pool_class, charge, pc_stats = N4.derive_pool_class(
        R / "results/02_master_molecules/identity_manifests.jsonl",
        node00_csv=R / "results/00_scope/all_candidates_v3_start.csv", pandas_mod=pd)
    print("闸门:", gate["assertion"])

    wm = OUT / "node04_work_manifest.jsonl"
    rows = [json.loads(l) for l in wm.read_text(encoding="utf-8").splitlines() if l.strip()]
    print(f"工作清单 {len(rows)} 行 | sha256 {N4.file_sha256(wm)[:16]}…")

    # ── 2 · 汇总（不完整即 raise，绝不在残缺数据上算 BBS）─────────
    bbs = N4.collect(rows, pool_class=pool_class, charge=charge, seeds=N4.SEEDS)
    print(f"BBS 表: {len(bbs)} 个配体")

    # ── 3 · 分层报告（设门槛【之前】必须先看这张表）───────────────
    strata = N4.stratify(bbs)
    print("\npool_class × charge 分层:")
    for s in strata:
        print(f"  {s['pool_class']:<10} {s['charge']:<8} n={s['n']:<5} "
              f"中位 {s['BBS_median']:>7} 均值 {s['BBS_mean']:>7} 最好 {s['BBS_best']:>7}")

    # ── 3b · 限制性靶点分布（诊断：区分真实广谱难度 vs 逐靶点系统偏移）──
    from collections import Counter
    lim = Counter(r["BBS_limiting_target"] for r in bbs if r.get("pool_class") == "candidate")
    n_c = sum(lim.values())
    print("\nBBS 限制性靶点分布（candidate）:")
    for t, n in lim.most_common():
        print(f"  {t}: {n:>5} ({n/n_c*100:>5.1f}%)")
    dominant = lim.most_common(1)[0]
    if dominant[1] / n_c > 0.85:
        print(f"  ⚠ {dominant[0]} 占 {dominant[1]/n_c*100:.0f}% —— 这更像【逐靶点系统偏移】"
              "而非真实广谱难度。BBS 若几乎总由同一靶点决定，它就退化为该单靶点分数。")
        print("     诊断参照：自对接对照的天然物参考分（X77 -8.52 / 0EN -7.68 / R30 -7.98，"
              "极差 0.84 kcal/mol；但三者是不同配体，只能作提示不能作校准）。")
    else:
        print(f"  混合分布（最高 {dominant[0]} {dominant[1]/n_c*100:.0f}%）—— BBS 确实在跨靶点取舍。")

    # ── 4 · 施加预先登记的规则 ───────────────────────────────────
    cand = [r for r in bbs if r.get("pool_class") == "candidate"]
    neu = [r["BBS"] for r in cand if r.get("charge_stratum") == "neutral"]
    chg = [r["BBS"] for r in cand if r.get("charge_stratum") == "charged"]
    sig = [r["sigma_seed_max"] for r in cand]
    shift = abs(statistics.median(chg) - statistics.median(neu)) if (neu and chg) else 0.0
    noise = statistics.median(sig) if sig else 0.0
    stratify_needed = shift > noise
    decision = {
        "median_BBS_neutral": round(statistics.median(neu), 4) if neu else None,
        "median_BBS_charged": round(statistics.median(chg), 4) if chg else None,
        "shift_kcal_per_mol": round(shift, 4),
        "median_sigma_seed": round(noise, 4),
        "rule_outcome": "threshold WITHIN charge strata" if stratify_needed
                        else "single pooled threshold",
        "n_candidate_neutral": len(neu), "n_candidate_charged": len(chg),
    }
    print(f"\n预先登记的门槛规则 → {decision['rule_outcome']}")
    print(f"  偏移 {shift:.3f} vs 噪声 σ_seed {noise:.3f} kcal/mol")
    if stratify_needed:
        print("  ⚠ 带电与中性分布偏移超过对接噪声：合并门槛会系统性排挤中性候选。")

    # ── 5 · 落盘：BBS 表 + manifest ─────────────────────────────
    csv_path = OUT / "node04_bbs_table.csv"
    cols = list(bbs[0].keys())
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(bbs)

    manifest = {
        "module": "node04_docking", "version": N4.NODE04_VERSION,
        "created_at_utc": N4.utc_now(),
        "upstream": gate,
        "pool_class": pc_stats,
        "docking": {
            "engine": "AutoDock Vina", "vina_path": str(VINA),
            "vina_sha256": N4.file_sha256(VINA),
            "exhaustiveness": N4.EXHAUSTIVENESS, "num_modes": N4.NUM_MODES,
            "cpu_per_process": N4.VINA_CPU, "seeds": list(N4.SEEDS),
            "n_runs": len(rows), "n_ligands": len(bbs),
            "work_manifest_sha256": N4.file_sha256(wm),
        },
        "bbs": {"rule": N4.BBS_RULE, "note": N4.BBS_NOTE,
                "aggregate": "per (ligand,target) = mean of top-mode affinity over seeds; "
                             "sigma = stdev of the same quantity across seeds"},
        "threshold_rule": THRESHOLD_RULE,
        "threshold_decision": decision,
        "strata": strata,
        "bbs_limiting_target_distribution": dict(lim),
        "outputs": {"bbs_table": str(csv_path), "bbs_table_sha256": N4.file_sha256(csv_path)},
    }
    mp = OUT / "Node04_manifest.json"
    mp.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n✓ {csv_path}\n✓ {mp}")
    print("下一步：core/node04_verify.py 独立复核，再谈 MD 门槛。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
