# -*- coding: utf-8 -*-
"""Node 04 → Node 05 · MD 选择集
    预览（只出计数与算力，【不】出分子名单）：
        python core/node04_md_selection.py preview
    定案（需显式给出百分比）：
        python core/node04_md_selection.py commit --pct 1.5

═══ 已登记的规则（在看到任何 BBS 之前写死）═══
1. 仅 pool_class == 'candidate'（碎片全部对接但不进聚类/选择）。
2. 电荷分层：预先登记的规则已触发（偏移 0.658 ≫ σ_seed 0.066），
   故在【各电荷分层内】分别取前 X%，比例与池构成一致（约 33:67）。
   理由：Vina 无静电项、无去溶剂化项，无法判定该偏移是真是假；
   MM-PBSA 二者皆有 —— 把电荷问题交给能回答它的那个阶段裁决。
3. 偏移稳健性：同时按 BBS_raw 与 BBS_centred 各取前 X%，取【并集】。
   理由：逐受体偏移（约 0.2 kcal/mol）解释了 63.4% 的限制性靶点分布，
   但对接无法证明它是假象。取并集使入选集不依赖该项校正。
4. N 由【算力预算】决定，不由分布决定 —— node04_bbs_gap_scan.py 已确认
   四个（定义 × 分层）组合均无 ≥3×σ_p95 的断层，分布连续。
   Methods 必须写作"在可用算力预算内选取各层前 X%"。

═══ 反循环保护 ═══
preview 模式【只输出计数与 GPU 成本，不输出任何分子标识】。
百分比必须在看到名单之前选定；commit 会把该次选择的时间与参数一并记录。
"""
import argparse, csv, json, statistics, sys
from pathlib import Path

R = Path("/mnt/hpc/home/wanghongmei/lang2")
TG = ["6W63", "3V3M", "4YLU"]
GPU_LO, GPU_HI, REPS = 24, 48, 3          # 每条 100 ns 轨迹的 GPU 时（V2 实测）


def load():
    rows = [r for r in csv.DictReader(
        open(R/"results/04_docking/node04_bbs_table.csv", encoding="utf-8"))
        if r.get("pool_class") == "candidate"]
    med = {t: statistics.median(float(r[f"mean_{t}"]) for r in rows) for t in TG}
    for r in rows:
        r["_raw"] = float(r["BBS"])
        r["_cen"] = max(float(r[f"mean_{t}"]) - med[t] for t in TG)
    return rows, med


def select(rows, pct):
    """返回 {stratum: {'raw':[...], 'cen':[...], 'union':[...]}}（分子名）。"""
    out = {}
    for cs in ("neutral", "charged"):
        sub = [r for r in rows if r["charge_stratum"] == cs]
        k = max(1, round(len(sub) * pct / 100))
        a = [r["ligand"] for r in sorted(sub, key=lambda r: r["_raw"])[:k]]
        b = [r["ligand"] for r in sorted(sub, key=lambda r: r["_cen"])[:k]]
        out[cs] = {"k_per_definition": k, "raw": a, "centred": b,
                   "union": sorted(set(a) | set(b)), "shared": sorted(set(a) & set(b))}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["preview", "commit"])
    ap.add_argument("--pct", type=float, default=None)
    a = ap.parse_args()
    rows, med = load()
    print(f"候选 {len(rows)}  (neutral {sum(1 for r in rows if r['charge_stratum']=='neutral')}"
          f" / charged {sum(1 for r in rows if r['charge_stratum']=='charged')})")

    if a.mode == "preview":
        print("\n【预览：只出计数与算力，不出分子名单 —— 防止先看赢家再定百分比】")
        print(f"{'各层前':>7}{'中性':>7}{'带电':>7}{'单定义':>8}{'并集':>7}"
              f"{'轨迹':>7}{'GPU 时':>16}")
        for p in (0.5, 1.0, 1.5, 2.0, 3.0):
            s = select(rows, p)
            per = s["neutral"]["k_per_definition"] + s["charged"]["k_per_definition"]
            uni = len(set(s["neutral"]["union"]) | set(s["charged"]["union"]))
            tr = uni * REPS
            print(f"{p:>6.1f}%{s['neutral']['k_per_definition']:>7}"
                  f"{s['charged']['k_per_definition']:>7}{per:>8}{uni:>7}{tr:>7}"
                  f"{tr*GPU_LO:>8,}–{tr*GPU_HI:<7,}")
        print("\n选定百分比后：commit --pct <值>")
        return 0

    if a.pct is None:
        raise SystemExit("commit 必须显式给出 --pct（该值须在看到名单之前选定）")
    s = select(rows, a.pct)
    uni = sorted(set(s["neutral"]["union"]) | set(s["charged"]["union"]))
    by = {r["ligand"]: r for r in rows}
    print(f"\n选中 {len(uni)} 个分子（各层前 {a.pct}%，raw ∪ centred）")
    for cs in ("neutral", "charged"):
        d = s[cs]
        print(f"  {cs:<8} 每定义 {d['k_per_definition']:>2} | 共有 {len(d['shared']):>2}"
              f" | 并集 {len(d['union']):>2} | 仅raw {len(set(d['raw'])-set(d['centred'])):>2}"
              f" | 仅centred {len(set(d['centred'])-set(d['raw'])):>2}")
    print(f"\n{'分子':<18}{'BBS_raw':>9}{'BBS_cen':>9}{'限制靶点':>10}{'电荷':>6}{'来源':>12}")
    for L in sorted(uni, key=lambda x: by[x]["_raw"]):
        r = by[L]
        src = ("both" if L in s[r["charge_stratum"]]["shared"] else
               "raw" if L in s[r["charge_stratum"]]["raw"] else "centred")
        print(f"{L:<18}{r['_raw']:>9.3f}{r['_cen']:>9.3f}"
              f"{r['BBS_limiting_target']:>10}{r['formal_charge']:>6}{src:>12}")

    tr = len(uni) * REPS
    man = {"module": "node04_md_selection", "created_at_utc": __import__("datetime")
           .datetime.now(__import__("datetime").timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "rule": {"stratum": "equal percentile within charge strata (~33:67, matches pool)",
                    "pct_per_stratum": a.pct,
                    "definitions": "union of top-X% under BBS_raw and BBS_centred",
                    "N_basis": "COMPUTE BUDGET — gap scan found no break ≥3×sigma_p95; "
                               "distribution is continuous. Must be reported as a budget cut.",
                    "pre_registered": True},
           "baseline_medians": {t: round(med[t], 4) for t in TG},
           "n_selected": len(uni), "n_trajectories": tr,
           "gpu_hours_estimate": [tr*GPU_LO, tr*GPU_HI],
           "per_stratum": {cs: {k: v for k, v in s[cs].items()} for cs in s},
           "selected": [{"ligand": L, "BBS_raw": round(by[L]["_raw"], 4),
                         "BBS_centred": round(by[L]["_cen"], 4),
                         "limiting_target": by[L]["BBS_limiting_target"],
                         "formal_charge": by[L]["formal_charge"],
                         "charge_stratum": by[L]["charge_stratum"]} for L in uni]}
    p = R/"results/04_docking/Node04_md_selection.json"
    p.write_text(json.dumps(man, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n{tr} 条轨迹 | {tr*GPU_LO:,}–{tr*GPU_HI:,} GPU-h\n✓ {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
