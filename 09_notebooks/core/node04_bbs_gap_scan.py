# -*- coding: utf-8 -*-
"""Node 04 · BBS 分布断层扫描（在设定 N 之前运行，程序【预先登记】如下）
    /mnt/hpc/.../md_env/bin/python core/node04_bbs_gap_scan.py

登记的判读程序（先写死，再看数据）：
  1. 只用 pool_class == 'candidate'；在【各电荷分层内】分别排序。
  2. 在每层前 60 名内扫描相邻名次的 BBS 间隔。
  3. 若某处间隔 ≥ 3 × σ_seed(p95)（= 3 × 0.270 = 0.810 kcal/mol），
     判为【真实断层】，可据此设 N 并在 Methods 中如实称之为分布断层。
  4. 若无处满足，则判为【无断层】—— 分布连续。此时 N 必须由算力预算决定，
     且 Methods 必须写作"在可用算力预算内选取各层前 X%"，
     【不得】把预算切割说成分布切割。

为何用 p95 而非中位：门槛处的一次误判代价是整分子的 MD，故取保守噪声上限。
"""
import csv, json, statistics
from pathlib import Path

R = Path("/mnt/hpc/home/wanghongmei/lang2")
SIGMA_P95 = 0.270          # 试点实测（生产 5 种子口径）
GAP_MULT  = 3.0
GAP_THR   = SIGMA_P95 * GAP_MULT
SCAN_TOP  = 60
TG = ["6W63", "3V3M", "4YLU"]

rows = [r for r in csv.DictReader(
    open(R/"results/04_docking/node04_bbs_table.csv", encoding="utf-8"))
    if r.get("pool_class") == "candidate"]
med = {t: statistics.median(float(r[f"mean_{t}"]) for r in rows) for t in TG}
for r in rows:
    r["_raw"] = float(r["BBS"])
    r["_cen"] = max(float(r[f"mean_{t}"]) - med[t] for t in TG)

print(f"候选 {len(rows)} | 断层判据 ≥ {GAP_MULT}×σ_p95 = {GAP_THR:.3f} kcal/mol\n")
report = {"gap_threshold": GAP_THR, "sigma_p95": SIGMA_P95, "strata": {}}

for defn, key in (("raw", "_raw"), ("centred", "_cen")):
    print(f"── 定义 {defn} ──")
    for cs in ("neutral", "charged"):
        sub = sorted((r for r in rows if r["charge_stratum"] == cs), key=lambda r: r[key])
        vals = [r[key] for r in sub[:SCAN_TOP]]
        gaps = sorted(((vals[i+1]-vals[i], i+1) for i in range(len(vals)-1)),
                      reverse=True)[:3]
        hit = [(g, n) for g, n in gaps if g >= GAP_THR]
        print(f"  {cs:<8} n={len(sub):<5} 最优 {vals[0]:+.3f}  "
              f"前3大间隔: " + " ".join(f"{g:.3f}@N={n}" for g, n in gaps))
        print(f"           → {'✅ 断层 N=' + str(hit[0][1]) if hit else '❌ 无断层（分布连续）'}")
        report["strata"][f"{defn}/{cs}"] = {
            "n": len(sub), "best": round(vals[0], 4),
            "top_gaps": [{"gap": round(g, 4), "at_N": n} for g, n in gaps],
            "verdict": ("gap" if hit else "continuous"),
            "gap_at_N": (hit[0][1] if hit else None)}
    print()

any_gap = any(v["verdict"] == "gap" for v in report["strata"].values())
print("=" * 62)
if any_gap:
    print("存在断层 —— 可据分布设 N，Methods 如实称之为分布断层。")
else:
    print("四个（定义 × 分层）组合均【无断层】—— 分布连续。")
    print("→ N 必须由算力预算决定，且 Methods 必须写明这是预算切割，不是分布切割。")
report["overall"] = "gap_found" if any_gap else "continuous_budget_cut_required"

p = R/"results/04_docking/node04_bbs_gap_scan.json"
p.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"\n✓ {p}")
