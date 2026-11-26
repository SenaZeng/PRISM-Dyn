# -*- coding: utf-8 -*-
"""Node 04 · BBS 基线校正对照 —— 偏移是否改变【最终选谁】
    /mnt/hpc/.../md_env/bin/python core/node04_bbs_centred.py

背景：置换检验显示"3V3M 限制 63.4% 候选"几乎全部由约 0.2 kcal/mol 的逐受体
偏移解释。于是 BBS 在多数分子上等于 mean_3V3M —— 广谱分数退化为单靶点分数。

但"退化"是否【实际影响结论】，取决于它是否改变最终进入 MD 的那批分子。
本脚本不替换 BBS，只做对照：
    BBS_raw      = max(各靶点均值)                        （现行定义）
    BBS_centred  = max(各靶点均值 − 该靶点全池中位)        （去除系统偏移）
然后比较两者的 top-N 集合重合度。

  重合度高 → 偏移是表述问题，不是选择问题；Methods 报告即可，排名沿用 raw。
  重合度低 → 偏移【改变了】选谁，必须与张老师决定用哪个定义，并写明理由。

⚠ 本脚本【不】主张 centred 更正确。对接无法证明偏移是假象；centred 只是
   一个对照量，用于量化该偏移对下游决定的影响。
"""
import csv, json, statistics
from collections import Counter
from pathlib import Path

R = Path("/mnt/hpc/home/wanghongmei/lang2")
TG = ["6W63", "3V3M", "4YLU"]
TOPN = (28, 50, 100, 200)

rows = [r for r in csv.DictReader(
    open(R / "results/04_docking/node04_bbs_table.csv", encoding="utf-8"))
    if r.get("pool_class") == "candidate"]
print(f"候选 {len(rows)} 个")

vals = {r["ligand"]: {t: float(r[f"mean_{t}"]) for t in TG} for r in rows}
med = {t: statistics.median(v[t] for v in vals.values()) for t in TG}
print("\n逐靶点全池中位（校正基线）:")
for t in TG: print(f"  {t}: {med[t]:+.3f}")

raw = {L: max(v.values()) for L, v in vals.items()}
cen = {L: max(v[t] - med[t] for t in TG) for L, v in vals.items()}
lim_raw = Counter(max(v, key=v.get) for v in vals.values())
lim_cen = Counter(max(TG, key=lambda t: vals[L][t] - med[t]) for L in vals)

print("\n限制性靶点分布:")
print(f"{'靶点':<8}{'raw':>12}{'centred':>12}")
for t in TG:
    print(f"{t:<8}{lim_raw[t]/len(vals)*100:>11.1f}%{lim_cen[t]/len(vals)*100:>11.1f}%")

print("\ntop-N 集合重合（按各自定义升序取最优 N 个）:")
o_raw = sorted(raw, key=raw.get)
o_cen = sorted(cen, key=cen.get)
print(f"{'N':>5}{'共有':>7}{'仅raw':>8}{'仅centred':>11}{'Jaccard':>10}")
res = {}
for n in TOPN:
    a, b = set(o_raw[:n]), set(o_cen[:n])
    j = len(a & b) / len(a | b)
    res[n] = {"shared": len(a & b), "only_raw": len(a - b),
              "only_centred": len(b - a), "jaccard": round(j, 3)}
    print(f"{n:>5}{len(a&b):>7}{len(a-b):>8}{len(b-a):>11}{j:>10.3f}")

print("\n两种排名的 Spearman 秩相关:")
rk_r = {L: i for i, L in enumerate(o_raw)}
rk_c = {L: i for i, L in enumerate(o_cen)}
n = len(vals)
d2 = sum((rk_r[L] - rk_c[L]) ** 2 for L in vals)
rho = 1 - 6 * d2 / (n * (n * n - 1))
print(f"  ρ = {rho:.4f}")

worst = min(res[n]["jaccard"] for n in TOPN)
print("\n判读:")
if worst >= 0.85:
    print("  ✅ top-N 高度重合 —— 偏移是【表述】问题，不改变选谁。")
    print("     排名沿用 BBS_raw；在 Methods 报告偏移与本对照即可。")
elif worst >= 0.6:
    print("  ⚠ top-N 部分变动 —— 偏移【影响边缘分子】。建议两种排名并列报告，")
    print("     并把差集分子单列（它们的入选与否取决于定义选择）。")
else:
    print("  ❌ top-N 显著不同 —— 偏移【决定了】选谁。必须先定义、再看结果，")
    print("     且该决定须与张老师明确记录在案。")

out = R / "results/04_docking/node04_bbs_centred_comparison.json"
out.write_text(json.dumps({
    "purpose": "quantify whether the per-receptor offset changes the selected set",
    "baseline_medians": {t: round(med[t], 4) for t in TG},
    "limiting_target_raw": dict(lim_raw), "limiting_target_centred": dict(lim_cen),
    "topN_overlap": res, "spearman_rho": round(rho, 4),
    "note": ("centred is a CONTROL, not a proposed replacement; docking cannot prove "
             "the offset is artefactual"),
}, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"\n✓ {out}")
