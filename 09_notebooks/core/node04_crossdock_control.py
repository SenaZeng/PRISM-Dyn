# -*- coding: utf-8 -*-
"""Node 04 · 交叉对接对照（3×3 矩阵）—— 判定"逐受体系统偏移" vs "真实生物学"
    cd /mnt/hpc/home/wanghongmei/lang2
    srun -p cpuPartition -n1 -c4 /mnt/hpc/.../md_env/bin/python core/node04_crossdock_control.py

问题：BBS 有 63.4% 的候选由 3V3M 限制。这是因为 SARS-CoV-1 口袋确实更难结合，
还是因为 3V3M 受体在【所有】配体上都系统性打低一点？自对接对照无法区分 ——
那里每个受体只见过【自己的】配体，三个配体化学上并不相同。

判定实验：把【同一个】配体对接进【三个】受体。同一配体、不同受体 →
分数差异只能归于受体。若某受体在三个配体上一致地打低，即为系统偏移。

复用 node04_redock_control.py 已生成的配体 PDBQT 与参考分子；
对角线（自对接）已算过，幂等跳过，故实际只新增 30 次运行。
"""
import json, sys
from pathlib import Path
R = Path("/mnt/hpc/home/wanghongmei/lang2"); sys.path.insert(0, str(R))
from core import node04_docking as N4                      # noqa: E402

OUT  = R/"results"/"04_docking"/"redock_control"
VINA = Path("/mnt/hpc/home/wanghongmei/apps/vina/vina_1.2.7_linux_x86_64")
LIG_HOME = {"X77":"6W63", "0EN":"3V3M", "R30":"4YLU"}

man = json.loads((R/"results/04_docking/receptors/Node04a_receptor_manifest.json")
                 .read_text(encoding="utf-8"))
grid, rows = {}, []
for het, home in LIG_HOME.items():
    lig = OUT/home/f"{het}.pdbqt"
    if not lig.is_file():
        raise SystemExit(f"缺配体 PDBQT：{lig}（先跑 node04_redock_control.py）")
    for pid, rc in sorted(man["receptors"].items()):
        aff = []
        for seed in N4.SEEDS:
            stem = f"{het}__{pid}__s{seed}"
            row = {"ligand": het, "target": pid, "seed": seed, "row_index": seed,
                   "ligand_pdbqt": str(lig), "ligand_pdbqt_sha256": N4.file_sha256(lig),
                   "receptor_pdbqt": rc["receptor_pdbqt"],
                   "receptor_pdbqt_sha256": rc["conversion"]["pdbqt_sha256"],
                   "box_center": rc["box_center"], "box_size": rc["box_size"],
                   "out_pdbqt": str(OUT/"cross"/pid/f"{stem}.pdbqt"),
                   "out_json":  str(OUT/"cross"/pid/f"{stem}.json")}
            # 自对接结果已在 <home>/ 下，路径不同 → 这里统一重算到 cross/（幂等）
            rec = N4.run_one(row, VINA)
            if rec.get("status") != "ok":
                raise SystemExit(f"{het}→{pid} seed {seed} 失败: {rec.get('error')}")
            aff.append(rec["affinity_best"])
        grid[(het, pid)] = sum(aff)/len(aff)
        rows.append({"ligand": het, "target": pid, "mean": round(grid[(het,pid)],3),
                     "self": het in LIG_HOME and LIG_HOME[het] == pid})

tg = sorted(man["receptors"])
print(f"\n{'配体':<6}" + "".join(f"{t:>10}" for t in tg) + "     (对角线=自对接)")
for het in LIG_HOME:
    print(f"{het:<6}" + "".join(
        f"{grid[(het,t)]:>9.2f}" + ("*" if LIG_HOME[het]==t else " ") for t in tg))

print("\n逐受体均值（跨三个配体 —— 同一批配体，故差异归于受体）:")
base = {t: sum(grid[(h,t)] for h in LIG_HOME)/len(LIG_HOME) for t in tg}
for t in tg: print(f"  {t}: {base[t]:>7.3f}")
spread = max(base.values()) - min(base.values())
worst  = max(base, key=base.get)
print(f"\n极差 {spread:.3f} kcal/mol；最不利受体 = {worst}")
if spread > 0.3:
    print(f"  ⚠ 存在逐受体系统偏移（{spread:.2f} kcal/mol）。若 {worst} 恰是 BBS 的主要"
          "限制性靶点，则该'广谱瓶颈'至少部分是打分偏移而非生物学。")
    print("     建议：在 Methods 中报告该矩阵；考虑以逐受体基线校正后再设 MD 门槛。")
else:
    print("  ✅ 三受体基线相当（<0.3 kcal/mol）—— 限制性靶点分布反映真实结合难度。")

p = OUT/"crossdock_control_manifest.json"
p.write_text(json.dumps({
    "module":"node04_crossdock_control","created_at_utc":N4.utc_now(),
    "purpose":"distinguish per-receptor scoring offset from genuine binding difficulty",
    "vina_sha256": N4.file_sha256(VINA), "seeds": list(N4.SEEDS),
    "exhaustiveness": N4.EXHAUSTIVENESS, "num_modes": N4.NUM_MODES,
    "matrix": [{"ligand":k[0],"target":k[1],"mean_affinity":round(v,4),
                "self_dock": LIG_HOME[k[0]]==k[1]} for k,v in grid.items()],
    "per_receptor_baseline": {t: round(v,4) for t,v in base.items()},
    "baseline_spread_kcal_per_mol": round(spread,4),
    "least_favourable_receptor": worst,
}, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"\n✓ {p}")
