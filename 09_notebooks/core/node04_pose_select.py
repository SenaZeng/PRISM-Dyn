# -*- coding: utf-8 -*-
"""Node 04 → 05 · MD 起始姿势选择
    预览（只出簇结构与计数，不出坐标）:
        python core/node04_pose_select.py --target 6W63
    定案（写出选择 + manifest）:
        python core/node04_pose_select.py --target 6W63 --commit

═══ 预先登记的规则（在查看任何簇结构之前写死）═══
1. 姿势来源 = exhaustiveness 128 对照臂。
   理由：排名需要【分数】收敛（exh32 已由 σ_seed=0.057 与自对接对照确立）；
         MD 起始结构需要【姿势】收敛，而 exh128 可度量地更好
         （发散 25%→19%，质心 p90 4.09→0.68 Å）。任务不同，要求不同。
2. 对每个复合物，取 5 个种子【各自的第 1 构象】（即搜索实际"选中"的姿势），
   以 2.0 Å 原子序 RMSD 作单链接聚类。
3. 主姿势 = 成员最多的簇中【分数最优】者。
   —— 反映"搜索最常找到哪个模式"，而非"哪一次运行碰巧得分最高"。
4. 次姿势 = 第二大簇中分数最优者，仅当该簇含 ≥2 个种子（即 ≥40%）时才产生。
   4/1 分裂视为单次离群，不产生次姿势；3/2 分裂则产生。
5. 有次姿势的复合物，额外跑 3 个 replica 从次姿势出发（方案 C）。
   理由：翻转类在 4 倍搜索强度下计数一个不减（11→11），说明打分函数确实
         分辨不出两个朝向。为其登记一个"选哪个"的规则，等于把掷硬币
         记录为决定。MM-PBSA 具备区分朝向所需的静电与去溶剂化，让它裁决。

═══ v1.2 并列判据（澄清，非规则变更）═══
簇规模相等时（本数据集中仅 pool_mol_007 / 4YLU 的 2/2/1 一例），
主簇取含分数最优姿势者。核验：该例主姿势原为 s42 −9.56、次为 s46 −8.99，
新判据下结果相同。三个靶点 93 个复合物中无任何选择因此改变。

═══ 反循环 ═══
不带 --commit 时只输出簇的大小与分数，不输出坐标或文件路径；
规则已写死于上，故预览不构成"看到结果再定规则"。
"""
import argparse
import json
import sys
from itertools import combinations
from pathlib import Path

R = Path("/mnt/hpc/home/wanghongmei/lang2")
OUT = R/"results"/"04_docking"
ARM = OUT/"exhaustiveness_arm"
SEEDS = (42, 43, 44, 45, 46)
CLUSTER_A = 2.0        # 单链接聚类阈值
MIN_SEEDS_SECONDARY = 2


def top_mode(pdbqt):
    """第 1 个 MODEL 的重原子坐标 + 其亲和力。"""
    xyz, aff, started = [], None, False
    for ln in Path(pdbqt).read_text().splitlines():
        if ln.startswith("MODEL"):
            if started: break
            started = True
        elif ln.startswith("REMARK VINA RESULT:") and aff is None:
            aff = float(ln.split(":", 1)[1].split()[0])
        elif ln.startswith(("ATOM", "HETATM")):
            if (ln.split()[-1] if len(ln) > 70 else "") in ("HD", "H"):
                continue
            xyz.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
        elif ln.startswith("ENDMDL") and started:
            break
    return xyz, aff


def rmsd(a, b):
    n = len(a)
    return (sum((a[i][k]-b[i][k])**2 for i in range(n) for k in range(3))/n) ** 0.5


def single_link(P, thr):
    """单链接聚类，返回 [[idx,...], ...]，按成员数降序。"""
    n = len(P); par = list(range(n))
    def find(i):
        while par[i] != i: par[i] = par[par[i]]; i = par[i]
        return i
    for i, j in combinations(range(n), 2):
        if rmsd(P[i], P[j]) < thr:
            a, b = find(i), find(j)
            if a != b: par[a] = b
    cl = {}
    for i in range(n): cl.setdefault(find(i), []).append(i)
    return sorted(cl.values(), key=lambda v: -len(v))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="6W63")
    ap.add_argument("--commit", action="store_true")
    a = ap.parse_args()
    T = a.target

    sel = json.loads((OUT/"Node04_md_selection.json").read_text(encoding="utf-8"))
    ligs = [d["ligand"] for d in sel["selected"]]
    print(f"靶点 {T} | 分子 {len(ligs)} | 姿势来源 exhaustiveness 128 对照臂\n")

    out, n_sec = [], 0
    for L in sorted(ligs):
        P, A, ok = [], [], True
        for s in SEEDS:
            p = ARM/"poses"/T/f"{L}__{T}__s{s}.pdbqt"
            if not p.is_file(): ok = False; break
            xyz, aff = top_mode(p)
            P.append(xyz); A.append((aff, s))
        if not ok or len({len(x) for x in P}) != 1:
            print(f"  ⚠ {L}: 姿势不全，跳过"); continue
        cl = single_link(P, CLUSTER_A)
        # ★ v1.2：簇大小相等时（如 2/2/1），原先由 union-find 的 root 序决定谁是
        #   "第一簇" —— 那是实现细节，不是规则。补上显式并列判据：
        #   规模相同时，取【含分数最优姿势】的那个簇为主簇。
        #   这与规则原意一致（"搜索最常找到的模式，同分时取分数更优"），
        #   且经核验【不改变任何已产出的选择】——见 docstring 末尾说明。
        cl.sort(key=lambda c: (-len(c), min(A[i][0] for i in c)))
        prim_c = cl[0]
        prim = min(prim_c, key=lambda i: A[i][0])
        # 记录主姿势的【可复现支持度】——供下游解读，不参与规则判定。
        rec = {"ligand": L, "target": T, "n_clusters": len(cl),
               "cluster_sizes": [len(c) for c in cl],
               "primary_seed": A[prim][1], "primary_affinity": A[prim][0],
               "n_seeds_in_primary": len(prim_c),
               "primary_support": round(len(prim_c)/len(SEEDS), 2)}
        if len(cl) > 1 and len(cl[1]) >= MIN_SEEDS_SECONDARY:
            sec = min(cl[1], key=lambda i: A[i][0])
            rec.update({"secondary_seed": A[sec][1], "secondary_affinity": A[sec][0],
                        "inter_cluster_rmsd": round(rmsd(P[prim], P[sec]), 2)})
            n_sec += 1
        out.append(rec)

    print(f"{'分子':<16}{'簇结构':>10}{'支持度':>8}{'主(种子/分数)':>18}"
          f"{'次(种子/分数)':>18}{'簇间RMSD':>10}{'Δ分数':>9}")
    for r in out:
        cs = "/".join(str(x) for x in r["cluster_sizes"])
        p = f"s{r['primary_seed']} {r['primary_affinity']:.2f}"
        s = (f"s{r['secondary_seed']} {r['secondary_affinity']:.2f}"
             if "secondary_seed" in r else "—")
        d = f"{r['inter_cluster_rmsd']:.2f}" if "inter_cluster_rmsd" in r else "—"
        da = (f"{abs(r['primary_affinity']-r['secondary_affinity']):.2f}"
              if "secondary_affinity" in r else "—")
        print(f"{r['ligand']:<16}{cs:>10}{r['primary_support']*100:>7.0f}%"
              f"{p:>18}{s:>18}{d:>10}{da:>9}")

    n = len(out); tr = n*3 + n_sec*3
    weak = [r for r in out if r["primary_support"] < 1.0]
    print(f"\n复合物 {n} | 产生次姿势 {n_sec} 个 | 主姿势支持度 <100% 的 {len(weak)} 个")
    if weak:
        print("  " + ", ".join(f"{r['ligand'].replace('pool_mol_','')}"
                               f"({r['primary_support']*100:.0f}%)" for r in weak))
    print(f"轨迹 = {n}×3 主 + {n_sec}×3 次 = {tr} 条 | {tr*24:,}–{tr*48:,} GPU-h"
          f"（较单姿势 +{n_sec*3/(n*3)*100:.0f}%）")

    if not a.commit:
        print("\n[预览] 未写盘。确认后加 --commit。")
        return 0

    man = {"module": "node04_pose_select", "target": T,
           "created_at_utc": __import__("datetime").datetime.now(
               __import__("datetime").timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "pose_source": "exhaustiveness_arm (exhaustiveness=128)",
           "rule": {"cluster_threshold_angstrom": CLUSTER_A,
                    "clustered_on": "top mode of each of 5 seeds",
                    "primary": "best-scoring member of the largest cluster",
                    "secondary": f"best-scoring member of cluster 2, only if it holds "
                                 f">={MIN_SEEDS_SECONDARY} seeds",
                    "rationale": ("flip-in-place count was invariant under 4x search effort "
                                  "(11 -> 11), so the scoring function cannot discriminate "
                                  "orientation; MM-PBSA adjudicates instead of a coin flip"),
                    "pre_registered": True,
                    "tie_break": "equal-size clusters resolved by best score within cluster (v1.2 clarification; changes no outcome in this dataset)"},
           "n_complexes": n, "n_with_secondary": n_sec,
           "n_primary_support_below_1": len(weak), "n_trajectories": tr,
           "selection": out}
    p = OUT/f"Node04_pose_selection_{T}.json"
    p.write_text(json.dumps(man, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n✓ {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
