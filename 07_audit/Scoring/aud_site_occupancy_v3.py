#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""aud_site_occupancy_v3.py — 配体是否仍在催化位点？（跨臂 30 条 / 主臂 99 条）

**本脚本只计算一个新量**：配体到催化二元体的最小镜像距离，逐帧。
phase4 的四个轴一个也不重算 —— 需要对照时**引用** phase4 自己写的
`analysis.json`，并在列名上注明来源（`*_from_phase4`）。

---------------------------------------------------------------------------
**v2 修订（2026-10-19）：删除本脚本自算的 `rmsd_lig` 列。**

v1 曾自算 `rmsd_lig` 并在注释中声称「与 phase4 的 groupselections 口径一致」。
**那句话是错的，且我未曾核验。** phase4 用 `rms.RMSD(..., ref_frame=0)`，参考系是
**轨迹第 0 帧**（生产起点）；v1 取的参考是**第 5000 帧**。两者量的是不同的东西
（相对初始位姿的位移 vs 后 50 ns 之内的游走），实测偏离 0.21×–1.39×、毫无规律。

依 D16：磁盘上多一份对不上的副本，正是日后被误当权威的那类东西。故**删除**，
而非「修正后保留」—— phase4 已经拥有那一轴，本脚本不该有第二份。

---------------------------------------------------------------------------
**v3 修订（2026-10-20）：三项。**

**① 删除全部 `PREREG` 引用 —— 该字典从未存在。**
v1/v2 在 `SITE_CUT` 注释与收尾打印中引用了一个 `PREREG['ligand_escape_rule']`
（「3 复本中 ≥2 逃逸 → 该候选降级」）。**实测：该键名在溯源记录中出现 0 次。**
v2 删掉了自算的 `rmsd_lig` 列，却留下了那段打印，于是每次运行都把一条**伪造的
预登记引用**写进日志。本版全部删除，不作「修正后保留」—— 一条不存在的规则
没有正确的引法。逃逸普查改为**纯描述性报告**，不援引任何判据。

**② 受体可配置，且催化二元体不得猜。**
`--arm` / `--p4-root` / `--out` 改为命令行参数，使同一脚本可指向主臂
（`md_solv_v3` / `results/step6e_solv_v3`）而无须改代码。
**6W63 刻意不写入 `DYAD`** —— 主臂的二元体编号须对着**产出那 99 条轨迹的 MD
拓扑**实测，不得从文献或 Node 04 的对接受体沿用（§J2 的 3V3M GLN306、
3V3M 一栏的 `Dropped residues: A:306` 都是先例）。未经 `--dyad` 显式提供即硬失败。

**③ 阈值来源逐次打印。** `site_occupancy` 的 5 Å 判据是**事后**选定的；
本版在每次运行的首行与 CSV 旁的 `.meta.json` 中都记下这一点，使该声明
随产物走，而不是只活在某一份文档里。

---------------------------------------------------------------------------
**v3.1 修订（2026-10-29）：施加已冻结预登记的两项要求。**

依 `claude/PRISM_Dyn_V3_B81_prereg.md`（冻结 2026-10-29 09:42 UTC）。

**① §6 第 2 项 —— 指向冻结文件的断言（不只是注释）。**
主臂（精确等于 `md_solv_v3`）运行时，`--dyad 6W63` 必须提供且必须等于
该登记 §D.1 所记的实测值 **(41, 145)**；不符即显式失败。另逐体系断言
resid 41 属组氨酸族、resid 145 属半胱氨酸族 —— 原版只核验「命中两个残基」，
若拓扑被重编号则编号仍命中两个残基而种类已错，原断言抓不住。
**精确相等，不得用 startswith** —— `md_solv_v3_ref` 以 `md_solv_v3` 为前缀（D28）。

**② §B 第 7 项 —— 到 His163 的逐帧最小距离，与 `d_dyad` 并列。**
**仅作描述，不附阈值，不得改变门禁判定。** 以 `--extra-resid 6W63=163` 启用，
**默认关闭**：跨臂复跑输出因此完全不变，且不替任何受体猜编号
（His163 是 6W63 拓扑中的实测编号；3V3M / 4YLU 的对应残基尚未实测）。

---------------------------------------------------------------------------
**v3.2 修订（2026-10-29）：两项报告层改动，不触及任何判据。**

**① `--only` 子集过滤器。** 已冻结预登记 §3 第 7 项的报告义务是「九条轨迹
逐条」（`057` / `279` / `074`）。按全臂枚举重跑要重读 99 条，为一个描述性列
付 11 倍读盘代价。`--only` 按体系名子串过滤，可重复，且把过滤器写入
`.meta.json` 的 `only_filter`，使子集产物不致日后被误当作全臂产物。

**② `n_rep` 分组改为按（分子, 受体, 姿势）。** v3.1 按（分子, 受体）分组，
而 `pool_mol_011` 与 `pool_mol_020` 各有 primary 与 secondary 两个姿势，
故各报 `n_rep = 6` 并被列为「门禁不可评估」——**伪警报**：各为两套完整的
三复本体系。§1.4 登记的对象是「某分子可用复本少于 3 条」，而承载三复本的
单位是**体系**。此改动使代码与已登记条款更相符，不是放宽它，且因两者皆不在
门禁范围内而**不改变任何判读**。与 D16 同类：检查的作用域须与其对象同源。

---------------------------------------------------------------------------
该项的由来：2026-10-29 提出的定义层面问题 —— 「占据位点」指留在功能性二元体
附近，还是留在配体最初所处的口袋内（其锚定可能在 His163 一侧）？判据保持以
**二元体**为参照（实测：三个分子在第 0 帧同时与 His41 2.81–3.07 Å、
His163 1.97–2.65 Å、Cys145 2.38–3.24 Å 接触，皆远在 5 Å 内，故三者皆非只结 S1 者，
两种定义在起点处无从分歧）；而本列使该问题**可由数据回答**，而非靠事先假定。

---------------------------------------------------------------------------

**要回答的问题**：phase4 的 `rmsd_lig_mean` 在若干条上很大（跨臂最大 34.74 Å），
但那个量**不能**被读作「配体离开了位点」。实测 `pool_mol_004__3V3M/r1`：
phase4 报 10.28 Å，而配体全程待在距二元体 2.4–2.6 Å 处、位点占有 0.998。
原因是该 RMSD **不对配体自身作叠合**，于是把刚体位移与内部构象变化混在一起 ——
一条柔性尾巴摆动、而锚定基团不动，就会产出一个看似逃逸的大数。

**本脚本给出直接的答案**：配体有多少比例的帧仍落在催化二元体 5 Å 之内。
用 `distance_array(..., box=ts.dimensions)` 做最小镜像，**构造上对 PBC 免疫**，
不依赖任何参考帧。

**同时是 §K10 / G24 所缺的那一维**：phase4 的 `hbond_occupancy` 只问「有没有氢键」，
不问「结在哪」。`site_occupancy` 是位点分辨的度量。

催化二元体（依 §B12.3 实测，不猜）：
  3V3M (SARS-CoV-1) : His41 + Cys145
  4YLU (MERS-CoV)   : His41 + Cys148
  6W63 (SARS-CoV-2) : **须实测后以 --dyad 6W63=41,145 形式提供**

输入: <arm>/*/*/06_md/{md.tpr, md_fit.xtc}
      <p4_root>/<sys>/<rep>/analysis.json   （仅引用，不重算）
输出: <out>.csv                             逐轨迹汇总
      <out>.meta.json                       运行参数与阈值来源声明
      <out>_series/<sys>__<rep>.npz         逐帧 d_dyad

用法:
  cd /mnt/hpc/home/wanghongmei/lang2

  # 跨臂（已完成，复跑用）
  python3 Scoring/aud_site_occupancy_v3.py \
      --arm cross_solv_v3 \
      --p4-root results/step6e_cross_v3 \
      --out results/aud_cross_site_occupancy

  # 主臂（须先实测 6W63 二元体，再填入 --dyad）
  python3 Scoring/aud_site_occupancy_v3.py \
      --arm md_solv_v3 \
      --p4-root results/step6e_solv_v3 \
      --out results/aud_main_site_occupancy \
      --dyad 6W63=41,145

  # 先探一条
  python3 Scoring/aud_site_occupancy_v3.py --one <06_md目录> --dyad 6W63=41,145
"""
import argparse
import glob
import json
import pathlib
import sys
import warnings

warnings.filterwarnings("ignore")
import numpy as np

ROOT = pathlib.Path("/mnt/hpc/home/wanghongmei/lang2")

LIG_CANDIDATES = ("MOL", "UNL")   # V3 = MOL / V2 = UNL；按 index.ndx 探测，不写死
STRIDE = 10                       # 与 phase4 的取样间隔一致

# 已实测的受体二元体（§B12.3）。**6W63 刻意缺席** —— 见 v3 修订 ②。
DYAD = {"3V3M": (41, 145), "4YLU": (41, 148)}

SITE_CUT = 5.0   # Å。用于统计「仍在位点」的帧占比。
# ⚠ **该阈值系事后选定** —— 记录中没有任何预登记条目在数值上定义过「逃逸」。
#    报告任何由它导出的数字时，须同时声明该阈值是事后的（依 D13）。
THRESHOLD_PROVENANCE = (
    "site_occupancy 判据 = 配体任一原子到催化二元体任一原子的最小镜像距离 "
    f"<= {SITE_CUT} Å 的帧占比。**该 {SITE_CUT} Å 阈值为事后选定**，"
    "记录中无任何预登记条目在数值上定义「逃逸」（依 D13 须随结果声明）。"
)

# ── 已冻结的预登记（依 §6 第 2 项）────────────────────────────────────────
# 该冻结文件是 Project 文档、不在集群上，故此处不作哈希核验，而是把它的
# **可执行约束**写成断言。标识随产物走（打印 + .meta.json）。
REGISTRATION = ("claude/PRISM_Dyn_V3_B81_prereg.md · 冻结 2026-10-29 09:42 UTC · "
                "§A 可执行条款 · §B §3 第 7 项（His163 列，仅作描述）")

# **精确相等，不得用 startswith** —— md_solv_v3_ref 以 md_solv_v3 为前缀（D28）。
MAIN_ARM = "md_solv_v3"
MAIN_ARM_REC = "6W63"      # 主臂受体；REGISTERED_DYAD 的键，不另写字面量

# §D.1 实测（2026-10-29，三个分子一致；NE2–SG 3.694 / 4.112 / 3.694 Å，
# 次近对 His163–Cys145 相隔 0.5–1.4 Å；His 集与 Node05 冻结记录逐一相符）。
REGISTERED_DYAD = {"6W63": (41, 145)}

# 残基族 —— 用于「编号对、种类错」的断言（重编号探测）。
HIS_RESNAMES = ("HIS", "HID", "HIE", "HIP", "HISD", "HISE", "HISH")
CYS_RESNAMES = ("CYS", "CYM", "CYX", "CYS2")


def parse_dyad_args(items):
    """--dyad 6W63=41,145 [--dyad XXXX=h,c ...] —— 合并进 DYAD，不覆盖已实测项除非显式同值。"""
    out = dict(DYAD)
    for it in items or []:
        if "=" not in it:
            sys.exit(f"--dyad 格式应为 REC=his,cys，收到: {it!r}")
        rec, nums = it.split("=", 1)
        rec = rec.strip()
        try:
            h, c = (int(x) for x in nums.split(","))
        except ValueError:
            sys.exit(f"--dyad {it!r}: 残基编号须为两个整数，形如 6W63=41,145")
        if rec in out and out[rec] != (h, c):
            sys.exit(f"--dyad {rec}={h},{c} 与已实测值 {out[rec]} 冲突 —— "
                     "若确要改，请先更新 §B12.3 并说明依据，不要在命令行上悄悄覆盖。")
        out[rec] = (h, c)
    return out


def parse_extra_args(items):
    """--extra-resid REC=163 —— **仅作描述**的附加残基，按受体。

    默认为空 dict，故不启用时本脚本的行为与 v3 逐字相同。
    不设任何默认编号 —— 不替任何受体猜。
    """
    out = {}
    for it in items or []:
        if "=" not in it:
            sys.exit(f"--extra-resid 格式应为 REC=resid，收到: {it!r}")
        rec, num = it.split("=", 1)
        try:
            out[rec.strip()] = int(num)
        except ValueError:
            sys.exit(f"--extra-resid {it!r}: resid 须为整数，形如 6W63=163")
    return out


def lig_resname(ndx_path):
    import re
    groups = [m.group(1) for m in
              (re.match(r"\s*\[\s*(.+?)\s*\]", ln) for ln in open(ndx_path)) if m]
    hit = [c for c in LIG_CANDIDATES if c in groups]
    if len(hit) != 1:
        raise AssertionError(f"{ndx_path}: 配体组不唯一或缺失，命中 {hit}；现有组 {groups}")
    return hit[0]


def analyse(src: pathlib.Path, dyad_map, p4_root: pathlib.Path, out_ser: pathlib.Path,
            extra_map=None):
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances

    sysname = src.parents[1].name                 # pool_mol_XXX__REC__pose 或 pool_mol_XXX
    rep = src.parents[0].name
    parts = sysname.split("__")
    # 主臂目录名可能不带 __REC__ 段（单受体）。受体须可判定，否则不猜。
    if len(parts) >= 2 and parts[1] in dyad_map:
        rec = parts[1]
    elif len(dyad_map) == 1:
        rec = next(iter(dyad_map))
    else:
        raise AssertionError(
            f"{sysname}: 无法从目录名判定受体，且 --dyad 提供了多个受体 {sorted(dyad_map)}。"
            " 催化二元体编号不可猜 —— 请按受体拆分运行，或在目录名中带上受体段。")
    h_resid, c_resid = dyad_map[rec]

    lig = lig_resname(src / "index.ndx")
    u = mda.Universe(str(src / "md.tpr"), str(src / "md_fit.xtc"))
    n_total = len(u.trajectory)
    start = n_total // 2                          # 与 phase4 同：10001 -> 5000

    ag_lig = u.select_atoms(f"resname {lig}")
    assert len(ag_lig) > 0, f"{sysname}/{rep}: 未选到配体 {lig}"
    ag_dyad = u.select_atoms(f"protein and resid {h_resid} {c_resid}")
    assert len(ag_dyad) > 0, f"{sysname}/{rep}: 未选到 His{h_resid}/Cys{c_resid}"
    # 两个残基都须命中，否则「二元体距离」实际只量了一个残基且不会报错
    n_res = len(set(ag_dyad.resids))
    assert n_res == 2, (f"{sysname}/{rep}: resid {h_resid} {c_resid} 只命中 {n_res} 个残基 —— "
                        "该拓扑的编号与预期不符，须实测后再跑")
    # 残基**种类**断言（v3.1，依预登记 §6 第 2 项 / §4）。
    # 「命中两个残基」不足以保证编号未被平移：重编号后编号仍会命中两个残基，
    # 而种类已错，上一行的断言抓不住。
    _rn = {int(r.resid): str(r.resname).upper() for r in ag_dyad.residues}
    assert _rn.get(h_resid) in HIS_RESNAMES, (
        f"{sysname}/{rep}: resid {h_resid} 的残基名为 {_rn.get(h_resid)!r}，不属组氨酸族 "
        f"{HIS_RESNAMES} —— 该拓扑的编号与 {REGISTRATION} 所记不符，停止")
    assert _rn.get(c_resid) in CYS_RESNAMES, (
        f"{sysname}/{rep}: resid {c_resid} 的残基名为 {_rn.get(c_resid)!r}，不属半胱氨酸族 "
        f"{CYS_RESNAMES} —— 该拓扑的编号与 {REGISTRATION} 所记不符，停止")

    # ── §B 第 7 项：附加残基的逐帧最小距离（**仅作描述**）────────────────
    # 默认不启用。启用时亦不附任何阈值，不参与任何判定。
    ag_extra, extra_resid, extra_resname = None, None, None
    extra_resid = (extra_map or {}).get(rec)
    if extra_resid is not None:
        ag_extra = u.select_atoms(f"protein and resid {extra_resid}")
        assert len(ag_extra) > 0, f"{sysname}/{rep}: 未选到 resid {extra_resid}"
        assert len(set(ag_extra.resids)) == 1, (
            f"{sysname}/{rep}: resid {extra_resid} 命中 {len(set(ag_extra.resids))} 个残基")
        extra_resname = str(ag_extra.residues[0].resname).upper()

    d_dyad, _d_extra = [], []
    for ts in u.trajectory[start::STRIDE]:
        dm = distances.distance_array(ag_lig.positions, ag_dyad.positions, box=ts.dimensions)
        d_dyad.append(float(dm.min()))
        if ag_extra is not None:
            de = distances.distance_array(ag_lig.positions, ag_extra.positions,
                                          box=ts.dimensions)
            _d_extra.append(float(de.min()))
    d_dyad = np.array(d_dyad)
    d_extra = np.array(_d_extra) if _d_extra else None

    out_ser.mkdir(parents=True, exist_ok=True)
    _npz = dict(d_dyad=d_dyad, start=start, stride=STRIDE,
                dyad_resids=np.array([h_resid, c_resid]))
    if d_extra is not None:
        _npz["d_extra"] = d_extra
        _npz["extra_resid"] = np.array([extra_resid])
    np.savez_compressed(out_ser / f"{sysname}__{rep}.npz", **_npz)

    # phase4 的值：**引用**，不重算。取不到就留空，绝不自己算一个顶上。
    p4 = p4_root / sysname / rep / "analysis.json"
    q = {}
    if p4.is_file():
        try:
            j = json.load(open(p4))
            q = {f"{k}_from_phase4": j.get(k) for k in
                 ("rmsd_lig_mean", "rmsd_prot_mean", "mindist_mean", "hbond_occupancy")}
        except Exception:
            q = {}

    # ── §B 第 7 项的描述性列 ──────────────────────────────────────────────
    # 两个计数列**只复用已登记的 5.0 Å**，不引入任何新切点；以计数而非比例
    # 呈现，正是为了不让它被误读为第二道判据。它们回答的是：
    # 「离开二元体的那些帧里，有多少仍在该附加残基 5 Å 之内」。
    _extra_cols = {}
    if d_extra is not None:
        _off = d_dyad > SITE_CUT
        _extra_cols = dict(
            extra_resid=extra_resid,
            extra_resname=extra_resname,
            d_extra_min=round(float(d_extra.min()), 3),
            d_extra_mean=round(float(d_extra.mean()), 3),
            d_extra_max=round(float(d_extra.max()), 3),
            n_frames_offdyad=int(_off.sum()),
            n_frames_offdyad_near_extra=int((_off & (d_extra <= SITE_CUT)).sum()),
        )

    # v3.2：姿势段。目录名形如 pool_mol_057__6W63__primary；缺失时记为 primary。
    # 承载「三条复本」的单位是**体系**（分子 × 受体 × 姿势），§1.4 的分组须用它。
    pose = parts[2] if len(parts) > 2 else "primary"

    return dict(
        system=sysname, rep=rep, target=rec, mol_id=parts[0], pose=pose,
        ligand_resname=lig,
        dyad=f"His{h_resid}+Cys{c_resid}",
        n_frames=len(d_dyad),
        d_dyad_mean=round(float(d_dyad.mean()), 3),
        d_dyad_first=round(float(d_dyad[0]), 3),
        d_dyad_last=round(float(d_dyad[-10:].mean()), 3),
        d_dyad_max=round(float(d_dyad.max()), 3),
        site_occupancy=round(float((d_dyad <= SITE_CUT).mean()), 4),
        **_extra_cols,
        **q)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="cross_solv_v3",
                    help="轨迹臂目录名，相对 ROOT（如 cross_solv_v3 / md_solv_v3）")
    ap.add_argument("--p4-root", default="results/step6e_cross_v3",
                    help="phase4 产物根，相对 ROOT（只读引用）")
    ap.add_argument("--out", default="results/aud_cross_site_occupancy",
                    help="输出前缀，相对 ROOT；将写 <out>.csv / <out>.meta.json / <out>_series/")
    ap.add_argument("--dyad", action="append", default=[],
                    help="补充受体二元体，形如 6W63=41,145。主臂必须提供（须实测，不得沿用文献）")
    ap.add_argument("--extra-resid", action="append", default=[],
                    help="【仅作描述】额外报告配体到该残基的逐帧最小距离，形如 "
                         "6W63=163。不附阈值、不改变门禁判定（预登记 §B 第 7 项）。"
                         "默认不启用；不替任何受体猜编号。")
    ap.add_argument("--only", action="append", default=[],
                    help="只处理体系名含该子串者，可重复（如 --only pool_mol_057 "
                         "--only pool_mol_279）。用于履行只针对门禁范围的报告义务，"
                         "而无须重跑全臂。过滤器写入 .meta.json 的 only_filter。"
                         "**用它时请换一个 --out 前缀，不要覆盖全臂产物。**")
    ap.add_argument("--one", default=None, help="只跑一条（06_md 目录绝对路径）")
    ap.add_argument("--expect-reps", type=int, default=3,
                    help="每个分子×受体期望的复本数（默认 3）。数目不符即显式报错，"
                         "不得静默以更少的复本成表 —— 依 D16 与预登记 §1.4")
    a = ap.parse_args()

    dyad_map = parse_dyad_args(a.dyad)
    extra_map = parse_extra_args(a.extra_resid)

    # ── 已冻结预登记的闸（§6 第 2 项）─────────────────────────────────────
    # **精确相等**，不得用 startswith：md_solv_v3_ref 以 md_solv_v3 为前缀（D28）。
    if a.arm == MAIN_ARM:
        for _rec, _pair in REGISTERED_DYAD.items():
            if _rec not in dyad_map:
                sys.exit(f"主臂（{MAIN_ARM}）须显式提供 "
                         f"--dyad {_rec}={_pair[0]},{_pair[1]} —— 依 {REGISTRATION}")
            if dyad_map[_rec] != _pair:
                sys.exit(f"--dyad {_rec}={dyad_map[_rec][0]},{dyad_map[_rec][1]} 与已冻结登记 "
                         f"§D.1 所记的实测值 {_pair} 不符。该编号于 2026-10-29 实测、"
                         f"三个分子一致。若确要改，须先以带日期的修订更新 "
                         f"{REGISTRATION}，不得在命令行上静默改动。")
        print(f"[登记] {REGISTRATION}")
        print(f"[登记] 主臂催化二元体已核验 = "
              f"His{REGISTERED_DYAD[MAIN_ARM_REC][0]} + Cys{REGISTERED_DYAD[MAIN_ARM_REC][1]}"
              f"（§D.1 实测，三个分子一致）")
    arm = ROOT / a.arm
    p4_root = ROOT / a.p4_root
    out_csv = ROOT / (a.out + ".csv")
    out_meta = ROOT / (a.out + ".meta.json")
    out_ser = ROOT / (a.out + "_series")

    srcs = ([pathlib.Path(a.one)] if a.one else
            [pathlib.Path(p).parent
             for p in sorted(glob.glob(str(arm / "*/*/06_md/md_fit.xtc")))])
    assert srcs, f"未找到 md_fit.xtc（{arm}）—— phase2 是否已完成？"

    # v3.2：子集过滤。全臂枚举之后再筛，故「找不到 md_fit.xtc」与「过滤后为空」
    # 是两条不同的错误信息 —— 不把二者混成一句。
    if a.only:
        _before = len(srcs)
        srcs = [s for s in srcs
                if any(k in s.parents[1].name for k in a.only)]
        assert srcs, f"--only {a.only} 未匹配任何体系（过滤前 {_before} 条）"
        print(f"--only 过滤: {_before} -> {len(srcs)} 条   匹配 {a.only}")
        print("  ⚠ 这是一份**子集**产物。过滤器已写入 .meta.json 的 only_filter；"
              "请确认 --out 前缀与全臂产物不同，以免覆盖。")

    # 受体可判定性预检：在算之前就失败，而不是算到一半
    unknown = set()
    for s in srcs:
        parts = s.parents[1].name.split("__")
        if len(parts) >= 2 and parts[1] not in dyad_map and len(dyad_map) != 1:
            unknown.add(parts[1])
    if unknown:
        sys.exit(f"以下受体无已实测的催化二元体: {sorted(unknown)}\n"
                 f"已知: {dyad_map}\n"
                 "→ 须先对着**产出这些轨迹的 MD 拓扑**实测其 His/Cys 编号，"
                 "再以 --dyad REC=h,c 提供。不得从文献或对接受体沿用。")

    print(f"待分析 {len(srcs)} 条 | 臂 {a.arm} | 催化二元体 {dyad_map}")
    if extra_map:
        print(f"附加残基（**仅作描述，不附阈值，不改变门禁判定**）: {extra_map}")
        print("  新增列 d_extra_min/mean/max 与两个计数列；两个计数列只复用已登记的 "
              f"{SITE_CUT} Å，不引入新切点。逐帧 d_extra 另存于同一 .npz。")
    print(THRESHOLD_PROVENANCE + "\n", flush=True)

    rows = []
    for i, s in enumerate(srcs, 1):
        try:
            r = analyse(s, dyad_map, p4_root, out_ser, extra_map)
        except Exception as e:
            print(f"[{i}/{len(srcs)}] ERR {s.parents[1].name}/{s.parents[0].name}"
                  f"  {type(e).__name__}: {e}", flush=True)
            continue
        rows.append(r)
        v = ("位点内" if r["site_occupancy"] >= 0.9 else
             ("部分离开" if r["site_occupancy"] >= 0.5 else "已离开位点"))
        p4v = r.get("rmsd_lig_mean_from_phase4")
        p4s = f" | phase4 rmsd_lig {p4v:6.2f}" if isinstance(p4v, (int, float)) else ""
        print(f"[{i}/{len(srcs)}] {r['system'].replace('__primary',''):<22s}/{r['rep']:<3s}"
              f" 二元体 首 {r['d_dyad_first']:6.2f} → 末 {r['d_dyad_last']:6.2f} Å"
              f" | 位点占有 {r['site_occupancy']:.3f} ({v}){p4s}", flush=True)

    if not rows:
        sys.exit("无可用结果")

    import pandas as pd
    df = pd.DataFrame(rows).sort_values(["mol_id", "target", "rep"])
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)

    json.dump(dict(arm=a.arm, p4_root=a.p4_root, dyad={k: list(v) for k, v in dyad_map.items()},
                   site_cut_angstrom=SITE_CUT, stride=STRIDE,
                   n_trajectories=len(rows),
                   threshold_provenance=THRESHOLD_PROVENANCE,
                   registration=REGISTRATION,
                   only_filter=list(a.only),
                   is_subset=bool(a.only),
                   extra_resid={k: int(v) for k, v in (extra_map or {}).items()},
                   extra_resid_note=("descriptive only; no threshold attached; "
                                     "does not modify the gate verdict; the two count "
                                     "columns reuse the registered 5.0 A cut and "
                                     "introduce no new cut-point (prereg \u00a7B item 7)")),
              open(out_meta, "w"), ensure_ascii=False, indent=1)

    # ── 逐（分子, 受体）描述性普查 ────────────────────────────────────────
    # **纯描述**。不援引任何通过/不通过判据 —— 记录中不存在以数值定义「逃逸」
    # 的预登记条目，而在见过数据之后现编一个，只是同一种失效换件外衣。
    print("\n" + "=" * 76)
    print(f"逐（分子, 受体, 姿势）—— site_occupancy < 0.5 的复本数（{SITE_CUT} Å 事后阈值下）")
    print("=" * 76)
    # v3.2：按**体系**（分子 × 受体 × 姿势）分组。按（分子, 受体）分组会把
    # primary 与 secondary 合并，使 011 / 020 报 n_rep = 6 并触发伪警报。
    g = (df.assign(below_half=(df.site_occupancy < 0.5).astype(int))
           .groupby(["mol_id", "target", "pose"])
           .agg(n_rep=("rep", "count"), n_below_half=("below_half", "sum"),
                site_occ_min=("site_occupancy", "min"),
                d_dyad_med=("d_dyad_mean", "median")).round(3).reset_index())
    print(g.to_string(index=False))

    # ── 复本数完整性闸（依 D16 与预登记 §1.4）────────────────────────────
    # 逐条 try/except 会在失败时打印 ERR 并继续，故一条轨迹失败会使该分子
    # **静默地**以 n_rep = 2 成表。那正是「检查的作用域与其对象不同源」。
    # 此处显式比对期望值，并把不足者单列。
    short = g[g.n_rep != a.expect_reps]
    if len(short):
        print("\n" + "!" * 76)
        print(f"⚠ 以下 (分子, 受体) 的可用复本数 != {a.expect_reps} —— "
              "**门禁对这些分子不可评估**（预登记 §1.4）：")
        print(short[["mol_id", "target", "pose", "n_rep"]].to_string(index=False))
        print("措辞须取 0/3 档的留置措辞，并写明可用复本数；"
              "**不得以可用复本的比例代替 3/3**。")
        print("!" * 76)

    n_tot = int(g.n_rep.sum())
    n_below = int(g.n_below_half.sum())
    print(f"\n合计 {n_tot - n_below}/{n_tot} 条复本 site_occupancy >= 0.5"
          f"（{n_below} 条低于 0.5）。**以上为描述性报告，不构成任何判据的通过或不通过。**")

    print("\n⚠ phase4 的 `rmsd_lig` 不可读作「配体离开位点」——"
          "它不对配体自身叠合，混入了内部构象变化。")
    if "rmsd_lig_mean_from_phase4" in df:
        odd = df[(df["rmsd_lig_mean_from_phase4"].fillna(0) > 5) & (df.site_occupancy >= 0.9)]
        for _, r in odd.iterrows():
            print(f"   · {r.system.replace('__primary','')}/{r.rep}: "
                  f"phase4 rmsd_lig {r.rmsd_lig_mean_from_phase4:.2f} Å，"
                  f"但位点占有 {r.site_occupancy:.3f}、二元体距离 {r.d_dyad_mean:.2f} Å —— 从未离开")

    print(f"\n已写 {out_csv}")
    print(f"     {out_meta}")
    print(f"逐帧 d_dyad 序列 {out_ser}/  （{len(rows)} 个 .npz）")


if __name__ == "__main__":
    main()
