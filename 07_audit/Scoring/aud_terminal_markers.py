#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
aud_terminal_markers.py —— 规则 16 第三项判据（完成标记）的**独立重做**
2026-10-30 · 因 `aud_frame_census.py` 的 C3 只数一种终止标记而误报，故另立本工具

为何另立而非修改原工具（两条项目规约）：
  · `aud_frame_census.py` 以 sha256 锚定，不得编辑；须包装或另立。
  · 一道检查不得与其所检之物共享代码。本工具的遍历、标记解析与体系枚举
    皆独立重写，不 import 原工具的任何部分。

**原工具 C3 的全部判定作废**，包括其对跨臂的 30/30 与参考臂的 12/12——
那两个恰好为对，但它们从未对着出错的那条路径被检验过。
`Scoring/frame_census_main_20261030.log` 的 VERDICT 行（"NOT ALL THREE PASS"）
系由该缺陷产生，**不得在稿件中引用**；该日志的 C1 与 C2 仍有效。

本工具所核：每条轨迹是否存在**至少一种**合法终止标记。
合法终止状态有二，皆由生产运行器写出：
  [MD2 DONE]     —— 新跑完成
  [RESUME DONE]  —— 自检查点续跑完成
二者的路径写法不同：前者记至复本目录（…/<sys>/<rep>），
后者记至其 06_md 子目录（…/<sys>/<rep>/06_md）。**两种写法须分别匹配。**

用法
  python3 aud_terminal_markers.py --selftest
  python3 aud_terminal_markers.py --arm md_solv_v3
  python3 aud_terminal_markers.py --arm cross_solv_v3 --expect-systems 12 --expect-traj 36

退出码
  0 通过 · 2 未通过 · 3 未执行核查（自测失败或参数不足）
"""
import argparse
import collections
import glob
import os
import re
import shutil
import sys
import tempfile

ROOT_DEFAULT = "/mnt/hpc/home/wanghongmei/lang2"

# 注册表只认已查实的两臂；其余须显式给出期望数（与原工具同一政策，独立实现）
REGISTRY = {"md_solv_v3": (33, 99), "md_solv_v3_ref": (4, 12)}

SYS = r"([A-Za-z0-9_]+__[0-9A-Za-z]+__[a-z]+)"
# 两种终止标记，分别匹配其各自的路径写法
MARKERS = {
    "MD2":    re.compile(r"\[MD2 DONE\][^\n]*?/" + SYS + r"/(r\d+)(?![^\s]*06_md)"),
    "RESUME": re.compile(r"\[RESUME DONE\][^\n]*?/" + SYS + r"/(r\d+)/06_md"),
}


def scan(arm_dir):
    """返回 {system: {replica: set(marker kinds)}}，自 <arm>/logs/ 下全部普通文件。"""
    seen = collections.defaultdict(lambda: collections.defaultdict(set))
    logdir = os.path.join(arm_dir, "logs")
    if not os.path.isdir(logdir):
        return seen, 0
    nfiles = 0
    for name in sorted(os.listdir(logdir)):
        p = os.path.join(logdir, name)
        if not os.path.isfile(p):
            continue
        nfiles += 1
        try:
            with open(p, errors="replace") as fh:
                text = fh.read()
        except OSError:
            continue
        for kind, rx in MARKERS.items():
            for m in rx.finditer(text):
                seen[m.group(1)][m.group(2)].add(kind)
    return seen, nfiles


def systems_on_disk(arm_dir):
    out = set()
    for p in glob.glob(os.path.join(arm_dir, "*", "*", "06_md", "md.xtc")):
        if "abandon" in p.lower():
            continue
        out.add(os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(p)))))
    return sorted(out)


def replicas_on_disk(arm_dir, system):
    out = set()
    for p in glob.glob(os.path.join(arm_dir, system, "*", "06_md", "md.xtc")):
        if "abandon" in p.lower():
            continue
        out.add(os.path.basename(os.path.dirname(os.path.dirname(p))))
    return sorted(out)


def audit(root, arm, exp_sys, exp_traj, verbose=True):
    arm_dir = os.path.join(root, arm)
    if not os.path.isdir(arm_dir):
        return None, ["arm directory absent: %s" % arm_dir]
    seen, nfiles = scan(arm_dir)
    systems = systems_on_disk(arm_dir)
    problems = []

    if exp_sys is not None and len(systems) != exp_sys:
        problems.append("systems on disk %d, expected %d" % (len(systems), exp_sys))

    total_traj = 0
    total_marked = 0
    bykind = collections.Counter()
    for s in systems:
        reps = replicas_on_disk(arm_dir, s)
        total_traj += len(reps)
        for r in reps:
            kinds = seen.get(s, {}).get(r, set())
            if kinds:
                total_marked += 1
                bykind[",".join(sorted(kinds))] += 1
            else:
                problems.append("no terminal marker: %s/%s" % (s, r))

    if exp_traj is not None and total_traj != exp_traj:
        problems.append("trajectories on disk %d, expected %d" % (total_traj, exp_traj))

    # 标记指向磁盘上不存在的复本 —— 反向错误，同样须报
    for s, reps in seen.items():
        for r in reps:
            if s not in systems or r not in replicas_on_disk(arm_dir, s):
                problems.append("marker for absent trajectory: %s/%s" % (s, r))

    if verbose:
        print("  arm                : %s" % arm)
        print("  log files scanned  : %d" % nfiles)
        print("  systems on disk    : %d" % len(systems))
        print("  trajectories       : %d" % total_traj)
        print("  with a marker      : %d / %d" % (total_marked, total_traj))
        print("  by marker kind     : %s" % (dict(bykind) or "{}"))
    return (total_marked, total_traj, bykind), problems


# --------------------------------------------------------------------------
# 自测：一道闸门须被证明能够失败，且须在**出错的那条路径**上被证明
# --------------------------------------------------------------------------
def _single_marker_logic(text):
    """旧工具的做法：只数 [MD2 DONE]。此处重建之，只为证明它在夹具上会误报。"""
    rx = re.compile(r"\[MD2 DONE\][^\n]*?/" + SYS + r"/(r\d+)")
    return {(m.group(1), m.group(2)) for m in rx.finditer(text)}


def selftest():
    print("=" * 70)
    print("自测 —— 三个夹具，其中夹具 B 正是使原工具误报的那条路径")
    print("=" * 70)
    tmp = tempfile.mkdtemp(prefix="aud_tm_")
    ok = True
    try:
        arm = os.path.join(tmp, "fixture_arm")
        for s, reps in [("mol_A__6W63__primary", ["r1", "r2", "r3"]),
                        ("mol_B__6W63__primary", ["r1", "r2", "r3"])]:
            for r in reps:
                d = os.path.join(arm, s, r, "06_md")
                os.makedirs(d)
                open(os.path.join(d, "md.xtc"), "w").write("x")
        logs = os.path.join(arm, "logs")
        os.makedirs(logs)

        # A: 三条皆 MD2 完成
        with open(os.path.join(logs, "prod_1.out"), "w") as fh:
            for r in ("r1", "r2", "r3"):
                fh.write("[MD2 DONE] /x/fixture_arm/mol_A__6W63__primary/%s branch=staged 2026-01-01T00:00:00\n" % r)
        # B: r1 为 MD2，r2/r3 仅有 RESUME —— 真完整，而旧逻辑会漏掉两条
        with open(os.path.join(logs, "prod_2.out"), "w") as fh:
            fh.write("[MD2 DONE] /x/fixture_arm/mol_B__6W63__primary/r1 branch=staged 2026-01-01T00:00:00\n")
        with open(os.path.join(logs, "resume_2.out"), "w") as fh:
            for r in ("r2", "r3"):
                fh.write("[RESUME] /x/fixture_arm/mol_B__6W63__primary/%s/06_md @ gpu01 2026-01-02T00:00:00\n" % r)
                fh.write("[RESUME DONE] /x/fixture_arm/mol_B__6W63__primary/%s/06_md 2026-01-02T02:00:00\n" % r)

        alltext = "".join(open(os.path.join(logs, f)).read() for f in sorted(os.listdir(logs)))

        # ① 旧逻辑必须在夹具 B 上漏报，否则本自测毫无意义
        old = _single_marker_logic(alltext)
        old_B = {r for (s, r) in old if s == "mol_B__6W63__primary"}
        if old_B == {"r1"}:
            print("  [1] 旧的单标记逻辑在夹具 B 上只认出 r1 —— 误报已复现  PASS")
        else:
            print("  [1] 旧逻辑在夹具 B 上认出 %s，与预期 {'r1'} 不符 —— 自测无意义  FAIL" % sorted(old_B))
            ok = False

        # ② 新逻辑必须认出全部六条
        res, prob = audit(tmp, "fixture_arm", 2, 6, verbose=False)
        if res and res[0] == 6 and res[1] == 6 and not prob:
            print("  [2] 新逻辑认出 6/6，无问题 —— 两种标记并计有效      PASS")
        else:
            print("  [2] 新逻辑得 %s，问题 %s —— FAIL" % (res, prob))
            ok = False
        if res and dict(res[2]) == {"MD2": 4, "RESUME": 2}:
            print("  [3] 标记种类分布 {'MD2': 4, 'RESUME': 2} 正确          PASS")
        else:
            print("  [3] 标记种类分布 %s 不符 —— FAIL" % (dict(res[2]) if res else None))
            ok = False

        # ④ 删去一条 RESUME 标记后必须失败 —— 闸门须被证明能够失败
        with open(os.path.join(logs, "resume_2.out"), "w") as fh:
            fh.write("[RESUME] /x/fixture_arm/mol_B__6W63__primary/r2/06_md @ gpu01 2026-01-02T00:00:00\n")
            fh.write("[RESUME DONE] /x/fixture_arm/mol_B__6W63__primary/r2/06_md 2026-01-02T02:00:00\n")
        res2, prob2 = audit(tmp, "fixture_arm", 2, 6, verbose=False)
        if res2 and res2[0] == 5 and any("no terminal marker" in p for p in prob2):
            print("  [4] 抽掉 r3 的标记后判为 5/6 并指名该复本 —— 能够失败  PASS")
        else:
            print("  [4] 抽掉标记后仍得 %s / %s —— 闸门不能失败  FAIL" % (res2, prob2))
            ok = False

        # ⑤ 指向磁盘上不存在复本的标记必须被报出 —— 反向错误
        with open(os.path.join(logs, "stray.out"), "w") as fh:
            fh.write("[MD2 DONE] /x/fixture_arm/mol_C__6W63__primary/r1 branch=staged 2026-01-01T00:00:00\n")
        _, prob3 = audit(tmp, "fixture_arm", 2, 6, verbose=False)
        if any("marker for absent trajectory" in p for p in prob3):
            print("  [5] 指向不存在复本的标记被报出 —— 反向错误可察          PASS")
        else:
            print("  [5] 指向不存在复本的标记未被报出 —— FAIL")
            ok = False
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("=" * 70)
    print("自测 %s" % ("全部通过 —— 可以执行核查" if ok else "未通过 —— 不执行任何核查"))
    print("=" * 70)
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=ROOT_DEFAULT)
    ap.add_argument("--arm")
    ap.add_argument("--expect-systems", type=int)
    ap.add_argument("--expect-traj", type=int)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()

    if a.selftest:
        sys.exit(0 if selftest() else 3)

    if not a.arm:
        print("未指定 --arm。**未执行任何核查。**")
        sys.exit(3)

    if not selftest():
        print("\n自测未通过。**未执行任何核查。**")
        sys.exit(3)

    exp_sys, exp_traj = a.expect_systems, a.expect_traj
    if a.arm in REGISTRY:
        r_sys, r_traj = REGISTRY[a.arm]
        exp_sys = exp_sys if exp_sys is not None else r_sys
        exp_traj = exp_traj if exp_traj is not None else r_traj
    elif exp_sys is None or exp_traj is None:
        print("\n臂 %r 不在注册表内，须显式给出 --expect-systems 与 --expect-traj。" % a.arm)
        print("**未执行任何核查。**")
        sys.exit(3)

    print()
    print("=" * 70)
    print("C3（重做）· 完成标记，两种合法终止状态并计")
    print("=" * 70)
    res, prob = audit(a.root, a.arm, exp_sys, exp_traj)
    print()
    if prob:
        for p in prob:
            print("  PROBLEM: %s" % p)
        print()
        print("C3 FAIL —— 上列为所缺或所余者，逐条指名。")
        sys.exit(2)
    print("C3 PASS —— %d / %d 条轨迹各有至少一种合法终止标记。" % (res[0], res[1]))
    print("注：本工具只核第三项判据。规则 16 要求三项，帧数与结构产物数另核。")
    sys.exit(0)


if __name__ == "__main__":
    main()
