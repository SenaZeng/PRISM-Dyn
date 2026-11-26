# -*- coding: utf-8 -*-
"""aud_frame_census.py · 规则 16 的三项判据，按臂普查
   Rule-16 completion census, per arm

    自检 / self-test :  python Scoring/aud_frame_census.py --selftest
    主臂 / main arm  :  python Scoring/aud_frame_census.py --arm md_solv_v3
    参考臂           :  python Scoring/aud_frame_census.py --arm md_solv_v3_ref
    跨臂（须显式给数）:  python Scoring/aud_frame_census.py --arm cross_solv_v3 \
                            --expect-systems 12 --expect-traj 36

READ-ONLY。不写、不建、不提交。退出码 0 = 三项判据皆通过；2 = 任一项失败。

═══ 与被检之物不共享代码（原则 3）═══
本脚本【不 import】core/node05_md_prep.py，亦不读取任何驱动或 notebook 所产出的
清单、freeze 或 record。它自行走目录树、自行解析 `gmx check` 的输出、自行计数。
故生产侧与核验侧的缺陷不可能同时不可见。

═══ 三项判据（规则 16）═══
  C1 帧数   每条轨迹的 `gmx check -f md.xtc` 须报 10001 帧，且末帧时刻 100000.000 ps。
            缺 `Last frame` 行者判为不完整——被放弃的那三条残片正是此形态。
  C2 结构产物数  phase2 产物 `md_fit.xtc` 的数量须等于轨迹数。
  C3 完成标记数  `<arm>/logs/` 之中 `[MD2 DONE]` 的【出现次数】须等于轨迹数。
            **计出现次数，不计含该文本的文件数。** 后者会把运行器脚本自身算进来——
            跨臂因此返回 66 而完成数为 30。

═══ 三者不可互相替代 ═══
  一个文件可以存在、大小貌似合理，而只载有预定采样的一小部分（跨臂被放弃的三条为
  1752/1742/1790 帧对 10001）。反之，完成标记可因日志轮替而丢失，数据却完好
  （主臂即如此：`pool_mol_001__6W63` 三条完整复本中仅一条的标记尚存）。
  故 C1 取自数据，C2、C3 为佐证；任一项缺失须声明所缺者，不得以其余两项顶替。

═══ D28 ═══
  通配一律锚定双分隔符：`<arm>/*/*/06_md`，绝不以前缀放宽。
  默认排除路径中含 `abandon` 者，并【报告】排除数，不静默丢弃。
"""
import argparse, glob, os, re, subprocess, sys

EXPECT_FRAMES = 10001
EXPECT_LAST_PS = 100000.0
PS_TOL = 0.001
# 已登记的臂；未登记者须显式给数（与 ARM_N 同一纪律）
REGISTRY = {"md_solv_v3": (33, 99), "md_solv_v3_ref": (4, 12)}

RE_COORDS = re.compile(r"^\s*Coords\s+(\d+)\s+(\S+)", re.M)
RE_LAST = re.compile(r"Last frame\s+(\d+)\s+time\s+([0-9.]+)", re.M)


def parse_gmx_check(text):
    """gmx check 的输出 -> {nframes, timestep, last_index, last_ps} 或 None 值。
    纯函数，不触碰文件系统，故可由自检以构造输入驱动。"""
    out = {"nframes": None, "timestep": None, "last_index": None, "last_ps": None}
    m = RE_COORDS.search(text or "")
    if m:
        out["nframes"] = int(m.group(1))
        out["timestep"] = m.group(2)
    m = RE_LAST.search(text or "")
    if m:
        out["last_index"] = int(m.group(1))
        out["last_ps"] = float(m.group(2))
    return out


def verdict(p):
    """解析结果 -> (ok, 原因)。亦为纯函数。"""
    if p["nframes"] is None:
        return False, "no Coords line — gmx check produced no frame count"
    if p["last_ps"] is None:
        return False, "no 'Last frame' line — truncated trajectory (%d frames read)" % p["nframes"]
    if p["nframes"] != EXPECT_FRAMES:
        return False, "%d frames, expected %d" % (p["nframes"], EXPECT_FRAMES)
    if abs(p["last_ps"] - EXPECT_LAST_PS) > PS_TOL:
        return False, "last frame at %.3f ps, expected %.3f" % (p["last_ps"], EXPECT_LAST_PS)
    return True, "%d frames, last %.3f ps" % (p["nframes"], p["last_ps"])


# ───────────────────────────────────────────────────────────── self-test
def selftest():
    fails = []

    def chk(name, got, want):
        ok = got == want
        print("  %-54s %s" % (name, "OK" if ok else "FAIL got=%r want=%r" % (got, want)))
        if not ok:
            fails.append(name)

    GOOD = ("Reading frame       0 time    0.000\n"
            "Reading frame   10000 time 100000.000\n"
            "Item        #frames Timestep (ps)\n"
            "        Step         10001    10\n"
            "        Coords       10001    10\n"
            "        Box          10001    10\n"
            "Last frame        10000 time 100000.000\n")
    TRUNC = ("Reading frame       0 time    0.000\n"
             "Reading frame    1700 time 17000.000\n"
             "        Coords        1752    10\n")
    SHORT = GOOD.replace("10001    10", " 9001    10").replace("time 100000.000", "time  90000.000")
    GARBAGE = "Fatal error:\nCannot open file md.xtc\n"
    EMPTY = ""

    print("A · parser and verdict, on constructed gmx-check output")
    chk("complete trajectory accepted", verdict(parse_gmx_check(GOOD))[0], True)
    chk("truncated (no Last frame) refused", verdict(parse_gmx_check(TRUNC))[0], False)
    chk("  and the reason names truncation",
        "truncated" in verdict(parse_gmx_check(TRUNC))[1], True)
    chk("wrong frame count refused", verdict(parse_gmx_check(SHORT))[0], False)
    chk("gmx error output refused", verdict(parse_gmx_check(GARBAGE))[0], False)
    chk("empty output refused", verdict(parse_gmx_check(EMPTY))[0], False)
    chk("frames parsed from complete output", parse_gmx_check(GOOD)["nframes"], 10001)
    chk("frames parsed from truncated output", parse_gmx_check(TRUNC)["nframes"], 1752)
    # 一条右数值、错时刻的轨迹：帧数对而末帧时刻不对，须拒绝
    WRONGT = GOOD.replace("Last frame        10000 time 100000.000",
                          "Last frame        10000 time  50000.000")
    chk("right frame count but wrong last time refused", verdict(parse_gmx_check(WRONGT))[0], False)

    print("\nB · tree walk and the two counting criteria, on a synthetic arm")
    import shutil, tempfile
    root = tempfile.mkdtemp()
    try:
        arm = os.path.join(root, "fake_arm")
        os.makedirs(os.path.join(arm, "logs"))
        for s in ("sysA", "sysB"):
            for r in ("r1", "r2", "r3"):
                d = os.path.join(arm, s, r, "06_md")
                os.makedirs(d)
                open(os.path.join(d, "md.xtc"), "w").write("x")
                open(os.path.join(d, "md_fit.xtc"), "w").write("x")
        # 6 条轨迹；logs 中写 6 次标记，另在一个 run_md2.sh 中也写一次（须不被计入）
        open(os.path.join(arm, "logs", "prod_1.out"), "w").write("[MD2 DONE] a\n" * 6)
        open(os.path.join(arm, "sysA", "r1", "run_md2.sh"), "w").write('echo "[MD2 DONE]"\n')
        chk("trajectory count from anchored glob", count_traj(arm)[0], 6)
        chk("md_fit count equals trajectory count", count_fit(arm), 6)
        chk("markers counted as events in logs/ only", count_markers(arm), 6)
        os.remove(os.path.join(arm, "sysB", "r3", "06_md", "md_fit.xtc"))
        chk("a missing md_fit is detected", count_fit(arm), 5)
        open(os.path.join(arm, "logs", "prod_2.out"), "w").write("[MD2 DONE] extra\n")
        chk("an extra marker is detected", count_markers(arm), 7)
        # localisation: markers that name their run dir are attributed; others are not
        open(os.path.join(arm, "logs", "prod_3.out"), "w").write(
            "[MD2 DONE] /x/sysA__6W63__primary/r1 branch=staged 2026\n"
            "[MD2 DONE] /x/sysA__6W63__primary/r2 branch=staged 2026\n"
            "[MD2 DONE] /x/sysQ__6W63__primary/r1 branch=staged 2026\n")
        by, unatt = markers_by_system(arm)
        chk("markers attributed to the system they name", sorted(by.get("sysA__6W63__primary", [])),
            ["r1", "r2"])
        chk("  a short system is visible as 2 of 3", len(by.get("sysA__6W63__primary", [])), 2)
        chk("  a marker naming a foreign system is kept separate",
            "sysQ__6W63__primary" in by, True)
        chk("  markers with no run dir are counted as unattributed", unatt, 7)
        ab = os.path.join(arm, "_abandoned_x", "sysC", "r1", "06_md")
        os.makedirs(ab)
        open(os.path.join(ab, "md.xtc"), "w").write("x")
        n, ex, lst = count_traj(arm)
        chk("abandoned path excluded from the arm count", n, 6)
        chk("  and counted by the unbounded walk, not assumed absent", ex, 1)
        chk("  and enumerated so it can be inspected", len(lst), 1)
    finally:
        shutil.rmtree(root)

    print("\n%s" % ("SELFTEST PASSED" if not fails else "SELFTEST FAILED: %r" % fails))
    return 0 if not fails else 2


# ───────────────────────────────────────────────────── counting (own logic)
def count_traj(arm):
    """(锚定计数, 树中别处的 md.xtc 数)。

    BF1 的重现，由自检在接触真实数据之前捕获：锚定通配 `<arm>/*/*/06_md/md.xtc`
    【按深度】就已排除嵌套于 `<arm>/_abandoned_*/...` 之下者，故对 hits 再作
    "abandon" 过滤恒得 0，而那一行报告会宣称「已排除 0 条」——读者会据此以为树中
    别无他物。本函数因此改以【两次独立的遍历】作答：锚定通配给出该臂的计数，
    无界遍历给出树中别处之数，二者之差即须被报告者。
    与 2026-10-30 上午 recon v1.0 中 C 节与 D 节深度不一致的那处缺陷同类。"""
    anchored = [p for p in glob.glob(os.path.join(arm, "*", "*", "06_md", "md.xtc"))
                if "abandon" not in p.lower()]
    everywhere = []
    for dirpath, _, files in os.walk(arm):
        if "md.xtc" in files:
            everywhere.append(os.path.join(dirpath, "md.xtc"))
    aset = set(os.path.normpath(x) for x in anchored)
    elsewhere = [x for x in everywhere if os.path.normpath(x) not in aset]
    return len(anchored), len(elsewhere), sorted(elsewhere)


def count_fit(arm):
    return len([p for p in glob.glob(os.path.join(arm, "*", "*", "06_md", "md_fit.xtc"))
                if "abandon" not in p.lower()])


RE_MARKER_RUN = re.compile(r"\[MD2 DONE\][^\n]*?/([A-Za-z0-9_]+__[0-9A-Za-z]+__[a-z]+)/(r\d+)")


def count_markers(arm):
    n = 0
    for dirpath, _, files in os.walk(os.path.join(arm, "logs")):
        for f in files:
            try:
                with open(os.path.join(dirpath, f), errors="ignore") as fh:
                    n += fh.read().count("[MD2 DONE]")
            except OSError:
                pass
    return n


def markers_by_system(arm):
    """把每条完成标记归到它所报的 <system>/<replica>，使 C3 的差额可定位而非仅可计数。
    标记行形如 `[MD2 DONE] <run dir> branch=... <date>`，故运行目录可自标记自身读出——
    不依赖任何清单或 freeze（原则 3）。返回 {system: set(replicas)} 与无法归属的行数。"""
    by, unattributed = {}, 0
    for dirpath, _, files in os.walk(os.path.join(arm, "logs")):
        for f in files:
            try:
                txt = open(os.path.join(dirpath, f), errors="ignore").read()
            except OSError:
                continue
            for line in txt.splitlines():
                if "[MD2 DONE]" not in line:
                    continue
                m = RE_MARKER_RUN.search(line)
                if m:
                    by.setdefault(m.group(1), set()).add(m.group(2))
                else:
                    unattributed += 1
    return by, unattributed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/mnt/hpc/home/wanghongmei/lang2")
    ap.add_argument("--arm")
    ap.add_argument("--expect-systems", type=int)
    ap.add_argument("--expect-traj", type=int)
    ap.add_argument("--limit", type=int, help="只查前 N 条（冒烟用，不作为普查）")
    ap.add_argument("--gmx", default="gmx")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    if not a.arm:
        return ap.error("--arm is required (or --selftest)")

    arm = os.path.join(a.root, a.arm)
    if not os.path.isdir(arm):
        sys.exit("REFUSED: not a directory: %s" % arm)
    if a.arm in REGISTRY:
        exp_sys, exp_traj = REGISTRY[a.arm]
        src = "registered"
        if (a.expect_systems and a.expect_systems != exp_sys) or \
           (a.expect_traj and a.expect_traj != exp_traj):
            sys.exit("REFUSED: %s is registered as %d/%d; refusing a conflicting override"
                     % (a.arm, exp_sys, exp_traj))
    else:
        if a.expect_systems is None or a.expect_traj is None:
            sys.exit("REFUSED: %r is not a registered arm. Supply --expect-systems and "
                     "--expect-traj explicitly; this script will not guess a count."
                     % a.arm)
        exp_sys, exp_traj, src = a.expect_systems, a.expect_traj, "supplied on the command line"

    def hr(s):
        print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)

    hr("aud_frame_census.py · arm=%s · expected %d systems / %d trajectories (%s)"
       % (a.arm, exp_sys, exp_traj, src))
    print("READ-ONLY. gmx check reads each trajectory in full, so this is minutes, not seconds.")

    hr("C2 · structural products (phase2 md_fit.xtc)")
    n_traj, n_excl, excl_list = count_traj(arm)
    n_fit = count_fit(arm)
    n_sys = len(set(os.path.normpath(p).split(os.sep)[-4]
                    for p in glob.glob(os.path.join(arm, "*", "*", "06_md", "md.xtc"))
                    if "abandon" not in p.lower()))
    print("  md.xtc      = %-5d (expected %d)" % (n_traj, exp_traj))
    print("  systems     = %-5d (expected %d)" % (n_sys, exp_sys))
    print("  md_fit.xtc  = %-5d (expected %d)" % (n_fit, exp_traj))
    print("  md.xtc elsewhere in the tree, NOT at %s/*/*/06_md/ : %d" % (a.arm, n_excl))
    for x in excl_list[:12]:
        print("      %s" % os.path.relpath(x, a.root))
    if len(excl_list) > 12:
        print("      ... and %d more" % (len(excl_list) - 12))
    if n_excl:
        print("  （这些不计入本臂的三项判据。它们由无界遍历数出，而非由锚定通配"
              "——后者按深度就看不到它们。）")
    c2 = (n_traj == exp_traj and n_sys == exp_sys and n_fit == exp_traj)
    print("  C2 %s" % ("PASS" if c2 else "FAIL"))

    hr("C3 · completion markers, counted as EVENTS in logs/ only")
    n_mk = count_markers(arm)
    print("  [MD2 DONE] occurrences under %s/logs = %d (expected %d)" % (a.arm, n_mk, exp_traj))
    c3 = (n_mk == exp_traj)
    print("  C3 %s" % ("PASS" if c3 else "FAIL"))
    if not c3:
        by, unatt = markers_by_system(arm)
        sysdirs = sorted(set(os.path.normpath(p).split(os.sep)[-4]
                             for p in glob.glob(os.path.join(arm, "*", "*", "06_md", "md.xtc"))
                             if "abandon" not in p.lower()))
        print("\n  其差额的归属（自标记行自身读出，不依赖任何清单）:")
        print("  Localisation of the shortfall, read from the marker lines themselves:")
        short = []
        for sd in sysdirs:
            reps = by.get(sd, set())
            if len(reps) != 3:
                short.append((sd, sorted(reps)))
        if short:
            for sd, reps in short:
                print("    %-44s markers for %d of 3 replicas: %s"
                      % (sd, len(reps), ",".join(reps) if reps else "none"))
            print("    —— 合计缺 %d 条标记，分布于 %d 个体系"
                  % (exp_traj - n_mk, len(short)))
        else:
            print("    每个体系皆有 3 条标记——差额不在体系层面，须查重复或无法归属的行")
        extra = sorted(k for k in by if k not in sysdirs)
        if extra:
            print("    标记指向不在本臂的体系: %r" % extra)
        if unatt:
            print("    无法归属到 <system>/<replica> 的标记行: %d" % unatt)
        print("\n  **这是日志留存的缺口，不是数据的缺口**——除非 C1 在同一体系上也失败。")
        print("  规则 16：须在稿中声明所缺的判据，不得以其余两项顶替。")

    hr("C1 · frame count, read from each trajectory")
    paths = sorted(p for p in glob.glob(os.path.join(arm, "*", "*", "06_md", "md.xtc"))
                   if "abandon" not in p.lower())
    if a.limit:
        paths = paths[:a.limit]
        print("  --limit %d given: this is a smoke run, NOT a census\n" % a.limit)
    bad, good = [], 0
    for i, p in enumerate(paths, 1):
        rel = os.path.relpath(p, a.root)
        try:
            r = subprocess.run([a.gmx, "check", "-f", p], capture_output=True, text=True,
                               timeout=3600)
            txt = (r.stdout or "") + (r.stderr or "")
        except Exception as e:
            txt = "EXCEPTION %s" % e
        ok, why = verdict(parse_gmx_check(txt))
        if ok:
            good += 1
        else:
            bad.append((rel, why))
        print("  [%3d/%3d] %-62s %s" % (i, len(paths), rel[-62:], "ok" if ok else "FAIL: " + why))
        sys.stdout.flush()
    if a.limit:
        c1_state = "SMOKE ONLY"
    elif bad or len(paths) != exp_traj:
        c1_state = "FAIL"
    else:
        c1_state = "PASS"
    c1 = (c1_state == "PASS")
    print("\n  complete: %d / %d" % (good, len(paths)))
    if bad:
        print("  FAILURES:")
        for rel, why in bad:
            print("    %-66s %s" % (rel, why))
    print("  C1 %s" % c1_state)

    hr("VERDICT (rule 16 requires all three)")
    print("  %-22s %s" % ("C1 frames", c1_state))
    print("  %-22s %s" % ("C2 structural count", "PASS" if c2 else "FAIL"))
    print("  %-22s %s" % ("C3 markers", "PASS" if c3 else "FAIL"))
    if a.limit:
        print("\n*** NO CENSUS WAS PERFORMED. --limit %d was given, so C1 examined %d of %d\n"
              "    trajectories. A smoke run must not be reported, cited or recorded as a\n"
              "    census, and this exit code (3) is distinct from a failed census (2) for\n"
              "    that reason. Re-run without --limit for the census." % (a.limit, len(paths), exp_traj))
        return 3
    allp = c1 and c2 and c3
    print("\n%s" % ("ALL THREE PASS — the arm's completeness rests on three criteria."
                    if allp else
                    "NOT ALL THREE PASS — state in the manuscript which criterion is missing "
                    "(rule 16); do not substitute the others for it."))
    return 0 if allp else 2


if __name__ == "__main__":
    sys.exit(main())
