#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
audit_githubdyn.py -- read-only independent verification of the deposited tree
==============================================================================
    From the root of a fresh clone:
        python3 07_audit/Scoring/audit_githubdyn.py .

    If the tree sits in a subdirectory named `githubdyn`, which is the layout
    used while it was assembled, the argument may be omitted:
        python3 07_audit/Scoring/audit_githubdyn.py

**Read-only: there is no --apply and no write path anywhere in this file.**

**This script shares no code with the assembler that produced the tree.** That
is deliberate. Verifying a product with the producer's own functions
establishes only that the product agrees with itself. Everything here is
recomputed independently: SHA-256 is implemented locally, the manifest is
parsed locally, counts are taken locally, and notebooks are parsed locally.

Eight checks. Each guards against a failure mode that can quietly produce a
deposit that looks complete and is not.

  A  Manifest integrity. The paths in `MANIFEST.sha256` and the files actually
     present must agree in both directions, and every digest must match a hash
     recomputed from disk. What an assembler reports having copied and what is
     on disk are two different things.

  B  Notebook revisions. The six notebooks in the deposit must be the corrected
     revisions, identified by known SHA-256 prefixes. Were an earlier revision
     deposited, a mismatch between its prose and its code would travel with it.

  C  Arm-name residue and arm switching. The three cross-ortholog notebooks must
     no longer contain `ARMS = ['md_solv_v3']` in their markdown, and the three
     main-arm notebooks must set the arm to `md_solv_v3`.

  D  Per-frame series counts. The five `*_series/` directories must hold
     99 / 36 / 33 / 12 / 9 members in the tree, and the superseded
     `aud_cross_site_series` must be absent. A gate upstream governs what is
     generated; this check governs what was deposited.

  E  Per-replicate counts. `analysis.json` must appear 99 times for the main
     arm, 36 for the cross-ortholog arm and 12 for the reference arm.

  F  Exclusions documented. `EXCLUSIONS.md` must carry the section on
     directory-level exclusions, and that section must name
     `aud_cross_site_series`. Without it the README points at nothing.

  G  README declarations. `README.md` must carry the per-frame series section,
     the note on the directory name that never existed, the note on inherited
     registration text, and the note on arm switching.

  H  Whether git will in fact take the deposited `.npz` files. A `.gitignore`
     generated from a list of unwanted extensions contains `*.npz`, while this
     deposit deliberately includes per-frame series in that format. The tree can
     then be complete, the manifest can be complete, and check D can pass, while
     git takes none of them. **Check D inspects the filesystem; check H inspects
     git's own decision. Neither substitutes for the other.**
"""
import hashlib
import json
import os
import re
import sys

# Authoritative SHA-256 prefixes (first 12 hex digits) of the corrected
# notebooks, as reported by each correction run.
EXPECT_NB = {
    "phase2_cross_V3.ipynb":        ("f0e584e7b326", "arm-name residue corrected (3 sites)"),
    "phase3_cross_V3.ipynb":        ("0d74d7c797ee", "arm-name residue corrected (1 site)"),
    "phase4_cross_V3.ipynb":        ("83a9a5c02be0", "arm-name residue corrected (1 site)"),
    "phase3_mmpbsa_solv_V3.ipynb":  ("6ecc05f00ca6", "option A: ARM restored to md_solv_v3"),
    "phase4_analysis_solv_V3.ipynb": ("240bf0b9c9cc", "option A: ARM restored to md_solv_v3"),
    "phase2_pbc_correct_V3.ipynb":  ("2a9860a57783", "unchanged (prose and code already agreed)"),
}

SERIES_EXPECT = {
    "aud_main_site_occupancy_series": 99,
    "aud_cross_site_occupancy_series": 36,
    "aud_main_site_occupancy_his163_series": 33,
    "aud_b81_scope9_series": 9,
    # Reference arm, added 2026-11-06: phase 4 was run for 4 systems x 3 replicates.
    "aud_ref_site_occupancy_series": 12,
}
FORBIDDEN_DIR = "aud_cross_site_series"      # superseded 30-member generation
ANALYSIS_EXPECT = {"step6e_solv_v3": 99, "step6e_cross_v3": 36,
    # Reference arm, 12 replicates. While this entry was absent the check
    # covered only 99 and 36, so the loss of those 12 would have passed silently.
    "step6e_solv_v3_ref": 12,}

RESIDUE = "ARMS = ['md_solv_v3']"
CROSS_NB = ("phase2_cross_V3.ipynb", "phase3_cross_V3.ipynb", "phase4_cross_V3.ipynb")
MAIN_NB = ("phase2_pbc_correct_V3.ipynb", "phase3_mmpbsa_solv_V3.ipynb",
           "phase4_analysis_solv_V3.ipynb")

README_NEEDS = [
    ("Per-frame d_dyad series", "per-frame series section"),
    ("md_cross_solv_v3", "note on the directory name that never existed"),
    ("Inherited registration text", "note on inherited registration text in cross-arm sidecars"),
    ("ARM = EXPECTED_ARM", "note on arm switching in the main-arm notebooks"),
]

FAILS = []


def fail(tag, msg):
    FAILS.append((tag, msg))
    print("    **FAIL** %s" % msg)


def ok(msg):
    print("    PASS     %s" % msg)


def h256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def tree_files(root):
    out = {}
    for dp, dns, fns in os.walk(root):
        if ".git" in dns:
            dns.remove(".git")
        for fn in fns:
            p = os.path.join(dp, fn)
            out[os.path.relpath(p, root).replace(os.sep, "/")] = p
    return out


def nb_cells(p):
    nb = json.load(open(p, encoding="utf-8"))
    return [(c.get("cell_type", "?"),
             "".join(c.get("source", "")) if isinstance(c.get("source"), list)
             else str(c.get("source", "")))
            for c in nb.get("cells", [])]


def check_a(root, files):
    print("\nA. Manifest and tree agree in both directions, digest by digest")
    mp = os.path.join(root, "MANIFEST.sha256")
    if not os.path.isfile(mp):
        return fail("A", "MANIFEST.sha256 is absent")
    # Manifest line format, established by reading the assembler rather than
    # assumed:
    #     <sha256>  <path within the deposit>  <-  <source path>
    # A pattern that takes the whole remainder of the line as the path also
    # captures the `  <-  <source>` tail, so every key fails to match and a
    # sound deposit is reported as corrupt throughout. The format can be read;
    # it should not be asserted.
    listed, unparsed = {}, 0
    for ln in open(mp, encoding="utf-8", errors="replace"):
        ln = ln.rstrip("\n")
        if not ln.strip() or ln.lstrip().startswith("#"):
            continue
        m = re.match(r"\s*([0-9a-f]{64})\s\s+(.*)$", ln)
        if not m:
            unparsed += 1
            continue
        dest = m.group(2).split("  <-  ", 1)[0].strip()
        listed[dest.replace("\\", "/").lstrip("./")] = m.group(1)
    if not listed:
        return fail("A", "MANIFEST.sha256 parsed to 0 entries (format differs from expectation)")
    if unparsed:
        fail("A", "%d line(s) in the manifest could not be parsed" % unparsed)
    # Two numbers are printed separately, each labelled for what it is. A
    # correct number under a misleading label is a defect in its own right,
    # because the reader sees only the label.
    #
    # LICENSE belongs in this set: it is one of the files deposited without a
    # manifest entry.
    meta = {"MANIFEST.sha256", "EXCLUSIONS.md", "NOT_FOUND.md", "PROVENANCE.md",
            "README.md", ".gitignore", ".gitattributes", "TOO_LARGE_FOR_GIT.md",
            "LICENSE"}

    # `.ipynb_checkpoints/` is reported separately and kept out of the two-way
    # comparison, so that one stray checkpoint cannot bury the real comparison.
    # It should not be in the deposit tree at all.
    ckpt = sorted(k for k in files if ".ipynb_checkpoints/" in k + "/")
    if ckpt:
        gi = files.get(".gitignore")
        covered = (".ipynb_checkpoints" in
                   open(gi, encoding="utf-8", errors="replace").read()) if gi else False
        fail("A", "the deposit tree holds %d `.ipynb_checkpoints/` file(s) (e.g. `%s`); "
                  ".gitignore %s -- %s"
             % (len(ckpt), ckpt[0],
                "covers that directory" if covered else "**does not cover that directory**",
                "git will not take it, but tree and manifest disagree; remove it" if covered
                else "**git will take it**; it must be removed first"))

    actual = {k for k in files if k not in meta and
              not k.endswith("IDENTICAL.sha256") and
              ".ipynb_checkpoints/" not in k + "/"}
    print("    manifest lists %d; the tree holds %d files, of which %d are compared"
          " (the other %d are manifest, documentation and de-duplication records%s)"
          % (len(listed), len(files), len(actual), len(files) - len(actual),
             "" if not ckpt else ", plus %d checkpoint file(s)" % len(ckpt)))
    only_tree = sorted(actual - set(listed))
    only_list = sorted(set(listed) - actual)
    if only_tree:
        fail("A", "%d file(s) in the tree are absent from the manifest, e.g. %s"
             % (len(only_tree), only_tree[:3]))
    if only_list:
        fail("A", "%d manifest entries are absent from the tree, e.g. %s"
             % (len(only_list), only_list[:3]))
    bad = [k for k in sorted(set(listed) & actual) if h256(files[k]) != listed[k]]
    if bad:
        fail("A", "%d digest mismatch(es), e.g. %s" % (len(bad), bad[:3]))
    elif not only_tree and not only_list:
        ok("agreement in both directions; %d files match digest by digest" % len(actual))


def check_b(root, files):
    print("\nB. Whether the deposited notebooks are the corrected revisions")
    byname = {}
    for rel, p in files.items():
        if rel.endswith(".ipynb"):
            byname.setdefault(os.path.basename(rel), []).append((rel, p))
    for base, (want, why) in sorted(EXPECT_NB.items()):
        got = byname.get(base)
        if not got:
            fail("B", "%s is not in the deposit (%s)" % (base, why))
            continue
        for rel, p in got:
            h = h256(p)[:12]
            if h == want:
                ok("%-32s %s  %s" % (base, h, why))
            else:
                fail("B", "the sha256 of %s is %s, expected %s -- "
                     "**the deposit does not hold the corrected revision**"
                     % (rel, h, want))


def check_c(root, files):
    print("\nC. Arm-name residue and arm switching")
    byname = {}
    for rel, p in files.items():
        if rel.endswith(".ipynb"):
            byname.setdefault(os.path.basename(rel), []).append((rel, p))
    for base in CROSS_NB:
        for rel, p in byname.get(base, []):
            md = [t for k, t in nb_cells(p) if k == "markdown"]
            n = sum(t.count(RESIDUE) for t in md)
            (ok if n == 0 else lambda m: fail("C", m))(
                "%-28s %d main-arm residue(s) in markdown" % (base, n))
    for base in MAIN_NB:
        for rel, p in byname.get(base, []):
            lines = [ln for k, t in nb_cells(p) if k == "code"
                     for ln in t.splitlines() if ln.lstrip().startswith("ARM")
                     and "=" in ln and "ARM_N" not in ln]
            arm = None
            for ln in lines:
                m = re.search(r'=\s*"([\w]+)"', ln)
                if m:
                    arm = m.group(1); break
            if arm == "md_solv_v3":
                ok("%-28s switch line reads md_solv_v3 (option A)" % base)
            else:
                fail("C", "the switch line of %s reads %r, expected md_solv_v3" % (base, arm))


def check_d(root, files):
    print("\nD. Per-frame series counts in the tree")
    cnt = {}
    for rel in files:
        if not rel.endswith(".npz"):
            continue
        d = os.path.basename(os.path.dirname(rel))
        cnt[d] = cnt.get(d, 0) + 1
    for d, want in sorted(SERIES_EXPECT.items()):
        got = cnt.get(d, 0)
        (ok if got == want else lambda m: fail("D", m))(
            "%-44s %d present (registered %d)" % (d, got, want))
    if FORBIDDEN_DIR in cnt:
        fail("D", "**the superseded 30-member `%s/` entered the deposit, %d files**"
             % (FORBIDDEN_DIR, cnt[FORBIDDEN_DIR]))
    else:
        ok("the superseded 30-member `%s/` is absent from the deposit" % FORBIDDEN_DIR)
    stray = {d: n for d, n in cnt.items() if d not in SERIES_EXPECT}
    if stray:
        fail("D", "npz files under unregistered directories: %s" % stray)


def check_e(root, files):
    print("\nE. Per-replicate analysis.json counts")
    cnt = {}
    for rel in files:
        if os.path.basename(rel) != "analysis.json":
            continue
        for seg in rel.split("/"):
            if seg in ANALYSIS_EXPECT:
                cnt[seg] = cnt.get(seg, 0) + 1
    for seg, want in sorted(ANALYSIS_EXPECT.items()):
        got = cnt.get(seg, 0)
        (ok if got == want else lambda m: fail("E", m))(
            "%-24s %d present (registered %d)" % (seg, got, want))


def check_f(root, files):
    print("\nF. Whether the reasons for exclusion travel with the deposit")
    p = files.get("EXCLUSIONS.md")
    if not p:
        return fail("F", "EXCLUSIONS.md is not in the deposit")
    s = open(p, encoding="utf-8", errors="replace").read()
    (ok if "Directory-level exclusions on record" in s else lambda m: fail("F", m))(
        "directory-level exclusions section in EXCLUSIONS.md")
    if FORBIDDEN_DIR in s:
        ok("it names `%s`, so the reason travels with the deposit" % FORBIDDEN_DIR)
    else:
        fail("F", "**the reason for excluding `%s` is nowhere in the deposit** -- "
             "the pointer in the README leads nowhere" % FORBIDDEN_DIR)


def check_g(root, files):
    print("\nG. The four declarations in the README")
    p = files.get("README.md")
    if not p:
        return fail("G", "README.md is not in the deposit")
    s = open(p, encoding="utf-8", errors="replace").read()
    for needle, label in README_NEEDS:
        (ok if needle in s else lambda m: fail("G", m))("%s (searched %r)" % (label, needle))


def check_h(root, files):
    """H. Whether git will in fact take them.

    This check exists because a `.gitignore` generated from a list of unwanted
    extensions contains `*.npz`, while the deposit deliberately includes
    per-frame series in that format. The files are then present in the tree,
    present in the manifest, and counted by check D, while git takes none of
    them. **Check D inspects the filesystem; check H inspects git's own
    decision. Neither substitutes for the other.**

    `git check-ignore` is authoritative and is tried first. Where git is
    unavailable the check falls back to reading `.gitignore` as text and says
    plainly that this is weaker evidence.
    """
    print("\nH. Whether git will take the deliberately deposited npz files"
          " (check D inspects the filesystem, this one inspects git)")
    npz = sorted(k for k in files if k.endswith(".npz")
                 and os.path.basename(os.path.dirname(k)) in SERIES_EXPECT)
    if not npz:
        return fail("H", "no npz files under any registered directory were found in the tree")
    try:
        import subprocess
        # **Do not use `git -C`.** That option arrived in git 1.8.5, and the
        # cluster this project runs on carries an older git which reports
        # `Unknown option: -C`; the authoritative half of this check would then
        # be lost silently, leaving only the weaker textual evidence. `cwd=`
        # works on every version. `--no-index` is likewise absent from older
        # versions, so it is tried first and then without.
        def run(extra):
            return subprocess.run(["git", "check-ignore"] + extra + npz,
                                  cwd=root, capture_output=True, text=True,
                                  timeout=120)
        r = run(["--no-index"])
        if r.returncode not in (0, 1):
            r = run([])
        # check-ignore prints each ignored path; exit code 0 means at least one
        # path is ignored, 1 means none is.
        ignored = [x for x in r.stdout.splitlines() if x.strip()]
        if r.returncode not in (0, 1):
            raise RuntimeError(r.stderr.strip()[:200] or "git returned %d" % r.returncode)
        if ignored:
            fail("H", "**git would ignore %d of %d deposited npz files** -- they would be "
                      "missing from what is published. Example: `%s`. Negative patterns "
                      "`!**/<dir>/*.npz` are required in `.gitignore` for these directories"
                 % (len(ignored), len(npz), ignored[0]))
        else:
            ok("git check-ignore: all %d registered npz files are taken (authoritative)" % len(npz))
    except Exception as e:
        print("    (git check-ignore unavailable: %s -- falling back to reading"
              " .gitignore as text, which is weaker evidence)" % e)
        gi = files.get(".gitignore")
        if not gi:
            return fail("H", ".gitignore is not in the deposit and git is unavailable:"
                             " no determination can be made")
        s = open(gi, encoding="utf-8", errors="replace").read()
        if "*.npz" not in s:
            ok(".gitignore does not exclude `*.npz`")
            return
        miss = [d for d in SERIES_EXPECT
                if ("!**/%s/*.npz" % d) not in s and ("!%s" % d) not in s]
        if miss:
            fail("H", "`.gitignore` contains `*.npz` and these directories have no"
                      " negative pattern: %s" % miss)
        else:
            ok("`*.npz` is excluded, but every registered directory has a negative pattern")


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else "githubdyn"
    root = os.path.abspath(root)
    if not os.path.isdir(root):
        sys.exit("deposit directory does not exist: %s" % root)
    print("deposit tree : %s" % root)
    print("mode         : **read-only**; writes no file; shares no code with the assembler")
    files = tree_files(root)
    for fn in (check_a, check_b, check_c, check_d, check_e, check_f,
               check_g, check_h):
        fn(root, files)
    print("\n" + "=" * 70)
    if FAILS:
        print("**%d check(s) failed** -- listed below. Do not run git until they are resolved."
              % len(FAILS))
        for tag, msg in FAILS:
            print("  [%s] %s" % (tag, msg))
        sys.exit(1)
    print("All eight checks passed.")


if __name__ == "__main__":
    main()
