#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cohort31_provenance.py  v1.1  ·  2026-10-30  ·  READ-ONLY

Emit the generation-era provenance of the 31 simulated candidates -- the 31-cohort
equivalent of the 28-molecule table (pool name / BBS / original name / source run /
era).

WHY IT VALIDATES ITSELF FIRST
    The era attribution here is a SECOND IMPLEMENTATION of logic that already lives
    in Scoring/aud_al_bbs_recompute.py. A second implementation of an existing
    statistic is how aud_cross_site_occupancy.py v1 came to disagree with phase 4.
    So before emitting anything this script reproduces the audit's own published
    counts, and ABORTS if they do not match:

        attribution by field : primary_source 932 · run_types 674 · source_runs 0
        BBS cycles (parents) : cycle1 233 · cycle2 237 · cycle3 242 · cycle4 220

    Those are scope-level counts from the post-patch audit log of 2026-10-30
    (claude/AL_recompute_outcome_2026-10-29.md R1/G2). If this script's attribution
    reproduces them, it is equivalent to the audit's on the same population. If it
    does not, the table it would emit is not comparable with anything already
    recorded -- so it refuses to emit one.

WRITES NOTHING except the log you tee it to.

USAGE
    cd /mnt/hpc/home/wanghongmei/lang2
    python3 Scoring/cohort31_provenance.py 2>&1 | tee Scoring/cohort31_provenance_20261030.log
"""
import sys, csv, json, re, pathlib
from collections import Counter, defaultdict

ROOT = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()

SCOPE = ROOT / "results/00_scope/all_candidates_v3_start.csv"
BBSTB = ROOT / "results/04_docking/node04_bbs_table.csv"
IDENT = ROOT / "results/02_master_molecules/identity_manifests.jsonl"
MDDIR = ROOT / "md_solv_v3"

# identical to the post-patch CYCLE_ORDER in Scoring/aud_al_bbs_recompute.py
CYCLE_ORDER = ["pIC50_early", "pIC50_v1", "pIC50_v2",
               "BBS_cycle1", "BBS_cycle2", "BBS_cycle3", "BBS_cycle4"]
FIELDS = ["primary_source", "run_types", "source_runs"]   # same order as the audit

EXPECT_FIELD = {"primary_source": 932, "run_types": 674, "source_runs": 0}
EXPECT_CYCLE = {"BBS_cycle1": 233, "BBS_cycle2": 237, "BBS_cycle3": 242, "BBS_cycle4": 220}

POOL = re.compile(r"pool_mol_\d+")


def die(msg):
    print("\n[ABORT] " + msg)
    sys.exit(1)


def name_column(header, rows, where):
    """Find the column whose values contain pool_mol_NNN. Reported, not assumed.

    v1.0 used fullmatch and aborted on the docking table, whose identifiers carry a
    microstate suffix (pool_mol_001__m0). search is correct for both files: the scope
    holds bare parent names, the docking table holds suffixed microstate names.
    The abort message also named the scope file regardless of which file was scanned.
    """
    best, hits = None, -1
    for c in header:
        n = sum(1 for r in rows[:400] if POOL.search((r.get(c) or "").strip()))
        if n > hits:
            best, hits = c, n
    if hits <= 0:
        die("no column in %s holds pool_mol_NNN values; headers were: %s" % (where, header))
    print("    name column in %s: %r (%d/%d sampled rows contain pool_mol_NNN)"
          % (pathlib.Path(where).name, best, hits, min(400, len(rows))))
    return best


def main():
    print("cohort31_provenance.py v1.1 · READ-ONLY · root=%s\n" % ROOT)

    # ---------- 1. the simulated systems, from the trajectory tree ----------
    print("=" * 78)
    print("1 · SIMULATED SYSTEMS (from md_solv_v3/, not from a manifest)")
    print("=" * 78)
    if not MDDIR.is_dir():
        die("%s not found" % MDDIR)
    prim, sec, other = [], [], []
    for d in sorted(p.name for p in MDDIR.iterdir() if p.is_dir()):
        if not d.startswith("pool_mol_"):
            continue
        (prim if d.endswith("__primary") else sec if d.endswith("__secondary") else other).append(d)
    print("    primary systems   : %d" % len(prim))
    print("    secondary systems : %d  %s" % (len(sec), sec if sec else ""))
    if other:
        print("    OTHER (unexpected): %d  %s" % (len(other), other[:10]))
    want = sorted({POOL.search(d).group(0) for d in prim})
    print("    distinct primary molecules: %d" % len(want))
    if len(want) != 31:
        print("    [WARN] expected 31 primary molecules, found %d -- continuing, but the"
              % len(want))
        print("           manuscript's '31 molecules' needs reconciling with this count.")

    # ---------- 2. input headers, printed before use ----------
    print("\n" + "=" * 78)
    print("2 · INPUT FILES (headers printed before any column is used)")
    print("=" * 78)
    for p in (SCOPE, BBSTB, IDENT):
        print("    %-58s %s" % (p.relative_to(ROOT), "OK" if p.exists() else "MISSING"))
        if not p.exists():
            die("%s is required" % p)

    with open(SCOPE, newline="") as fh:
        rd = csv.DictReader(fh)
        scope_hdr = rd.fieldnames or []
        scope_rows = list(rd)
    print("\n    scope header (%d cols): %s" % (len(scope_hdr), scope_hdr))
    print("    scope rows: %d" % len(scope_rows))
    missing = [f for f in FIELDS if f not in scope_hdr]
    if missing:
        die("scope lacks attribution field(s) %s -- the audit used %s" % (missing, FIELDS))
    name_c = name_column(scope_hdr, scope_rows, SCOPE)

    with open(BBSTB, newline="") as fh:
        rd = csv.DictReader(fh)
        bbs_hdr = rd.fieldnames or []
        bbs_rows = list(rd)
    print("\n    bbs table header (%d cols): %s" % (len(bbs_hdr), bbs_hdr))
    print("    bbs table rows: %d" % len(bbs_rows))

    # ---------- 3. attribution, then SELF-VALIDATION ----------
    print("\n" + "=" * 78)
    print("3 · ERA ATTRIBUTION — validated against the audit before use")
    print("=" * 78)
    era, via = {}, Counter()
    for r in scope_rows:
        nm = (r.get(name_c) or "").strip()
        lab = used = None
        for f in FIELDS:
            raw = (r.get(f) or "")
            hits = [c for c in CYCLE_ORDER if c.lower() in raw.lower()]
            if hits:
                lab, used = hits[0], f
                break
        era[nm] = lab
        if used:
            via[used] += 1

    print("    attribution by field:")
    ok = True
    for f in FIELDS:
        exp = EXPECT_FIELD[f]
        good = via[f] == exp
        ok &= good
        print("      %-16s measured %5d   audit %5d   %s" % (f, via[f], exp, "OK" if good else "MISMATCH"))

    cyc = Counter(v for v in era.values() if v)
    print("    BBS-cycle parents:")
    for c, exp in EXPECT_CYCLE.items():
        good = cyc[c] == exp
        ok &= good
        print("      %-16s measured %5d   audit %5d   %s" % (c, cyc[c], exp, "OK" if good else "MISMATCH"))

    if not ok:
        die("this script's attribution does NOT reproduce the audit's counts.\n"
            "        Any table emitted from it would not be comparable with the recorded\n"
            "        figures, so none is emitted. Reconcile against\n"
            "        Scoring/aud_al_bbs_recompute.py before using this script.")
    print("\n    [OK] attribution reproduces the audit on all 7 checks.")

    # ---------- 4. BBS per molecule ----------
    bbs_col = None
    for cand in ("BBS", "bbs", "BBS_raw", "bbs_raw", "broad_spectrum_score", "Broad_Spectrum_Score"):
        if cand in bbs_hdr:
            bbs_col = cand
            break
    bname_c = name_column(bbs_hdr, bbs_rows, BBSTB) if bbs_rows else None
    best_bbs = {}
    if bbs_col and bname_c:
        print("    BBS column detected: %r" % bbs_col)
        for r in bbs_rows:
            nm = (r.get(bname_c) or "").strip()
            m = POOL.search(nm)
            if not m:
                continue
            try:
                v = float(r[bbs_col])
            except (TypeError, ValueError):
                continue
            k = m.group(0)
            if k not in best_bbs or v < best_bbs[k]:
                best_bbs[k] = v
        print("    BBS read for %d molecules (best = most negative microstate)" % len(best_bbs))
    else:
        print("    [WARN] no recognised BBS column in %s -- BBS column will read '-'" % BBSTB.name)
        print("           header was: %s" % bbs_hdr)

    # ---------- 5. the table ----------
    print("\n" + "=" * 78)
    print("4 · THE 31-COHORT PROVENANCE TABLE")
    print("=" * 78)
    scope_by = {(r.get(name_c) or "").strip(): r for r in scope_rows}
    print("%-16s %8s  %-14s  %-28s %s"
          % ("POOL", "BBS", "ERA", "primary_source", "run_types"))
    print("-" * 110)
    tally = Counter()
    absent = []
    for k in want:
        r = scope_by.get(k)
        if r is None:
            absent.append(k)
            print("%-16s %8s  %-14s  %s" % (k, "-", "NOT-IN-SCOPE", "—"))
            tally["NOT-IN-SCOPE"] += 1
            continue
        e = era.get(k) or "UNATTRIBUTED"
        tally[e] += 1
        b = "%8.2f" % best_bbs[k] if k in best_bbs else "       -"
        print("%-16s %8s  %-14s  %-28s %s"
              % (k, b, e, (r.get("primary_source") or "")[:28], (r.get("run_types") or "")[:40]))

    print("\n" + "=" * 78)
    print("5 · TALLY BY ERA")
    print("=" * 78)
    for e in CYCLE_ORDER + ["UNATTRIBUTED", "NOT-IN-SCOPE"]:
        if tally[e]:
            print("    %-16s %3d" % (e, tally[e]))
    print("    %-16s %3d" % ("TOTAL", sum(tally.values())))

    print("""
READ THIS AGAINST THE MANUSCRIPT
  Any candidate in pIC50_v2 (03VAE / 07VAE / 08VAE / 09VAE) means §2.1 must declare
  the unresolved loader question: those runs paired the original decoder (33,505-
  fragment vocabulary) with the deployed 10,173-fragment vocabulary. If the tally
  shows zero there, the question does not touch this cohort and §2.1 needs no note.

  UNATTRIBUTED or NOT-IN-SCOPE rows are findings, not noise -- a simulated candidate
  with no scope record would mean the cohort and the scope file disagree. Report the
  count; do not drop such rows.""")


if __name__ == "__main__":
    main()
