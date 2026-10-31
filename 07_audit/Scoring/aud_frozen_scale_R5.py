#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""chk_frozen_scale_R5.py — R5 自检：冻结标尺与已发表 Z 列的独立核验
                            Independent verification of the frozen COMPASS scale

READ-ONLY. Writes nothing. Exit 0 = all checks pass, 2 = a check failed.

**Why an independent script rather than re-running the notebook.**
HANDOFF §1.7: an independent checker shares no code with what it checks. This file
reimplements the R5 check of the Scope B registration from the registered definitions
alone; it imports nothing from reference_panel_COMPASS_benchmark_V3.ipynb. If the two
agree, that agreement means something. If they disagree, one of them is wrong and the
disagreement is the finding.

**What it checks, in order**
  A · the frozen scale recomputed from step8d_input_combined_SOLV_V3.csv (31 rows)
      equals the values registered in R1, to < 1e-5
  B · element-wise: Z_i = (x_i - median31)/(1.4826*MAD31 + 1e-9) reproduces every
      published Z_* column to < 1e-6
  B'· reverse assertion: bare MAD (no 1.4826) MUST fail check B. If it passes,
      check B is not sensitive to the MAD convention and is therefore worthless.
  C · COMPASS_v3 = sum(W_i * Z_i) reproduces the published score (tolerance 0.05;
      the residual is the 4-dp rounding of the weights stored in the CSV)
  D · Tier-1 = [057, 279, 074]; pool_mol_001 -> rank 22; pool_mol_022 -> rank 31

**Why check B matters and a rank check does not.**
A MAD-convention error is a uniform scale factor. It changes F(Gamma) rank 1 from
13.9908 to 20.7428 -- exactly 1.4826x -- and leaves the rank ORDER completely
unchanged. Any check based on ranks is structurally blind to it. Only an
element-wise comparison against the published Z columns catches it.

Usage:
  cd /mnt/hpc/home/wanghongmei/lang2
  /mnt/hpc/home/wanghongmei/anaconda3/envs/md_env/bin/python Scoring/chk_frozen_scale_R5.py
"""
import hashlib
import pathlib
import sys

import numpy as np
import pandas as pd

ROOT = pathlib.Path("/mnt/hpc/home/wanghongmei/lang2")
RES = ROOT / "results"
CAND_IN = RES / "step8d_input_combined_SOLV_V3.csv"
CAND_RANK = RES / "step8d_compass_ranking_SOLV_V3.csv"

AXES = ["mindist_mean", "hbond_occupancy", "rmsd_prot_mean",
        "dG_PBSA_mean", "sigma_PBSA", "rmsd_lig_mean"]

Z_COL = {"mindist_mean": "Z_d_min", "hbond_occupancy": "Z_H_occ",
         "rmsd_prot_mean": "Z_RMSD_prot", "dG_PBSA_mean": "Z_dG_PBSA",
         "sigma_PBSA": "Z_sigma_PBSA", "rmsd_lig_mean": "Z_RMSD_lig"}

# R3 of the Scope B registration -- invariant, never refitted
W_V3 = {"mindist_mean": -2.0000, "hbond_occupancy": +1.5078,
        "rmsd_prot_mean": -1.2227, "dG_PBSA_mean": -1.9317,
        "sigma_PBSA": -1.2738, "rmsd_lig_mean": -1.1590}

# R1 of the Scope B registration -- (median, 1.4826 * MAD) on the frozen 31
REGISTERED = {
    "mindist_mean":    (0.193500, 0.007710),
    "hbond_occupancy": (0.909500, 0.103634),
    "rmsd_prot_mean":  (2.286200, 0.268944),
    "dG_PBSA_mean":    (0.783300, 4.072257),
    "sigma_PBSA":      (2.334700, 0.995269),
    "rmsd_lig_mean":   (4.322000, 2.459485),
}

EPS = 1e-9
TOL_SCALE, TOL_Z, TOL_F = 1e-5, 1e-6, 0.05
EXPECT_TIER1 = ["pool_mol_057", "pool_mol_279", "pool_mol_074"]
EXPECT_RANKS = {"pool_mol_001": 22, "pool_mol_022": 31}

fail = []


def banner(s):
    print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)


def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def robust_z(x, med, scale):
    return (np.asarray(x, dtype=float) - med) / (scale + EPS)


banner("0 · inputs")
for p in (CAND_IN, CAND_RANK):
    if not p.is_file():
        sys.exit("MISSING: %s" % p)
    print("%-44s %8s B  sha256 %s" % (p.name, format(p.stat().st_size, ","), sha(p)))

cand = pd.read_csv(CAND_IN)
rank = pd.read_csv(CAND_RANK)
if len(cand) != 31:
    sys.exit("REFUSED: input table has %d rows, expected 31 (28 would mean a V2 file)" % len(cand))
if len(rank) != 31:
    sys.exit("REFUSED: ranking table has %d rows, expected 31" % len(rank))
m = cand.merge(rank, on="mol_id", suffixes=("", "_pub"))
if len(m) != 31:
    sys.exit("REFUSED: mol_id did not align, merged to %d/31" % len(m))
print("31 primary candidates, mol_id aligned 31/31")

banner("A · frozen scale recomputed vs registered R1")
FROZEN = {}
print("%-18s %12s %12s   %12s %12s   %s" %
      ("axis", "median", "1.4826*MAD", "registered", "registered", "max |diff|"))
worst_scale = 0.0
for a in AXES:
    v = cand[a].astype(float).values
    med = float(np.median(v))
    mad = float(np.median(np.abs(v - med)))
    if mad <= 0:
        fail.append("A: MAD = 0 on axis %s -- cannot standardise" % a)
    scale = 1.4826 * mad
    FROZEN[a] = (med, scale)
    rm, rs = REGISTERED[a]
    d = max(abs(med - rm), abs(scale - rs))
    worst_scale = max(worst_scale, d)
    print("%-18s %12.6f %12.6f   %12.6f %12.6f   %.2e" % (a, med, scale, rm, rs, d))
if worst_scale >= TOL_SCALE:
    fail.append("A: recomputed scale differs from registered R1 by %.2e (tol %.0e)"
                % (worst_scale, TOL_SCALE))
print("\nmax deviation %.2e   tolerance %.0e   -> %s"
      % (worst_scale, TOL_SCALE, "PASS" if worst_scale < TOL_SCALE else "FAIL"))

banner("B · element-wise Z_* vs published columns  (sensitive to the MAD convention)")
missing = [Z_COL[a] for a in AXES if Z_COL[a] not in m.columns]
if missing:
    sys.exit("REFUSED: ranking table lacks published Z columns: %s" % missing)
worst_z = 0.0
for a in AXES:
    med, scale = FROZEN[a]
    z_mine = robust_z(m[a].values, med, scale)
    z_pub = m[Z_COL[a]].astype(float).values
    d = float(np.max(np.abs(z_mine - z_pub)))
    worst_z = max(worst_z, d)
    print("   %-18s -> %-14s max deviation %.3e" % (a, Z_COL[a], d))
if worst_z >= TOL_Z:
    fail.append("B: Z columns differ by %.3e (tol %.0e). Most likely cause: the MAD "
                "convention. A deviation near 48%% means 1.4826 was dropped or applied twice."
                % (worst_z, TOL_Z))
print("\nmax deviation %.3e   tolerance %.0e   -> %s"
      % (worst_z, TOL_Z, "PASS" if worst_z < TOL_Z else "FAIL"))

banner("B' · reverse assertion -- bare MAD MUST fail check B")
a0 = AXES[0]
med0, scale0 = FROZEN[a0]
z_bad = robust_z(m[a0].values, med0, scale0 / 1.4826)
d_bad = float(np.max(np.abs(z_bad - m[Z_COL[a0]].astype(float).values)))
print("   bare-MAD deviation on %s: %.4f   (must be >> %.0e)" % (a0, d_bad, TOL_Z))
if d_bad <= 1e-3:
    fail.append("B': bare MAD also passes check B -- check B is NOT sensitive to the "
                "MAD convention and is worthless as written. Stop and investigate.")
print("   -> %s" % ("PASS (check B is genuinely sensitive)" if d_bad > 1e-3 else "FAIL"))

banner("C · COMPASS_v3 reproduced from the invariant weights")
F_mine = np.zeros(len(m))
for a in AXES:
    med, scale = FROZEN[a]
    F_mine = F_mine + W_V3[a] * robust_z(m[a].values, med, scale)
F_pub = m["COMPASS_v3"].astype(float).values
dF = float(np.max(np.abs(F_mine - F_pub)))
print("   max deviation %.4f   tolerance %.2f   (residual = 4-dp rounding of the stored weights)"
      % (dF, TOL_F))
if dF >= TOL_F:
    fail.append("C: COMPASS_v3 differs by %.4f (tol %.2f)" % (dF, TOL_F))
print("   -> %s" % ("PASS" if dF < TOL_F else "FAIL"))

banner("D · Tier-1 membership and two cross-check ranks")
order = m.assign(F=F_mine).sort_values("F", ascending=False).reset_index(drop=True)
tier1 = order.mol_id.head(3).tolist()
print("   Tier-1 recomputed: %s" % tier1)
if tier1 != EXPECT_TIER1:
    fail.append("D: Tier-1 is %s, expected %s" % (tier1, EXPECT_TIER1))
for name, want in EXPECT_RANKS.items():
    idx = order.index[order.mol_id == name]
    if len(idx) != 1:
        fail.append("D: %s not found exactly once in the ranking" % name)
        continue
    got = int(idx[0]) + 1
    print("   %-14s rank %2d   expected %2d   %s" % (name, got, want, "ok" if got == want else "MISMATCH"))
    if got != want:
        fail.append("D: %s is rank %d, expected %d" % (name, got, want))
print("   -> %s" % ("PASS" if not [f for f in fail if f.startswith('D')] else "FAIL"))

banner("RESULT")
if fail:
    print("R5 SELF-TEST FAILED -- %d check(s):\n" % len(fail))
    for f in fail:
        print("  * " + f)
    print("\nDo NOT project any reference row onto this scale until these are resolved.")
    sys.exit(2)
print("R5 SELF-TEST PASSED.")
print("  The frozen scale in R1 is reproduced from the real 31-candidate table.")
print("  The element-wise Z check is sensitive to the MAD convention (B' confirms it).")
print("  The weights reproduce the published COMPASS_v3 and the published Tier-1.")
print("\nThe scale is verified. Reference rows may be projected onto it under R2, and")
print("re-normalising on any enlarged table remains prohibited.")
print("\n(read-only: nothing was written)")
