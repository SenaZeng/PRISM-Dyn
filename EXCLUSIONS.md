# Deliberate exclusions and their reasons

Generated 2026-11-06 20:12:14 +0800

## 1 · Hard exclusions by name (hazardous files)

- **`step8c_compass_v3_ranking_SOLV.csv`** — named v3 but containing 28 rows (the V2 cohort). A reader would take it from the filename to be this manuscript's cohort. To deposit it, it would have to be renamed to …_main28_… with a README in the same directory stating that it is a superseded cohort.
- **`permutation_test_report_SOLV.txt`** — contains a sentence marked "PAPER-READY SENTENCE" quoting n = 28 and r = 0.116. This manuscript's values are ρ = +0.078, p = 0.676, n = 31; the report would actively hand a reader a sentence this manuscript has excluded.
- **`phase3_cross_dG_newTier1_SOLV.csv`** — **a V2 cross-arm product, not part of this manuscript.** Six rows, 2026-06-23, molecules `005` / `027` / `028` (confirmed by content in `Cross_arm_decision_2026-10-29.md` §1.2). This manuscript's cross-arm comprises six molecules — `004` / `012` / `057` / `074` / `279` / `001` — in 12 systems and 36 trajectories. Both begin with `phase3_cross_dG_`, and placed side by side under `05_cross_ortholog/` the older one would be taken for part of this work. No claim in this manuscript rests on it.

## 2 · The third layer: generative-side models, not in this deposit

Reason: this manuscript makes no claim about the efficacy of the generative stage. Depositing the generative models alongside it would read the active-learning material into this manuscript, and these components are the evidential basis of the next paper; placing them under this manuscript's DOI would confuse which paper's claims they support.

- `model_best_V2.pth` — decoder weights — generative side, evidential basis of the next paper (PRISM-Gen 2.0)
- `activity_predictor.pkl` — pIC50 surrogate — as above
- `activity_predictor_BBS_v*.pkl` — the BBS surrogate versions — as above, and the comparison between versions belongs to that paper
- `train_corrected.py` — training script — as above

**One exception:** `fragments_V2.csv` is deposited. §2.1.1 identifies it by content and hash (10,173 fragments / 319,428 bytes / MD5 `bc8816cc`), and an identification made by hash is of no use to a reader if the file it identifies is not deposited.

## 2b · Single-file exceptions that override a whole-class exclusion

Excluding a whole class is a blunt instrument; each of the files below falls inside an excluded tree and is nonetheless deposited, with its reason attached. **An unrecorded exception is of the same kind as a silent gap.**

- **`release/PRISM_Dyn/scripts/step7_pbsa_fix_v4.in`** — its parameters agree item by item with those stated in the manuscript (startframe 5001 / endframe 10000 / interval 25 / radiopt 1), and agree independently with `same_frames_as_pbsa: (5001,10000,25)` and `frames_per_replica_mmpbsa: 200` as recorded in `results/step7_gbsa_solv_v3/gbsa_v3_build.json`. **This establishes it as an input whose parameters agree with the run as recorded, not as the file that run actually read** — the latter would require locating the step7 driver for the v3 arm, and the README must be worded accordingly.

## 3 · Whole-class exclusions by directory

Several generations of superseded product trees exist under lang2. Excluding them file by file does not hold them back, so they are excluded as classes, by **path segment**. Among them, the `step7` and `step7_gbsa_solv` generations belong to the **28-candidate cohort**, which May directed on 2026-10-27 to be excluded from this manuscript; released with the deposit, they would give a reader endpoint energies incompatible with this manuscript's 31 candidates.

- **The release tree is the assembly for a different release, with a different scope from this deposit** — 1 file, for example `release/FRATTVAE_V2/fragments_V2.csv`

### 3a · Register of whole-class exclusions (those with no match in this assembly)

The directories below are registered as excluded by path segment. **No pattern in this assembly pointed at them**, so the table above does not list them — but the reason for excluding them must remain visible with the deposit, or the claim that they were excluded is untestable.

- **`15_全文修改前备份/`** — a backup taken before the full-text revision
- **`_abandoned_20261015/`** — the abandoned first attempt at the cross-ortholog arm (the deviation is recorded in §3.6.1)
- **`_archive_prebugfix_20260621/`** — products from before the bug fix
- **`_archive_release_old/`** — a superseded release archive
- **`_s1_authoritative_inputs/`** — ⚠ **its name says authoritative, and its contents do not belong to this manuscript.** 70 files were read directly on 2026-10-30: it is not a set of MD inputs (my earlier characterisation was wrong and is corrected here) but the record of the **ligand parameterisation (S1) stage** — per-molecule `ligand_authoritative.mol` with `parameterization_manifest.json` and `failure_manifest.json`, together with `_s1_bijection_audit.csv`, `_s1_provenance_frozen.json` and `_S1_FAILED.json`. It holds 22 `pool_mol_*` directories (dated 19–20 July), of which six — 002 / 018 / 023 / 025 / 026 / 028 — are **not** among this manuscript's 31 molecules, while fifteen of this manuscript's — 003 / 017 / 022 / 029 / 030 / 031 / 034 / 057 / 061 / 073 / 074 / 120 / 146 / 243 / 279 — are **not** among its; the intersection is only 16. It therefore belongs to a superseded cohort. **It is precisely because the name is misleading that this exclusion must be recorded**: released with the deposit, a reader would take it from the directory name to be this manuscript's authoritative inputs.
- **`aud_cross_site_series/`** — **the cross-arm per-frame series from the 30-trajectory generation (2026-10-19 18:28).** This manuscript's cross-arm has 36 trajectories (AMD-2, 2026-10-31), whose series are in `aud_cross_site_occupancy_series/` (36 files). The two differ by one segment of the directory name, and deposited side by side a reader would have no way to tell which corresponds to §3.6 of the manuscript.
- **`step7/`** — step7 products of 2026-06-01, belonging to the 28-candidate cohort (May directed on 2026-10-27 that it be excluded from this manuscript)
- **`step7_gbsa_solv/`** — step7 products of 2026-06-21, belonging to the same 28-candidate cohort
- **`test_single/`** — a single-system trial run
- **`tmp/`** — temporary products

## 3b · Deduplication by content (a representative plus hash proof)

- `02_md_inputs` / `md_solv_v3/*/r*/03_min/em.mdp`: 627 files holding only **9 distinct contents**; a representative of each is kept, with every path and its sha256 in `02_md_inputs/IDENTICAL.sha256`
- `02_md_inputs` / `cross_solv_v3/*/r*/03_min/em*.mdp`: 304 files holding only **9 distinct contents**; a representative of each is kept, with every path and its sha256 in `02_md_inputs/IDENTICAL.sha256`
- `02_md_inputs` / `md_solv_v3/*/r*/06_md/ndx_groups.json`: 145 files holding only **1 content**; a representative is kept, with every path and its sha256 in `02_md_inputs/IDENTICAL.sha256`

## 4 · Exclusions by extension

Trajectories and binaries: .cpt, .edr, .gz, .h5, .log, .npz, .pkl, .pth, .tar, .tng, .trr, .xtc, .zip

Note that `.npz` appears in that list as a default rule, and the per-frame `d_dyad` series under `07_audit/` are nonetheless deposited under an explicit exception; the reasons are given in the README. They are not trajectories and contain no coordinates.
