# What was deliberately excluded, and why

Originally generated 2026-11-06 20:12:14 +0800; translated and extended
2026-11-23 (see `PROVENANCE.md`).

## 1 · Excluded by name (files that would mislead)

- **`step8c_compass_v3_ranking_SOLV.csv`** — named v3 but containing 28 rows
  (the V2 cohort). A reader would take it from its filename for this
  manuscript's cohort. To deposit it at all it would have to be renamed
  `…_main28_…` with a README in the same directory stating that it is a
  superseded cohort.
- **`permutation_test_report_SOLV.txt`** — contains a sentence marked
  "PAPER-READY SENTENCE" citing n = 28, r = 0.116. This manuscript's values are
  ρ = +0.078, p = 0.676, n = 31. The report would hand a reader a sentence this
  manuscript has excluded.
- **`phase3_cross_dG_newTier1_SOLV.csv`** — **a V2 cross-arm product, not this
  manuscript's.** Six rows, 2026-06-23, molecules `005` / `027` / `028`
  (confirmed by content in `Cross_arm_decision_2026-10-29.md` §1.2). This
  manuscript's cross-ortholog arm is six molecules — `004` / `012` / `057` /
  `074` / `279` / `001` — in 12 systems and 36 trajectories. Both filenames
  begin `phase3_cross_dG_`, so placed together under `05_cross_ortholog/` a
  reader would take it for part of this work. No claim here rests on it.

## 2 · Third layer: generative-side models, not deposited here

Reason: this manuscript makes no claim about the efficacy of the generative
stage. Depositing the generative models alongside it would read the
active-learning material into this work, and those components are the evidential
basis of the next paper; filing them under this manuscript's DOI would confuse
which paper's claims they support.

- `model_best_V2.pth` — decoder weights — generative side, evidential basis of the next paper (PRISM-Gen 2.0)
- `activity_predictor.pkl` — pIC50 surrogate — as above
- `activity_predictor_BBS_v*.pkl` — BBS surrogate versions — as above, and the comparison between versions belongs to that paper
- `train_corrected.py` — training script — as above

**One exception:** `fragments_V2.csv` **is** deposited. §2.1.1 identifies it by
content and hash (10,173 fragments / 319,428 bytes / MD5 `bc8816cc`), and an
identification made by hash whose file is not deposited is not an
identification a reader can use.

## 2a · Single-file exceptions that override a class exclusion

A class exclusion is a blunt instrument. Each file below falls inside an
excluded class and was deposited anyway, with its reason. **An unrecorded
exception is the same thing as a silent gap.**

- **`04_structural/results/aud_*_series/*.npz` · `05_cross_ortholog/results/aud_cross_site_occupancy_series/*.npz` · `06_controls/results/aud_ref_site_occupancy_series/*.npz`**
  — **189 per-frame ligand–dyad distance arrays, about 1.0 MB in total.**
  Section 4 excludes `.npz` by extension; these 189 are deposited regardless.
  The reason is §6.2: the coordinate files of the three arms total about 470 GB
  and are not released, and **these 189 arrays are deposited in their place**, so
  that a reader can recompute occupancy at thresholds other than 5.0 Å. A claim
  that rests on arrays, whose arrays are not deposited, is not a claim a reader
  can check. The distribution is 99 main arm, 36 cross-ortholog, 12 reference,
  plus two recomputations over main-arm trajectories (33 for the His163 subset,
  9 for b81_scope9), totalling 189. **The `.npz` entry in section 4 should be
  read together with this one.**

- **`release/PRISM_Dyn/scripts/step7_pbsa_fix_v4.in`** — its parameters agree
  item by item with those reported in the manuscript (startframe 5001 /
  endframe 10000 / interval 25 / radiopt 1), and agree independently with
  `same_frames_as_pbsa: (5001,10000,25)` and `frames_per_replica_mmpbsa: 200`
  recorded in `results/step7_gbsa_solv_v3/gbsa_v3_build.json`. **This is "an
  input whose parameters agree with the run as recorded", not something
  established to be "the file that run read"** — the latter would require
  locating the v3 arm's step7 driver, and `README.md` states it in those terms.

## 3 · Excluded as whole directories

Several generations of superseded product trees exist under `lang2`. Excluding
them file by file does not hold them back, so they are excluded by **path
segment**. Among them, the `step7` and `step7_gbsa_solv` generations belong to
the **28-candidate cohort**, which May directed on 2026-10-27 be excluded from
this manuscript; released with the deposit, they would hand a reader endpoint
energies incompatible with this manuscript's 31 candidates.

- **the `release` tree is another release's assembly, with a different boundary from this deposit** — 1 file, e.g. `release/FRATTVAE_V2/fragments_V2.csv`

- **`tablefig/`** — the output directory for figures and tables, regenerated on
  every run of `09_notebooks/Scoring/make_tablefig_V3.ipynb`, and not deposited.
  Its contents (each table's `.md` and `.csv`, `CAPTIONS.md`, each figure's
  `.png` and `.pdf`, and the `MANIFEST.txt` recording library versions at run
  time) are all recomputable from the deposited data and that notebook.
  **The footnotes to Table 3 and Table S9 in the manuscript are worded
  accordingly — "generated by that notebook", not "deposited".**

## 3a · Directory-level exclusions on record (no matches in this assembly)

The directories below are registered as excluded by path segment. **No pattern
in this assembly reached them**, so they do not appear in the table above — but
the reason for excluding them must remain visible with the deposit, or the claim
that they were excluded cannot be checked.

- **`15_全文修改前备份/`** (*backup taken before the full-text revision*) — kept in its original Chinese because it is a real path on the compute system; translating it would point at a directory that does not exist. A path is not prose.
- **`_abandoned_20261015/`** — the abandoned first attempt at the cross-ortholog arm (the deviation is recorded in §3.6.1)
- **`_archive_prebugfix_20260621/`** — products from before a defect was fixed
- **`_archive_release_old/`** — a superseded release archive
- **`_s1_authoritative_inputs/`** — ⚠ **its name says authoritative, and its contents do not belong to this manuscript.** 70 files were read on 2026-10-30: it is not an MD input set (my earlier characterisation was wrong and is corrected here) but the record of the **ligand parameterisation (S1) stage** — per-molecule `ligand_authoritative.mol` with `parameterization_manifest.json` and `failure_manifest.json`, plus `_s1_bijection_audit.csv`, `_s1_provenance_frozen.json` and `_S1_FAILED.json`. It holds 22 `pool_mol_*` directories (dated 19–20 July), of which six — 002 / 018 / 023 / 025 / 026 / 028 — are **not** among this manuscript's 31 molecules, while fifteen of this manuscript's — 003 / 017 / 022 / 029 / 030 / 031 / 034 / 057 / 061 / 073 / 074 / 120 / 146 / 243 / 279 — are **not** among its; the intersection is only 16. It is therefore a superseded cohort. **Precisely because its name misleads, this exclusion must leave a record**: released with the deposit, a reader would take it from its directory name for this manuscript's authoritative inputs.
- **`aud_cross_site_series/`** — **the cross-arm per-frame series of the 30-trajectory generation (2026-10-19 18:28).** This manuscript's cross-ortholog arm is 36 trajectories (AMD-2, 2026-10-31), whose series are in `aud_cross_site_occupancy_series/` (36 files). The two names differ by one fragment, and deposited side by side a reader could not tell which corresponds to §3.6 of the manuscript.
- **`step7/`** — step7 products of 2026-06-01, belonging to the 28-candidate cohort (May directed their exclusion from this manuscript on 2026-10-27)
- **`step7_gbsa_solv/`** — step7 products of 2026-06-21, the same cohort
- **`test_single/`** — a single-system trial run
- **`tmp/`** — temporary products

## 3b · De-duplicated by content (one representative plus a hash proof)

- `02_md_inputs` / `md_solv_v3/*/r*/03_min/em.mdp`: 627 files holding only **9 distinct contents**; one representative of each is deposited. Every path and its sha256 is in `02_md_inputs/IDENTICAL.sha256`.
- `02_md_inputs` / `cross_solv_v3/*/r*/03_min/em*.mdp`: 304 files holding only **9 distinct contents**; representatives only, same proof file.
- `02_md_inputs` / `md_solv_v3/*/r*/06_md/ndx_groups.json`: 145 files holding only **1 distinct content**; representative only, same proof file.

## 4 · Excluded by extension

Trajectories and binaries: `.cpt`, `.edr`, `.gz`, `.h5`, `.log`, `.npz`,
`.pkl`, `.pth`, `.tar`, `.tng`, `.trr`, `.xtc`, `.zip`.

**`.npz` carries a recorded exception — see section 2a.** 189 per-frame distance
arrays are deposited despite this rule.
