# PRISM-Dyn V3 — data and code accompanying the manuscript

This repository holds the code and derived data required to check the claims of
the PRISM-Dyn V3 manuscript: a three-level audit of post-docking
molecular-dynamics rescoring against three coronavirus main proteases
(SARS-CoV-2 6W63, SARS-CoV-1 3V3M, MERS-CoV 4YLU).

Assembled by `stage_githubdyn.py` v1.7. Released under the MIT License
(`LICENSE`).

## Boundary

What is deposited is what is needed to check the manuscript's claims.
Generative-side components — decoder weights, surrogate models, training
scripts — are **not** here; the reasons are in `EXCLUSIONS.md`. Trajectory files
are excluded for size and are available from the corresponding author on
request.

Names that were searched for and not found are recorded in `NOT_FOUND.md`,
together with the patterns used, so that an absence can be distinguished from an
omission.

## What a clone contains

A clone contains **691 files totalling 59.38 MB** (decimal megabytes, 10⁶
bytes, summed over the files `git ls-files` enumerates).

**682** of them are recorded in `MANIFEST.sha256`, each with its SHA-256 digest
and the path it was taken from, so that a clone can be checked against the
archive it was built from rather than against a description of it.

The remaining **nine** are the repository's own documentation, configuration and
digest records, and are not listed in the manifest because they are not products
of the pipeline:

```
.gitattributes                  .gitignore
EXCLUSIONS.md                   LICENSE
MANIFEST.sha256                 NOT_FOUND.md
PROVENANCE.md                   README.md
02_md_inputs/IDENTICAL.sha256
```

`02_md_inputs/IDENTICAL.sha256` records every source path and content digest
behind the de-duplicated `.mdp` and `.json` inputs, of which only one
representative per distinct content is deposited.

## Directory layout

| directory | holds |
|---|---|
| `01_scoring/` | composite weights, rankings, bootstrap and charge-gate outputs |
| `02_md_inputs/` | simulation inputs for the three arms, plus the de-duplication record |
| `03_endpoint/` | MM-PBSA and MM-GBSA endpoint products |
| `04_structural/` | main-arm structural audit, occupancy, per-frame series |
| `05_cross_ortholog/` | cross-ortholog arm products |
| `06_controls/` | reference panel, registered criteria, docking controls |
| `07_audit/` | audit scripts run against the analysis products, and the audit of this deposit |
| `08_pool_and_prep/` | candidate pool and preparation records |
| `09_notebooks/` | pipeline notebooks and the helper modules they import |

`09_notebooks/core/` holds the preparation and docking notebooks together with
the Python modules those notebooks import directly, so that the notebooks run
from a fresh clone without a `sys.path` adjustment.

## The nine checks defined on this deposit

Eight run against the assembled tree and are reported by
`audit_githubdyn.py`. The ninth is by construction deferred to a clone.

The audit script is deposited at `07_audit/Scoring/audit_githubdyn.py`. It is
read-only, with no write path, and it shares no code with the assembler that
produced the tree, so that it cannot establish merely that the tree agrees with
the functions that built it: digests, manifest parsing, counts and notebook
parsing are all implemented independently within it. From the root of a clone:

```
python3 07_audit/Scoring/audit_githubdyn.py .
```

Run against a clone it consults `git` directly, which is what makes check H
below authoritative there.

**A · Manifest completeness.** `MANIFEST.sha256` and the tree must agree in both
directions, and every digest must match the file. What the assembler says it
copied and what the tree actually holds are two different things.

**B · Patches landed.** The six deposited notebooks must be the corrected
versions, compared against known SHA-256 values. A deposit carrying the
pre-correction copies would ship prose and code that disagree.

**C · Arm-name residue and arm switching.** The three cross-ortholog notebooks
must carry no `ARMS = ['md_solv_v3']` in their markdown; the three main-arm
notebooks' switch line must read `md_solv_v3`.

**D · Series counts in the tree.** The per-frame series directories must hold
99, 36, 33, 9 and 12 files, and the superseded `aud_cross_site_series/` must be
absent. This checks the file system.

**E · Per-replica counts.** `analysis.json` must number 99 for the main arm, 36
for the cross-ortholog arm and 12 for the reference arm.

**F · Exclusions traceable.** `EXCLUSIONS.md` must carry its registry of
directory-level exclusions and must name `aud_cross_site_series` within it, so
that a pointer from elsewhere does not lead to a blank.

**G · Declarations present.** This file must carry the per-frame series section,
the explanation of the directory name that never existed, the note on inherited
registration text, and the note on arm switching.

**H · Whether git will actually take them.** The deliberately deposited `.npz`
files must survive `.gitignore`. This checks git's verdict, not the file system,
and the two cannot substitute for one another: a tree can be complete, a
manifest can be complete, check D can pass, and git can still take none of them.
Where `git` cannot be consulted — on the compute cluster, where the deposit is
assembled — this check falls back to reading `.gitignore` as text, which the
audit reports as the weaker evidence.

**The ninth check belongs to the reader.** Clone the repository and verify every
digest in `MANIFEST.sha256` against the clone. `.gitattributes` sets `* -text`,
disabling line-ending conversion, precisely so that this check holds on any
platform.

## A limit on what the archive establishes

Two caveats concern how an instrument is identified, and both would otherwise
surface only on attempting a reproduction.

**A naming discrepancy.** One verification instrument appears in the execution
log under a filename differing from that of the archived copy, while the
archived copy's hash is exactly the value the log records. One artefact, two
names. The file has not been renamed to match the log, because the log records
what was run and the hash is what identifies it. A reproduction should resolve
an instrument by hash rather than by name.

**An archived hash identifies the archived copy, not necessarily the copy that
ran.** For two instruments of the cross-ortholog arm, the hash recorded in the
code archive does not match the hash of the copy present on the compute system.
For one of the two — the arm's construction driver — the divergence is
substantive rather than incidental: the archived copy names the registered
three-compound set, while the copy on the system names the five compounds the
arm actually executed. That difference reconciles the arm's registered and
executed trajectory counts at the level of source rather than of documents, and
it locates the substitution precisely: it was made by editing one constant, not
by writing a different program.

**The scope of that comparison is two instruments.** We report it for those two
because we performed it for those two. **We do not claim to have performed it
for every instrument in the archive.** Where an archived copy and a running copy
can be compared they should be, and where they differ the difference is itself a
finding.

## The cross-ortholog six-axis table: why it is not here, and how to rebuild it

`cross_six_axis_SOLV_V3.csv` is derived, and so is not deposited separately. It
can be rebuilt value-for-value from three files under `05_cross_ortholog/`:

- `phase3_cross_dG_aggregated_SOLV_V3.csv` — ΔG and σ
- `phase4_cross_structural_SOLV_V3.csv` — distance, occupancy, RMSD
- `aud_cross_site_occupancy.csv` — per-replica site occupancy

Seven of its columns (`dG_PBSA_mean`, `sigma_PBSA`, `n_rep`, `rmsd_prot_mean`,
`rmsd_lig_mean`, `mindist_mean`, `hbond_occupancy`) are the first two files
joined on (`mol_id`, `target`), **value-for-value identical** (measured). The
remaining three are recomputed from the third file under these definitions, also
value-for-value identical:

- `occ_min` — the **minimum** per-replica `site_occupancy`
- `d_dyad_max` — the **maximum** of the per-replica `d_dyad_mean`, that is, the
  largest of the per-replica means. It is **not** the true per-frame maximum.
- `escaped` — the **number of replicas** with `site_occupancy` < 0.50

**These three definitions must be cited with the data.** Without them the claim
of rebuildability cannot be checked, and `d_dyad_max` in particular invites
being read as a true per-frame maximum.

## Per-frame d_dyad series (`*_series/*.npz`): why they are released

Site occupancy is defined as the fraction of frames in which the minimum image
distance from any ligand atom to any catalytic-dyad atom is at most **5.0 Å**,
and **that 5.0 Å threshold was chosen after the fact** — no pre-registration
entry defines escape numerically. This is recorded in each `.meta.json` under
`threshold_provenance`.

Since the threshold was chosen after the fact, publishing only the occupancies
computed under it would leave a reader able to do nothing but accept the number
we chose. The per-frame series are therefore released: **they let a reader
recompute occupancy at any other threshold, and so measure how far the
conclusions depend on that choice.**

Each `.npz` carries `d_dyad` (a per-frame scalar, in nm), `start`, `stride` and
`dyad_resids`. **They are not trajectories** and contain no coordinates.

| directory | files | contents |
|---|---|---|
| `04_structural/results/aud_b81_scope9_series/` | 9 | per-frame d_dyad for the nine §B8.1 gate replicas |
| `04_structural/results/aud_main_site_occupancy_his163_series/` | 33 | per-frame d_dyad and d_extra for the His163 subset |
| `04_structural/results/aud_main_site_occupancy_series/` | 99 | main-arm per-frame d_dyad |
| `05_cross_ortholog/results/aud_cross_site_occupancy_series/` | 36 | cross-ortholog per-frame d_dyad |
| `06_controls/results/aud_ref_site_occupancy_series/` | 12 | reference-panel per-frame d_dyad |
| **total** | **189** | |

**A further `aud_cross_site_series/` (30 files, 2026-10-19) is not here.** It is
the generation of cross-ortholog series that preceded the arm reaching 36
trajectories, and its CSV and sidecar no longer exist. The grounds for excluding
it are the file count and the mtime; the reason is in `EXCLUSIONS.md`.

## Three things declared with the deposit

### 1 · Inherited registration text

The `registration` field of `aud_cross_site_occupancy.meta.json` cites a
pre-registration clause concerning the **His163 column** (prereg §B item 7),
while that product's own `extra_resid` is empty and it has no His163 column.
**That is inherited text from a script default, and must not be read as
establishing that the cross-ortholog audit was pre-registered under that
clause.** The clause actually refers to `aud_main_site_occupancy_his163.*`,
which does carry that column and does fall within its scope. The sidecar is not
rewritten in place: it is a dated product, and re-running 36 trajectories to
correct one text field is not a trade worth making.

### 2 · A directory name that never existed

The `ALLOWED_ARMS` list in the three cross-ortholog notebooks contains
**`md_cross_solv_v3`**, and that directory **has never existed**. The MD stage
went through several attempts and no directory was ever deleted, but none was
ever given this name. It was written from an assumption and propagated through
two generators.

**It is deliberately not removed from the code.** Deleting an element that no
notebook has ever used as `ARM` would not change the outcome of
`set(ARMS) <= ALLOWED_ARMS`; leaving it lets a reader see how an invented name
spread through a pipeline. The three real arm names are `md_solv_v3` (main),
`md_solv_v3_ref` (reference) and `cross_solv_v3` (cross-ortholog).

### 3 · Arm switching in the main-arm notebooks

`phase2_pbc_correct_V3.ipynb`, `phase3_mmpbsa_solv_V3.ipynb` and
`phase4_analysis_solv_V3.ipynb` each serve more than one arm. Switching changes
one line in cell 1, `ARM = EXPECTED_ARM = "…"`, which carries its own explanatory
comment. **All three are delivered set to the main arm `md_solv_v3`** (corrected
2026-10-31: two of them had been left at `md_solv_v3_ref` while their prose and
retained outputs were the main arm's). To take the reference arm, change that
line to `md_solv_v3_ref`.

## Two hazards declared with the deposit

1. `phase4_analysis_solv_V3.ipynb` defines `STRIDE = 10` twice, in cell 2 and
   cell 5. The values agree, so there is no error at present; a constant defined
   twice becomes a hazard the moment one of the two is edited alone.
2. The C3 section of `aud_frame_census.py` is withdrawn: its C3 and VERDICT must
   not be cited, while C1 and C2 stand. For a terminal-marker census use
   `07_audit/Scoring/aud_terminal_markers.py`, which recognises both legal
   terminal states, `[MD2 DONE]` and `[RESUME DONE]`.

## Citation

Please cite the PRISM-Dyn V3 manuscript. A full citation will be added here on
publication.
