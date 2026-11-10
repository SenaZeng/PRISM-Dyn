# PRISM-Dyn V3 — Data and code accompanying the manuscript

> This directory was assembled by `stage_githubdyn.py` v1.7. **Still to be added**: the citation, and a map from the deposit to the sections of the manuscript.

## Scope

This deposit contains what is required to check the claims made in the manuscript. The generative-side models (decoder weights, surrogate models, training scripts) are **not here**; the reasons are given in `EXCLUSIONS.md`. Trajectory files are held separately on account of their size.

## Directories

- `01_scoring/`
- `02_md_inputs/`
- `03_endpoint/`
- `04_structural/`
- `05_cross_ortholog/`
- `06_controls/`
- `07_audit/`
- `08_pool_and_prep/`
- `09_notebooks/`

## The six-axis table for the cross-ortholog arm: why it is not here, and how to rebuild it

`cross_six_axis_SOLV_V3.csv` **is a derived product and is therefore not deposited separately.** It can be rebuilt value-for-value from three files under `05_cross_ortholog/` in this deposit:

- `phase3_cross_dG_aggregated_SOLV_V3.csv` (ΔG and σ)
- `phase4_cross_structural_SOLV_V3.csv` (distance · occupancy · RMSD)
- `aud_cross_site_occupancy.csv` (per-replica site occupancy)

Seven of its columns (`dG_PBSA_mean` · `sigma_PBSA` · `n_rep` · `rmsd_prot_mean` · `rmsd_lig_mean` · `mindist_mean` · `hbond_occupancy`) are the result of joining the first two files on (`mol_id`, `target`), and agree **value for value** (verified by direct comparison). The remaining three columns are recomputed from the third file under the definitions below, and likewise agree **value for value**:

- `occ_min` — the **minimum** of the per-replica `site_occupancy`
- `d_dyad_max` — the **maximum** of the per-replica `d_dyad_mean`, that is, the largest of the per-replica means. It is **not** the true per-frame maximum (the column heading in §3.6.2 follows this reading)
- `escaped` — the **number of replicas** with `site_occupancy` < 0.50

**These three definitions must be cited together with the data.** Without them, the claim that the table "can be rebuilt" is an assertion no reader can test; `d_dyad_max` in particular is easily misread as a true per-frame maximum.

## Per-frame d_dyad series (`*_series/*.npz`): why they are released with the deposit

The criterion for `site_occupancy` is the fraction of frames in which the minimum image distance from any ligand atom to any atom of the catalytic dyad is ≤ **5.0 Å**, and **that 5.0 Å threshold was chosen after the fact** — nothing in the record pre-registers a numerical definition of "escape". This is recorded in the `threshold_provenance` field of each `.meta.json`.

Since the threshold was chosen after the fact, giving only the occupancies computed under it would leave a reader no option but to adopt the one value we chose. The per-frame series are therefore released with the deposit: **they allow a reader to recompute the measure at any other threshold, and so to test how far the conclusions depend on that choice.**

Each `.npz` carries `d_dyad` (a per-frame scalar, in nm), `start`, `stride` and `dyad_resids`. **They are not trajectories** and contain no coordinates; the trajectories are held separately.

| Directory | Count | Contents |
|---|---|---|
| `aud_b81_scope9_series/` | 9 | per-frame d_dyad for the 9 replicas under the §B8.1 gate |
| `aud_cross_site_occupancy_series/` | 36 | per-frame d_dyad for the cross-ortholog arm (the 36 replicas after AMD-2) |
| `aud_main_site_occupancy_his163_series/` | 33 | per-frame d_dyad and d_extra for the His163 subset |
| `aud_main_site_occupancy_series/` | 99 | per-frame d_dyad for the main arm |
| `aud_ref_site_occupancy_series/` | 12 | per-frame d_dyad for the reference arm |

**A further directory, `aud_cross_site_series/` (30 files, 2026-10-19), is not here**: it holds the series from the generation preceding the cross-ortholog arm's 36 replicas, and its CSV and sidecar no longer exist. The grounds for excluding it are the file count and the mtime; the reasoning is in `EXCLUSIONS.md`.

## Three matters that must be stated alongside the deposit

### 1 · One piece of inherited text that requires a statement

The `registration` field of `aud_cross_site_occupancy.meta.json` cites a pre-registration clause concerning the **His163 column** (prereg §B, item 7), while the product itself has an empty `extra_resid` and no His163 column. **That is inherited text arising from a script default, and it must not be taken to show that the cross-ortholog arm's audit was pre-registered under that clause.** The clause in fact refers to `aud_main_site_occupancy_his163.*`, which does have that column and does fall within its scope. The sidecar is not rewritten in place here: it is a dated product, and re-running 36 trajectories to correct one text field is not a worthwhile exchange.

### 2 · A directory name that never existed

The `ALLOWED_ARMS` set in six of the notebooks contains an entry **`md_cross_solv_v3`**, and that directory **has never existed**: the MD stage went through many attempts and every directory was kept rather than deleted, but none was ever given that name. It was written from an assumption, and propagated through both `make_cross_notebooks.py` and `core/patch_phase_arm_derived.py`.

**This entry is deliberately not removed from the code.** Deleting an element that no notebook has ever used as `ARM` would not change the outcome of `set(ARMS) <= ALLOWED_ARMS`; leaving it in place lets a reader see how an invented name propagated through a pipeline. There are three real arm names: the main arm `md_solv_v3`, the reference arm `md_solv_v3_ref`, and the cross-ortholog arm `cross_solv_v3`.

### 3 · Arm switching in the main-arm notebooks

`phase2_pbc_correct_V3.ipynb` · `phase3_mmpbsa_solv_V3.ipynb` · `phase4_analysis_solv_V3.ipynb` each serve more than one arm; switching means changing a single line in cell 1, `ARM = EXPECTED_ARM = "…"`, and that line carries its own explanatory comment. **All three are delivered set to the main arm `md_solv_v3`** (corrected 2026-10-31: two of them had previously been left at `md_solv_v3_ref`, while their prose and retained outputs were those of the main arm). To take the reference arm, change that one line to `md_solv_v3_ref`.

## Two hazards that must be stated alongside the deposit

1. `phase4_analysis_solv_V3.ipynb` defines `STRIDE = 10` twice, once in cell 2 and once in cell 5. The two values agree, so nothing is wrong at present, but a constant defined twice becomes a hazard as soon as one of the definitions is changed on its own.
2. C3 of `aud_frame_census.py` is withdrawn: its C3 and VERDICT must not be cited, while C1 and C2 still stand. For a census of terminal markers use `07_audit/aud_terminal_markers.py`, which recognises both legitimate terminal states, `[MD2 DONE]` and `[RESUME DONE]`.
