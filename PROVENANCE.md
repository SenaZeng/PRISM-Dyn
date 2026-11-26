# How this deposit was produced

**status: COMPLETE**

This deposit was assembled once by a script and amended once by hand. Both are
recorded here, because a provenance file that records only the clean half is
not a provenance file.

## Assembly

| | |
|---|---|
| assembled | 2026-11-06 20:12:14 +0800 |
| host | gpu01 |
| lang2 root | `/mnt/hpc/home/wanghongmei/lang2` |
| assembler | `stage_githubdyn.py` v1.7, sha256 `7e8747e21187e50f12bde932289961eb16f0f1337b32f19bc9ab9958b6291805` |
| files written | 643 (12 newly copied, 631 skipped as identical in content) |

**How the count can be checked.** The assembler's read-only report prints, for
each of its patterns, the number of matches and the running total, so the file
count in the repository can be traced pattern by pattern back to its sources
under `lang2`.

## Amendment, 2026-11-23

The assembled tree was amended by hand, outside the assembler. The assembler was
deliberately not re-run, because re-running it would have overwritten the
documentation files and would not have included the modules added below.

| change | count |
|---|---|
| Python modules and submission scripts added under `09_notebooks/core/` | +26 |
| duplicate notebook removed: `02_md_inputs/core/Node05_md_prep.ipynb`, byte-identical to the copy retained at `09_notebooks/core/Node05_md_prep.ipynb` | −1 |
| `09_notebooks/Scoring/make_tablefig_V3.ipynb` brought into the staging tree | +1 |
| `LICENSE` brought into the staging tree | +1 |

The 26 modules were added because the deposited notebooks import six of them
directly, and a clone without them raises `ImportError` on the first cell. The
import closure was checked: among local modules those six import only
`master_molecule`, which is in the set.

`MANIFEST.sha256` was then regenerated with every digest recomputed from disk,
in the same format the assembler writes, and the documentation files were
translated to English and corrected.

| | after amendment |
|---|---|
| files in a clone | 690 |
| recorded in `MANIFEST.sha256` | 681 |
| not recorded (documentation, configuration, digest records) | 9 |
| total volume | 59.35 MB (decimal, 10⁶ bytes) |

## A known consequence

Because the amendment was made outside the assembler, `stage_githubdyn.py` v1.7
**can no longer reproduce this deposit**: re-running it would produce a smaller
tree, without the 26 modules, and would overwrite the documentation files with
their generated Chinese originals. Bringing the assembler back into step with
this deposit is outstanding work, not a property of the deposit as it stands.

Source paths and digests for every recorded file are in `MANIFEST.sha256`.
Deliberate exclusions are in `EXCLUSIONS.md`. Names searched for and not found
are in `NOT_FOUND.md`.
