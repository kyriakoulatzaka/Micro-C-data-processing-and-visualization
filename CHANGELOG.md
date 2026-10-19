## v4 — QC, balancing weights, and corrected chromosome subsetting

- Documented the author's real-data inspection of nc1–8 and nc14 balancing weights.
- Recorded nc1–8 measurements: 61,225 NaN weights in 360,163 bins (17.0%).
- Added `scripts/qc_and_subset.py` for per-chromosome weight summaries and safe, streaming subsetting with bin-ID remapping.
- Added `tests/test_qc_and_subset.py` and `docs/qc_weights_and_mainchroms.md`.
- Distinguished previous real-data analysis from tests of the reorganized code.

# Changelog

## v3 — Micro-C regional heatmaps

- Added `scripts/plot_heatmap.py` for locus-specific Micro-C matrix extraction and static heatmap export to PNG, PDF or SVG.
- Supports `.cool` paths and `.mcool::/resolutions/<bin-size>` URIs.
- Defaults to raw count contacts (`balance=False`) with log-scaled YlOrBr palette, white zeros and no axes.
- Added explicit `--balanced`, `--vmin`, `--vmax`, `--title`, `--cmap`, `--dpi` and `--max-bins` options.
- Added `examples/regions.tsv` containing the user-specified **ems** nc14, dm6 locus window.
- Added `matplotlib` dependency and six synthetic/offline plot tests.
- Updated README with reproduction steps and scientific interpretation caveats.

The packaged repository contains code and a figures README, **not** the user's original experimental Micro-C datasets or PNG figures.

## v2 — replicate merging and COOL→HIC conversion

- Includes a replicate-merging workflow, a data inspection script, Cooler-to-Juicer format conversion and unit tests.
