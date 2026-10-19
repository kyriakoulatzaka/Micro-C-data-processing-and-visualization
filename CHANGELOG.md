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
