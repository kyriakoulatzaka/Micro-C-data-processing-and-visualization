# Balancing-weight inspection and chromosome subsetting

## Previously performed analysis on real Micro-C data

The author inspected /bins/weight in merged 400-bp Micro-C matrices for **both nc1–8 and nc14** and compared raw and balanced Cooler outputs. The following values were reported for **nc1–8 only**:

| Metric | nc1–8 |
| --- | ---: |
| Total bins | 360,163 |
| NaN weights | 61,225 (17.0%) |
| Zero weights | 0 |
| Min finite weight | 0.000544426474887378 |
| Max finite weight | 9.79914356645668 |
| Mean finite weight | 0.025708628100725318 |

A 2,000 × 2,000 example block returned 4,000,000 finite raw entries versus 3,606,201 finite balanced entries, and `balance=True` did not raise an error. This demonstrates **readability** of existing weights, **not** that balancing convergence or reliability was established.

Some unplaced/random contigs were entirely NaN, but ranking by *fraction NaN* cannot establish that most NaNs came from those contigs. Compare **absolute NaN counts** across all chromosomes, particularly 2L, 2R, 3L, 3R, X, and 4.

Bins with missing balancing weights produce NaN contacts in balanced views. A raw matrix sum cannot be directly compared with a balanced matrix sum, which is on a different scale. Possible explanations for missing weights include masked, low-coverage, or otherwise excluded bins; determining the exact cause requires additional QC.

The author also attempted to restrict the contact matrix to primary chromosomes. **The original notebook code filtered pixels without remapping bin IDs**; the produced subset (if any) is not considered validated. This script repairs that issue while streaming input pixels. It does not alter or invalidate the original merged file.

## Inspect both developmental stages

```bash
python scripts/qc_and_subset.py inspect /path/to/nc1-8_merged_400bp.cool --window-bins 500
python scripts/qc_and_subset.py inspect /path/to/nc14_merged_400bp.cool --window-bins 500
```

The optional window uses a small square block to reduce memory load.

## Create a safe canonical-chromosome subset

```bash
python scripts/qc_and_subset.py subset \
  /path/to/nc1-8_merged_400bp.cool \
  /path/to/nc1-8_merged_400bp_mainchroms_corrected.cool
```

The subset retains 2L, 2R, 3L, 3R, X, and 4 (also accepting a `chr` prefix), drops other bins and their contacts, and **remaps both pixel bin IDs** to the new axis. It omits original weights; optionally calculate new weights after checking the output:

```bash
cooler balance /path/to/nc1-8_merged_400bp_mainchroms_corrected.cool
```

Use a new filename, keeping the original merged data unchanged. A separate rebalancing result must still be inspected for NaNs and convergence.

## Validation status

The author already performed the real-data weight checks reported above. The reorganized scripts in this repository have been unit-tested using small **synthetic** data, but have not yet been executed end-to-end on the author's real merged matrices. This is a software-validation limitation, not a claim that the author's work used synthetic data.
