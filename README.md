# Micro-C processing, `.cool` → `.hic` conversion, and locus-specific heatmaps

A **Bash + Python workflow** for processing public *Drosophila melanogaster* embryonic Micro-C contact matrices. It combines four replicate matrices at a common genomic resolution into a single Cooler (`.cool`) file and converts that merged file into Juicer (`.hic`) format for downstream exploration. It also creates regional Micro-C contact heatmaps directly from `.cool` or `.mcool` matrices with Python, Cooler, and Matplotlib. This is a binned-matrix workflow, **not raw FASTQ alignment or Micro-C read processing**.

> **Project status:** The contact conversion functions and the heatmap helpers have passed offline unit tests on synthetic data. The replicate-merging stage and complete `.cool` → `.hic` workflow **have not been executed end-to-end against the user’s original real data in this repository**. Treat the scripts as a documented draft requiring validation before making performance or biological claims.

## Background

The input is a binned Micro-C contact matrix, not a raw paired-end sequencing dataset. The pipeline exports stored (unbalanced) contact counts and reformats them as Juicer's `shortWithScore` representation before generating a `.hic` file. The resulting `.hic` is intended for visualization and inspection with compatible tools; it does not add new biological evidence or run contact calling.

**Example dataset:** NCBI GEO [GSE265818](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE265818), *3D Chromatin Structures Foreshadow the Drosophila Fate Map*. The nc1–8 samples are GSM8229625–GSM8229628 (two biological × two technical replicates); the original multi-resolution matrices include 400-bp bins and were aligned against the dm6 genome assembly (BDGP Release 6). Cite the dataset and publication in any derivative work.

## Tools

- Linux (or WSL), Bash, GNU sort, gzip
- Python 3.10, Cooler 0.10.4, h5py, NumPy, Matplotlib
- Java (OpenJDK 17 recommended for testing compatibility)
- [Juicer Tools](https://github.com/aidenlab/juicer) (the original workflow used `juicer_tools_1.22.01.jar`)
- `dm6.chrom.sizes` matching the chromosome naming scheme in the `.hic`

## Setup

Install [Miniforge](https://github.com/conda-forge/miniforge) first, then create the environment:

```bash
mamba env create -f environment.yml
mamba activate microc2hic
java -version
```

Install a Java runtime if necessary. Download your Juicer Tools JAR from its official distribution source. Download `dm6.chrom.sizes` from UCSC:

```bash
wget -O dm6.chrom.sizes \
  https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes
```

## Stage 1 — Merge the four nc1–8 replicate contact matrices

The GEO supplementary files are already **multi-resolution `.mcool` files**. If they contain `::/resolutions/400`, there is no need to run `cooler zoomify` to create that resolution again. The previous notebook used filenames ending in `.cool` for the output of `zoomify`, but this output normally has the **multi-resolution HDF5 structure**; the extension alone does not change that.

Copy or download the four files into a single directory, for example `~/microc_data/GSE265818/`. Then execute:

```bash
mkdir -p ~/microc_data/merged/nc1-8
bash scripts/merge_microc_replicates.sh \
  ~/microc_data/merged/nc1-8/nc1-8_merged_400bp.cool \
  400 \
  "$HOME/microc_data/GSE265818/GSM8229625_nc1-8_1_1.base_resolution_100.mcool::/resolutions/400" \
  "$HOME/microc_data/GSE265818/GSM8229626_nc1-8_1_2.base_resolution_100.mcool::/resolutions/400" \
  "$HOME/microc_data/GSE265818/GSM8229627_nc1-8_2_1.base_resolution_100.mcool::/resolutions/400" \
  "$HOME/microc_data/GSE265818/GSM8229628_nc1-8_2_2.base_resolution_100.mcool::/resolutions/400"
```

The script checks that every input uses the same **400-bp genomic bin axis** before merging, then combines raw count matrices with `cooler merge` (sum, **not mean**). It refuses to overwrite an existing output. To use a new location, adjust the quoted input paths. **Do not describe merging as an assessment of replicate concordance:** read-depth/QC comparisons and concordance would be a separate analytical module.

### Optional: balance the merged matrix

```bash
cooler balance ~/microc_data/merged/nc1-8/nc1-8_merged_400bp.cool
```

This is an optional downstream normalization step; its parameters and success must be checked for the specific dataset. `cooler merge` does not preserve `/bins/weight` from the inputs. Verify chromosome coordinate conventions and assembly using dataset metadata; editing a `genome-assembly` text attribute **does not lift over genomic coordinates**.

## Stage 2 — Inspect your merged `.cool`

```bash
python scripts/inspect_cool.py /path/to/nc1-8_merged_400bp.cool --expected-resolution 400
```

This reports the bin size, number of pixels and whether the file already contains `bins/weight` values. **The conversion uses raw `count` values, not the existing balanced weights.** Normalization vectors recalculated by Juicer are not equivalent to copying Cooler weights.

## Stage 3 — Convert merged `.cool` to `.hic`

```bash
bash scripts/convert_cool_to_hic.sh \
  /path/to/nc1-8_merged_400bp.cool \
  /path/to/nc1-8_merged_400bp.cisOnly.hic \
  /path/to/dm6.chrom.sizes \
  /path/to/juicer_tools_1.22.01.jar \
  400
```

By default the workflow includes **only cis contacts** (same-chromosome interactions). Add `--include-trans` at the end if you need genome-wide interchromosomal contacts instead.

For larger datasets you may need to adjust disk space and memory:

```bash
JAVA_HEAP=20g SORT_MEMORY=2G TMP_BASE=/path/to/large/temporary/disk \
  bash scripts/convert_cool_to_hic.sh INPUT.cool OUTPUT.hic dm6.chrom.sizes juicer_tools.jar 400
```

The workflow streams Cooler pixels to avoid writing a redundant uncompressed dump, but sorting still requires potentially substantial temporary disk space. It deletes its temporary folder after completion.

## Stage 4 — Plot locus-specific Micro-C heatmaps

This optional stage plots a square contact matrix for one genomic interval directly from a **merged** `.cool` or original `.mcool` resolution group. **It does not need a `.hic` file.** The initial example from the project author's original script targeted the **ems** locus at nc14 (3R:13778677–14027796, dm6), which is distinct from the **nc1–8** merge example described above. To reproduce the ems plot you need the separate `nc14_merged_400bp.cool` file; it is *not* provided in this repository.

Run from the repository root:

```bash
python scripts/plot_heatmap.py \
  --cool /path/to/nc14_merged_400bp.cool \
  --region '3R:13778677-14027796' \
  --out figures/ems_nc14_merged_400bp.png
```

You can also plot a locus from the **nc1–8** matrix; supply a real region matching that input:

```bash
python scripts/plot_heatmap.py \
  --cool /path/to/nc1-8_merged_400bp.cool \
  --region '3R:13778677-14027796' \
  --out figures/ems_nc1-8_merged_400bp.png
```

This stage recreates the earlier paper-inspired plotting style (YlOrBr colormap, white for missing and zero counts, logarithmic colorbar, no axes), but fixes the plot scale automatically based on positive values rather than imposing `vmin=1e-6` on unbalanced integer counts. `--balanced` will attempt to use existing Cooler `/bins/weight`; by default **unbalanced/raw counts** are used, exactly as in the original `balance=False` script. Nothing in this step implies those counts have been normalized. For comparable heatmaps across stages or loci, explicitly choose the *same* color scale and normalization mode:

```bash
python scripts/plot_heatmap.py \
  --cool /path/to/nc14_merged_400bp.cool \
  --region '3R:13778677-14027796' \
  --out figures/ems_nc14_fixedscale.png \
  --vmin 1 --vmax 100 --title 'ems | nc14 | dm6 | 400 bp (raw counts)'
```

The `--balanced` flag is only valid when the input stores suitable balancing weights; if missing, run a suitable normalization step first. The tool accepts `3R` or `chr3R` naming where an unambiguous alias exists in the file. It limits the dense extracted matrix to 2500 bins per axis by default to reduce accidental memory exhaustion (override with `--max-bins` if necessary).

**Scientific caution:** Fixed raw count color scales do not compensate for differences in sequencing depth between experiments; for cross-stage quantitative interpretation, evaluate library-depth and normalization choices. The `ems` window and its coordinates came from the author's input script; this repository includes **no experimental heatmap image** because the source `.cool` file and generated PNG have not been supplied. The figure-rendering test uses only synthetic matrix values and does not validate the experimental contacts.

## Run the small tests

No Micro-C file is required for the unit tests:

```bash
python -m unittest discover -s tests -v
```

You can also manually test the **conversion stage only** with a small example:

```bash
python scripts/convert_contacts.py --chrom-sizes examples/toy.chrom.sizes --cis-only \
  < examples/toy.contacts.tsv
```

Expected example result:

```text
0\tchr3R\t1\t0\t16\tchr3R\t401\t1\t5
0\tchr3R\t401\t0\t16\tchr3R\t1201\t1\t3
```

The printed output uses tab separators.

## Important considerations

1. **Chromosome names**: Cooler files may use `3R` while UCSC dm6 uses `chr3R`. The Python converter maps names to the provided chromosome sizes file and avoids creating invalid `chrchr3R` names.
2. **Coordinates**: Cooler bin starts use zero-based coordinates; Juicer records here use one-based positions (`start + 1`). This must be evaluated in an end-to-end comparison of the resulting contact map.
3. **Counts vs balancing**: This pipeline exports raw integer counts. Do not mistake them for ICE/KR-balanced values. Existing `/bins/weight` is inspected but not applied.
4. **Storage and resolution**: Cooler stores only the upper triangle for `symmetric-upper` files. At 400 bp, genome-wide matrices can be very large. Ensure enough temporary disk space and Java heap memory.
5. **`.hic` validation**: A successful Java process alone is insufficient. Check that the `.hic` opens, exposes the requested resolution, and retains expected cis interaction counts at selected bins. Check whether the required normalization vectors were generated.
6. **Input semantics**: This is a conversion of pre-binned aggregate contacts, not a raw read-based Hi-C processing pipeline. `shortWithScore` is being used to carry aggregate counts.

## Planned final validation

- Run a small test `.cool` through the **full** pipeline and open the resulting `.hic`.
- Compare contacts at known genomic bin pairs in Cooler and Juicer at 400 bp.
- Check total cis counts and verify how diagonals are represented.
- Confirm chromosome naming, norms and available resolutions.
- Document the exact original source of the Micro-C dataset, its accession and original processing pipeline.

## References

- [Original public dataset GSE265818 at GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE265818)
- [Cooler CLI documentation](https://cooler.readthedocs.io/en/latest/cli.html)
- [Juicer Tools: `pre` file formats and options](https://github-wiki-see.page/m/aidenlab/juicer/wiki/Pre)
- [UCSC dm6 chromosome sizes](https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes)

## Portfolio note

Skills illustrated: handling public 3D genomics datasets, command-line data processing, replicate matrix merging, resolution selection, Python scripting, HDF5/COOL inspection, chromosome naming harmonization, data-format interoperability, Juicer Tools, region-level matrix extraction, Matplotlib visualizations, logarithmic data display, and workflow documentation. All claims of successful `.hic` generation from real data should wait for the end-to-end validation above.
