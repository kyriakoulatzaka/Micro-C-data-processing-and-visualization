# ems locus: nc1–8 vs nc14 Micro-C figures

These figures come from the author's **real-data reanalysis** of publicly available *Drosophila melanogaster* embryonic Micro-C contact matrices (GSE265818), not from the synthetic test matrices used to unit-test the revised scripts.

## Expected image files

| Filename | Developmental stage | Resolution | Window annotation |
|---|---|---|---|
| `ems_nc1-8_400bp_viewA.png` | nc1–8 | 400 bp bins | original filename described approximately 150 kb; exact genomic coordinates to confirm |
| `ems_nc1-8_400bp_viewB.png` | nc1–8 | 400 bp bins | user reported 166,080 bp, exact genomic coordinates to confirm |
| `ems_nc14_400bp_viewA.png` | nc14 | 400 bp bins | original filename described approximately 150 kb; exact genomic coordinates to confirm |
| `ems_nc14_400bp_viewB.png` | nc14 | 400 bp bins | user reported 166,080 bp, exact genomic coordinates to confirm |

**Do not confuse** the 400-bp bin resolution with the genomic window width. A separate preserved plotting command used `region = '3R:13778677-14027796'` (249,120 bp inclusive) for an nc14 image. The exact association of that command with these four PNG versions still needs to be checked.

The previously shared heatmap script used `balance=False` and `LogNorm` with the `YlOrBr` palette. Confirm that the same settings were used for both plotted stages before asserting normalized or quantitative differences.

**Source:** Public Micro-C data, NCBI GEO [GSE265818](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE265818). These plots reflect the author's reanalysis of published data. Formal differential-contact analysis has not been claimed.

**Upload status:** The original four PNG files need to be added to this folder to render the portfolio gallery.
