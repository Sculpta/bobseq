# Differential splicing of the 24-plex (Fig. 7C)

**Question.** Which of the plate's perturbations changes splicing genome-wide relative to the control wells: cycloheximide
(CHX, 1 and 50 µg/mL), risdiplam (25 nM) or the splice-switching ASO (1 and 100 nM)? Five contrasts, each on its own wells
(3 controls + 3 per treated condition): CHX dose trend (control / 1 / 50), CHX 50 vs control, Ris 25 vs control, ASO dose trend
(control / 1 / 100), ASO 100 vs control. Fig. 7C = `figures/volcano_chx_dose.svg`.

## Method (`src/ds.py`)

Unit = fragments (S6: a read pair spanning a junction counts once, no deduplication), PSI per well = junction fragments /
LSV-side fragments. Preselection, condition-blind: LSV side with >= 100 fragments on average and > 0 in every well of the
contrast, PSI variance across the wells >= 0.005. Test: pairwise contrasts, Welch t-test on the per-well PSI; dose contrasts,
OLS of the per-well PSI on dose rank 0 / 1 / 2 (9 wells, 7 residual df), t-test on the slope; Benjamini-Hochberg within contrast.
Effect: dPSI = PSI(top dose or treated) - PSI(control) from the pooled fragment counts (`dpsi_dedup` = the same on the
UMI-deduplicated molecules of S4). Hit = padj < 0.10 and dPSI >= 0.10, positive direction (every LSV contributes a +/- pair whose
members sum to 1; the rising member represents the event). Top 10 per contrast = the hits ranked by dPSI, one row per LSV and
per junction, padded with the largest-dPSI non-hits (`hit = 0`) when a contrast has fewer than 10.

## Result (`results/ds_summary.tsv`)

| contrast | test | rows tested (LSV sides) | padj < 0.1, any direction | hits: rows / LSVs / genes |
|---|---|---|---|---|
| CHX dose trend | OLS slope | 530 (283) | 53 | 28 / 28 / 25 |
| CHX 50 µg/mL vs control | Welch | 584 (311) | 10 | 6 / 6 / 6 |
| Ris 25 nM vs control | Welch | 425 (239) | 0 | 0 |
| ASO dose trend | OLS slope | 448 (244) | 0 | 0 |
| ASO 100 nM vs control | Welch | 432 (241) | 0 | 0 |

CHX top 10 by dPSI (dose trend): CLASRP +0.46 (padj 0.006), EIF3C +0.39 (0.088), AKAP8 +0.34 (0.007), HRAS +0.28 (0.092), PSMG4
+0.27 (0.009), SNHG29 +0.25, CENPX +0.24, ENSG00000233461 +0.22, TXNL4A +0.22, USF2 +0.22. Risdiplam 25 nM and the ASO give no
row at padj < 0.1 in either direction (SMN2 exon 7 cannot appear here: its skipping reads are SMN1/SMN2 multimappers removed by
the MAPQ filter; see `../smn2_exon7/`).

Tables: `ds_<contrast>.tsv` (every tested row: PSI per condition and well, dPSI, test, padj, hit), `ds_hits.tsv`, `ds_top10.tsv`
(+ `.md`, `.xlsx`), `barchart_psi.tsv` (the bar-chart numbers), `top10_loci.tsv` and `top10_junctions.bed` (the coordinates the
track figures use).

## Run

```
python3 run.py            # ds, volcano, barcharts, loci; about 1 min; numpy, scipy, matplotlib, openpyxl
```

Figures: `volcano_<contrast>.svg` (and `_labelled`: the top 10 named), `barcharts_<contrast>.svg` (top-10 PSI over every
condition but Ris 500 nM, pooled fragment PSI with Clopper-Pearson 95 % interval, wells as dots). Caveats: 3 vs 3 wells make the
pairwise tests weak (6 CHX hits against 28 in the dose model); the dose model treats dose as a rank; the junction universe is the
catalogue the LSVs are defined on.
