# Retamap — Mapping *Cytisus scoparius* invasion in northern Patagonia with Sentinel-2 and Random Forest

This repository contains the Google Earth Engine (GEE) methods used to detect and
map the invasive shrub *Cytisus scoparius* (Scotch broom) across northern
Patagonia (Argentina) from annual **Sentinel-2** mosaics over **2017–2025**.

The approach exploits the strong temporal/phenological contrast of the species:
its conspicuous **yellow flowering** in late spring (November–December, high NDYI)
versus a non-flowering late-summer reference (February). A **two-stage Random
Forest** classifier is used: a first model (RF1) trained on photo-interpreted
reference polygons, and a second model (RF2) trained on temporally **stable
samples** derived from the multi-year RF1 series.

> **Sentinel-2 note:** Processing Baseline 04.00 (scenes after 2022-01-25)
> introduced a +1000 DN offset; bands B2/B3/B4/B8 are corrected by subtracting
> 1000.

## Citation



If you use this code, please cite the paper above.

## Pipeline overview

Scripts are sequential — each consumes the GEE assets exported by the previous
one. Assets use the prefix `projects/ee-licanemartinez/assets/Retamap/`; the
numeric prefix of each asset matches the script that creates it.

| # | Script | Purpose | Main output asset |
|---|--------|---------|-------------------|
| 01 | `01_sentinelMosaic.js` | Annual Sentinel-2 mosaics: Cloud Score+ masking (cs ≥ 0.6), NDYI, water & elevation (≤1100 m) masks; Nov–Dec NDYI quality mosaic + Feb median, merged into a 10-band stack | `01_MergedBands_YYYY` |
| 02 | `02_RF1fit.js` | First-stage Random Forest (RF1): trains on reference + complementary control polygons, predicts per year, applies a connectivity filter | `02_RF1_2compCTRL_prediction_connectedFilter_YYYY` |
| 03 | `03_stableRegions.js` | Temporal-stability classes from the 9-year RF1 stack (stable/quasi-stable retama & background) | `03_RF2_stable_categories_rf1compCtrl` |
| 04 | `04_RF2sampleExport.js` | Samples stable pixels per year, splits 80/20 train/val, extracts spectral features | `04_RF2_samples_YYYY<SUFFIX>_{train,val}` |
| 05 | `05_RF2train.js` | Trains the single temporally generalized RF2 on pooled `_train` samples | `05_RF2_classifier<SUFFIX>` |
| 06 | `06_RF2predict.js` | Applies the RF2 classifier to every annual mosaic (raw predictions) | `06_RF2_raw_prediction_YYYY<SUFFIX>` |
| 07 | `07_RF2gapFill.js` | Temporal gap-fill (8-of-9-year consistency rule) | `07_RF2_gapFill_YYYY<SUFFIX>` |
| 08 | `08_RF2patchFilter.js` | Connectivity filter → **final** annual maps | `08_RF2_prediction_YYYY<SUFFIX>` |
| 09 | `09_holdoutValidation.js` | Holdout validation on `_val` samples × 3 pipeline stages → metrics CSV | `Holdout_Summary_Metrics<SUFFIX>.csv` |
| 11 | `11_predsVisualizer.js` | Visualization of RF1/raw/gap-fill/final predictions per year | — |
| 11.2 | `11.2_predsVisualizer_finalSeries.js` | Final-map time series (color gradient per year) | — |
| 12 | `12_gridExport.js` | Regular analysis grid over the ROI | `12_grid_<SIZE>m[_clipped]` |

Run order: **01 → 02 → 03 → 04 → 05 → 06 → 07 → 08**. Script 09 (validation) can
run after 08. Scripts 11/11.2 (visualization) and 12 (grid) are auxiliary.

## Canonical configuration

This is the authoritative configuration that produced the published results
(`RUN_SUFFIX = _trimmedRF1comp_s1n10k_qs1n10k_s0n15k_qs0n15k`):

- **Imagery:** Sentinel-2 L2A harmonized (`COPERNICUS/S2_HARMONIZED`), 10 m;
  scenes > 40 % cloud discarded; Cloud Score+ clear-sky threshold 0.6. Ten-band
  feature stack: B2, B3, B4, B8 + NDYI for each of the two seasonal composites.
- **RF1:** 500 trees (`smileRandomForest`); stratified sampling of
  **10,000 retama** + **15,000 control** pixels; connectivity filter retaining
  class-1 patches of ≥ 5 connected pixels (50-pixel neighborhood). Retama samples
  drawn only from each polygon's digitization year; control samples drawn across
  all years.
- **Stable samples:** four temporal-stability classes over 9 years — stable
  retama (9/9), quasi-stable retama (8/9), stable background (0/9), quasi-stable
  background (1/9).
- **RF2:** 500 trees, same 10-feature stack and hyperparameters as RF1; a single,
  temporally generalized classifier. Sampling per year:
  **10,000 stable retama + 10,000 quasi-stable retama + 15,000 stable background
  + 15,000 quasi-stable background = 50,000 px yr⁻¹ (450,000 px total over
  2017–2025)**. Four-class labels collapsed to binary; 80/20 train/validation
  split with a fixed seed before spectral extraction.
- **Post-classification:** (1) temporal gap-fill using the 8-of-9-year
  consistency rule; (2) connectivity filter with a minimum patch size of
  **3 connected pixels** (50-pixel neighborhood).

## Validation

Accuracy was assessed by **holdout validation** on the 20 % of stable-sample
pixels withheld from RF2 training, evaluated at three pipeline stages
(A: raw RF2, B: gap-filled, C: connectivity-filtered final map). Full per-year ×
per-stage metrics are in [`validation/`](validation/), along with the R analysis
([`Retamap_validation.Rmd`](validation/Retamap_validation.Rmd), rendered HTML and
plots).

Final-map (stage C) performance, mean ± SD across 2017–2025:

| Metric | Value |
|--------|-------|
| Overall accuracy | 0.989 ± 0.003 |
| F1-score | 0.987 ± 0.003 |
| Cohen's κ | 0.978 ± 0.006 |
| Producer's accuracy (*C. scoparius*) | 0.979 ± 0.006 |
| User's accuracy (*C. scoparius*) | 0.994 ± 0.003 |
| Omission error (*C. scoparius*) | 0.021 ± 0.006 |
| Commission error (*C. scoparius*) | 0.006 ± 0.003 |

## Repository structure

```
.
├── 01_sentinelMosaic.js … 12_gridExport.js   # canonical pipeline + viz/grid
├── validation/                               # definitive holdout validation
│   ├── Holdout_Summary_Metrics_trimmedRF1comp_s1n10k_qs1n10k_s0n15k_qs0n15k.csv
│   ├── Retamap_validation.Rmd                # R analysis (self-contained)
│   ├── Retamap_validation.html               # rendered report
│   ├── validation_plots_byStage.png
│   ├── validation_plots_byYear.png
│   └── metricas.txt                          # metric definitions
├── r_analysis/                               # downstream R analysis (occupation, drivers, management)
└── README.md
```

## How to run

The scripts are written for the **GEE Code Editor** (JavaScript API). Each script
declares an `ASSET_PREFIX` and a `RUN_SUFFIX` in its config section. Run the
scripts in the order shown above; each stage exports its assets via GEE Batch
tasks that must complete before the next stage is run. The `EXPORT_TO_DRIVE`
toggle (default `false`) additionally exports outputs to a Google Drive folder.

Reproducing the analysis requires the input reference asset
`gt_polys_6_ctrlReduce_moreCtrls_refineRetamas` (manually digitized reference
polygons) and a GEE account with access to the `Retamap` asset folder.

## Data & assets

All pipeline assets live under
`projects/ee-licanemartinez/assets/Retamap/`. The R validation analysis can be
reproduced from the CSV in [`validation/`](validation/) (R ≥ 4.3;
packages: `dplyr`, `tidyr`, `ggplot2`, `readr`, `patchwork`).

## License

`[to be added]`
