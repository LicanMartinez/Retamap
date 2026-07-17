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
...

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
| 09v | `09v_valGenerator.js`, `09v_valApp.js` | Independent validation: map-stratified reference sampling (Olofsson) + blind photo-interpretation through a GEE App | `09v_valPoints_fc` |
| 11 | `11_predsVisualizer.js` | Visualization of RF1/raw/gap-fill/final predictions per year | — |
| 11.2 | `11.2_predsVisualizer_finalSeries.js` | Final-map time series (color gradient per year) | — |
| 12 | `12_gridExport.js` | Regular analysis grid over the ROI | `12_grid_<SIZE>m[_clipped]` |

Run order: **01 → 02 → 03 → 04 → 05 → 06 → 07 → 08**. Independent validation
(the `09v` family) runs after stage 08. Scripts 11/11.2 (visualization) and 12
(grid) are auxiliary.

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
  2017–2025)**. Four-class labels collapsed to binary before spectral
  extraction; the classifier is trained on the pooled training partition.
- **Post-classification:** (1) temporal gap-fill using the 8-of-9-year
  consistency rule; (2) connectivity filter with a minimum patch size of
  **3 connected pixels** (50-pixel neighborhood).

## Validation

Accuracy is assessed against an **independent reference dataset** built by human
photo-interpretation, using a **map-stratified sampling** design over the final
maps (stage C) following Olofsson et al. (2014). Two years spanning the series
are validated (2017 and 2025).

Three strata partition the study area, enabling design-based area and accuracy
estimation:

- **S1 — mapped retama:** pixels classified as *C. scoparius* in 2017 and/or 2025.
- **S2 — near non-retama:** non-S1 pixels within ~1 km of mapped retama
  (commission-prone zones).
- **S3 — far non-retama:** the remainder of the area, entering the estimators
  with its area weight.

Sample sizes per stratum are set from the target standard error of each stratum
(Olofsson design). Reference labels are collected through a blind GEE App
(`09v_valApp.js`) that navigates each validator's assigned points over the annual
mosaics, with per-validator capture in Google Sheets; the point set is served as
the blind asset `09v_valPoints_fc`. Results are in preparation.

## Repository structure

```
.
├── 01_sentinelMosaic.js … 08_RF2patchFilter.js   # canonical classification pipeline
├── 09v_val*.{js,R}                               # independent validation (generator, App, randomizer)
├── 09vNear_val*.{js,R}                            # near-only validation variant (2 strata)
├── 11_predsVisualizer.js, 11.2_…finalSeries.js   # visualization
├── 12_gridExport.js                              # analysis grid
├── r_analysis/                                   # downstream R analysis (occupation, drivers, management)
├── experimental/                                 # non-canonical exploratory variants
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
`projects/ee-licanemartinez/assets/Retamap/`, with the numeric prefix of each
asset matching the script that creates it. The independent validation point set
is served as the blind asset `09v_valPoints_fc` (point IDs only).

## License

`[to be added]`
