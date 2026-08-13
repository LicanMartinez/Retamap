# Retamap — Mapping *Cytisus scoparius* invasion in northern Patagonia with Sentinel-2 and Random Forest

This repository contains the Google Earth Engine (GEE) and R code used to detect
and map the invasive shrub *Cytisus scoparius* (Scotch broom) across northern
Patagonia (Argentina) from annual **Sentinel-2** mosaics over **2017–2025**, and
to assess the accuracy of the resulting maps.

The approach exploits the strong phenological contrast of the species: its
conspicuous **yellow flowering** in late spring (November–December, high NDYI)
versus a non-flowering late-summer reference (February). A **two-stage Random
Forest** classifier is used: a first model (RF1) trained on photo-interpreted
reference polygons, and a second model (RF2) trained on temporally **stable
samples** derived from the multi-year RF1 series.

> **Sentinel-2 note:** Processing Baseline 04.00 (scenes after 2022-01-25)
> introduced a +1000 DN offset. The pipeline applies **no** correction of its own:
> the `_HARMONIZED` collection already shifts post-baseline scenes back onto the
> pre-baseline scale.

## Citation

...

## Where each part of the paper lives

| Manuscript section | Code |
|---|---|
| 2.2.1 Satellite imagery and mosaic construction | `01_sentinelMosaic.js` |
| 2.2.2 Reference data and first-stage Random Forest (RF1) | `02_RF1fit.js` |
| 2.2.3 Stable-sample definition | `03_stableRegions.js` |
| 2.2.4 Second-stage Random Forest (RF2) | `04_RF2sampleExport.js` → `05_RF2train.js` → `06_RF2predict.js` |
| 2.2.5 Post-classification filtering | `07_RF2gapFill.js` → `08_RF2patchFilter.js` |
| 2.2.6 Validation | `validation/` (see below) |
| 2.3 / 2.4 Occupation, drivers and management | `tools/12_gridExport.js` builds the 1 km grid; the analysis itself is in `r_analysis/` |
| Reported invaded areas (ha) | `tools/13_areaAudit.js` |

## Pipeline overview

Scripts are sequential — each consumes the GEE assets exported by the previous
one, and **every arrow is a GEE batch task that has to finish before the next
script is run**. Assets live under `projects/ee-licanemartinez/assets/Retamap/`;
the numeric prefix of each asset matches the script that creates it.

| # | Script | Purpose | Main output asset |
|---|--------|---------|-------------------|
| 01 | `01_sentinelMosaic.js` | Annual Sentinel-2 mosaics: Cloud Score+ masking (cs ≥ 0.6), NDYI, water & elevation (≤ 1100 m) masks; Nov–Dec NDYI quality mosaic + Feb median, merged into a 10-band stack | `01_MergedBands_YYYY` |
| 02 | `02_RF1fit.js` | First-stage Random Forest (RF1): trains on the reference polygons plus the complementary control polygons, predicts per year, applies a connectivity filter | `02_RF1_2compCTRL_prediction_connectedFilter_YYYY` |
| 03 | `03_stableRegions.js` | Temporal-stability classes from the 9-year RF1 stack (stable / quasi-stable retama and background) | `03_RF2_stable_categories_rf1compCtrl` |
| 04 | `04_RF2sampleExport.js` | Samples the stable classes once, splits 80/20, extracts spectral features year by year | `04_RF2_samples_YYYY<SUFFIX>_{train,val}` |
| 05 | `05_RF2train.js` | Trains the single temporally generalized RF2 on the pooled `_train` samples | `05_RF2_classifier<SUFFIX>` |
| 06 | `06_RF2predict.js` | Applies the RF2 classifier to every annual mosaic (raw predictions) | `06_RF2_raw_prediction_YYYY<SUFFIX>` |
| 07 | `07_RF2gapFill.js` | Temporal gap-fill (8-of-9-year consistency rule) | `07_RF2_gapFill_YYYY<SUFFIX>` |
| 08 | `08_RF2patchFilter.js` | Connectivity filter → **final** annual maps | `08_RF2_prediction_YYYY<SUFFIX>` |

Run order: **01 → 02 → 03 → 04 → 05 → 06 → 07 → 08**. Everything in
`validation/` runs after stage 08; everything in `tools/` is auxiliary and can be
run at any time.

> **One deliberate loop in the numbering.** `02_RF1fit.js` reads an asset called
> `04_compBackgroundPolys`, which is exported by section 6 of
> `04_RF2sampleExport.js`. This is not an accident of naming: inspecting the
> stable samples of step 03 revealed background units that RF1 kept calling
> *C. scoparius*, so complementary control polygons were digitized over them and
> fed back into RF1 (the single refinement iteration described at the end of
> manuscript section 2.2.2). On a cold start, run section 6 of script 04 once to
> create that asset, or set `INCLUDE_COMP_POLYS = false` in script 02 to
> reproduce the unrefined RF1.

> **Note on stage-08 outputs.** `08_RF2_prediction_YYYY<SUFFIX>` is a
> **presence-only** raster: retama pixels hold 1 and everything else — including
> background — is *nodata*, not 0. Any area statistic must therefore reintroduce
> the valid-observation footprint from `01_MergedBands_YYYY`, and must be measured
> at the native 10 m resolution: coarsening a presence-only mask makes Earth
> Engine read it from an image pyramid where a single 10 m pixel fills the whole
> coarse cell. Exports also use `scale: 10` without an explicit `crs`, so the
> assets are in EPSG:4326, where the pixel measures ~75.1 m², not 100 m².
> `tools/13_areaAudit.js` exists to keep both of those honest.

## Canonical configuration

This is the authoritative configuration that produced the published results
(`RUN_SUFFIX = _trimmedRF1comp_s1n10k_qs1n10k_s0n15k_qs0n15k`):

- **Imagery:** Sentinel-2 harmonized **Level-1C / top-of-atmosphere**
  (`COPERNICUS/S2_HARMONIZED` — the surface-reflectance archive would be
  `COPERNICUS/S2_SR_HARMONIZED`), 10 m; scenes > 40 % cloud discarded; Cloud
  Score+ clear-sky threshold 0.6. Ten-band feature stack: B2, B3, B4, B8 + NDYI
  for each of the two seasonal composites. After the water and elevation masks
  the classified area is **950,119 ha**, 38 % of the study-area polygon.
- **Reference polygons:** 87 *C. scoparius* polygons (72.3 ha; median 0.32 ha)
  and 150 control polygons (2,047 ha; median 3.9 ha), digitized over the seven
  years with usable flowering-period imagery (2017–2023).
- **RF1:** 500 trees (`smileRandomForest`); stratified sampling of up to
  **10,000 retama** and **15,000 control** pixel locations. Retama samples are
  drawn only in each polygon's digitization year while control samples are drawn
  across all nine years, so the fitted table holds roughly **10,000 retama and
  135,000 control rows**. Connectivity filter: class-1 patches of fewer than
  **5** contiguous pixels are reclassified as background (patch extent evaluated
  up to a search limit of 50 pixels).
- **Stable samples:** four temporal-stability classes over 9 years — stable
  retama (9/9), quasi-stable retama (8/9), stable background (0/9), quasi-stable
  background (1/9). Everything else is discarded as ambiguous.
- **RF2:** 500 trees, same 10-feature stack and hyperparameters as RF1; a single,
  temporally generalized classifier. **50,000 pixel locations are drawn once**
  (10,000 stable retama + 10,000 quasi-stable retama + 15,000 stable background
  + 15,000 quasi-stable background), split 80/20, and the spectral features of
  the **40,000 training locations are extracted from each of the nine mosaics
  (≈ 360,000 training rows)**. Those rows are nine temporal replicates of the
  same sites, not independent samples. The 20 % partition is held out but is
  **not** used to assess accuracy.
- **Post-classification:** (1) temporal gap-fill using the 8-of-9-year
  consistency rule; (2) connectivity filter with a minimum patch size of
  **3** contiguous pixels.

## Validation

Accuracy is assessed against an **independent reference dataset** built by blind
human photo-interpretation, on a sample **stratified by the final map itself**,
following Olofsson et al. (2014). The first and last years of the series (2017
and 2025) are validated with the same set of points.

Three strata partition the classified area:

| Stratum | Definition | Area | Points |
|---|---|---|---|
| **S1** | mapped as retama in 2017 and/or 2025 | 1,244 ha (0.13 %) | 400 |
| **S2** | background within ~1 km of mapped retama | 70,344 ha (7.4 %) | 400 |
| **S3** | the remaining background | 878,530 ha (92.5 %) | 119 |

Reference labels are collected through a blind GEE App (`validation/09v_valApp.js`)
that walks each interpreter through their assigned points over the raw annual
mosaics and NDYI/NDVI time series — it never displays the classification —
with capture in per-interpreter Google Sheets. The point set is served as the
blind asset `09v_valPoints_fc` (point IDs only). S1 points are assigned to two
interpreters, S2 and S3 points to one.

Headline results (crude, with 95 % Wilson intervals): overall accuracy
**92.6 %** in 2017 and **92.9 %** in 2025; commission of the retama class ~9 % in
both years; omission 12.1 % and 16.0 %. Because the allocation is deliberately
disproportionate, these describe accuracy **within the mapped patches and their
surroundings**, not over the whole classified area — see the manuscript for the
full statement. Full tables, confidence intervals, the area-weighted estimators
and the HTML report are produced by `validation/09v_valAnalysis_compute.R` and
`09v_valAnalysis_summarized.Rmd`.

### Validation file map

| File | Role |
|---|---|
| `09v_valGenerator.js` | GEE — builds the three strata, computes the Olofsson sample sizes and the stratum weights, exports the raw point CSV |
| `09v_valRandomize.R` | R — shuffles and assigns the points, writes the blind FeatureCollection, the per-interpreter sheets and KMLs, the local master key, and regenerates `09v_valAssign.js` |
| `09v_valAssign.js` | **Generated data, not code**: point IDs per interpreter. Loaded by the App through `require()` |
| `09v_valApp.js` | GEE App — the blind interpretation interface |
| `09v_valSheetsIngest.R` | R — turns the pasted sheets into canonical CSVs and reconciles them against the master key |
| `09v_valAnalysis_compute.R` | R — the metrics engine (confusion matrices, Wilson intervals, Olofsson estimators, inter-observer agreement, figures) |
| `09v_valAnalysis_summarized.Rmd` | R — the report |
| `09v_valMapLabels.js` + `09v_valCompare.R` / `.Rmd` | Map-agnostic tooling to score **any** map version on the frozen validation points, with a paired McNemar test |
| `09v_valKmlRefresh.R` | Utility — regenerates only the KMLs, without re-shuffling the assignment |

> **Do not re-run `09v_valGenerator.js`** to evaluate a different version of the
> map. Its stratified sample would draw different points and orphan the reference
> labels already collected. Use `09v_valMapLabels.js` + `09v_valCompare.R`, which
> exist for exactly that.

## Repository structure

```
.
├── 01_sentinelMosaic.js … 08_RF2patchFilter.js   # the pipeline that produces the published maps
├── validation/                                   # independent accuracy assessment (see above)
├── tools/                                        # 11, 11.2 visualization · 12 grid · 13 area audit
├── r_analysis/                                   # downstream R analysis (occupation, drivers, management)
├── explorations/                                 # side explorations: run, validated, NOT adopted
├── experimental/                                 # scratch forks, never evaluated (untracked)
├── validation_data/                              # local working data of the R scripts (untracked)
└── README.md
```

Three different kinds of "not part of the pipeline", which are easy to confuse:

- **`explorations/`** — branches that were carried through to a measured
  validation result and are documented with their verdict. Versioned, so they
  stay reproducible. Currently `yearlyRF2/` (one Random Forest per year instead
  of the pooled one — rejected: it worsens omission in both validated years).
- **`experimental/`** — scratch forks that were never evaluated. Untracked.
- **`validation_data/`** — the working data of the validation R scripts. Untracked
  **on purpose**: the master key carries the reference labels, and publishing it
  would break the blindness of any future interpretation round.

`r_analysis/` currently holds only its README: the downstream statistical
analysis is maintained separately.

## How to run

The GEE scripts are written for the **Code Editor** (JavaScript API). Each
declares an `ASSET_PREFIX` and a `RUN_SUFFIX` in its config section. Run them in
the order shown above; each stage exports assets through GEE batch tasks that
must complete before the next stage starts. The `EXPORT_TO_DRIVE` toggle
(default `false`) additionally writes outputs to a Google Drive folder.

Reproducing the classification requires the reference asset
`gt_polys_6_ctrlReduce_moreCtrls_refineRetamas` (the manually digitized
polygons), the study-area polygon `3_study_area_retama`, and a GEE account with
access to the `Retamap` asset folder.

The R scripts currently carry an absolute `PROJECT_DIR` in their config section
and read from `validation_data/`, which is not distributed. They are published as
a record of how the reported numbers were computed rather than as a turnkey
pipeline; adapt `PROJECT_DIR` and supply the reference labels to re-run them.

## Data & assets

All pipeline assets live under `projects/ee-licanemartinez/assets/Retamap/`, with
the numeric prefix of each asset matching the script that creates it. The
independent validation point set is served as the blind asset
`09v_valPoints_fc` (point IDs only).

## License

`[to be added]`
