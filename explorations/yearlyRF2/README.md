# `yearlyRF2` — one RF2 per year instead of one pooled RF2

**Status: explored, measured, NOT adopted (2026-08-09).** The canonical pipeline
(`05_RF2train.js` → `08_RF2patchFilter.js`) is unchanged and remains the source
of the published maps. These scripts are kept because the question may be worth
re-opening with a different design — see *If this is ever revisited* below.

## The question

Canonical stage 05 pools the nine per-year `_train` sample assets (~360k
features) into a **single** RF2 and uses it for every year. This variant instead
trains **one RF2 per year**: the sampling *locations* are identical across years
(stage 04 runs `stratifiedSample` once, then `sampleRegions` per year), so the
only thing that differs between the per-year models is that year's spectral
signal. The hypothesis was that year-specific phenology, atmospheric conditions
and mosaic quality would be fit better by a year-specific model.

Run for **2017 and 2025 only** — the two years covered by the independent
validation.

## Scripts

| Script | Produces |
|---|---|
| `05y_RF2trainPerYear.js` | `05y_RF2_classifier_<year>_yearlyRF2` |
| `06y_RF2predictPerYear.js` | `06y_RF2_raw_prediction_<year>_yearlyRF2` |
| `08y_RF2patchFilter.js` | `08y_RF2_prediction_<year>_yearlyRF2_noGapFill` **and** `08y_RF2_prediction_<year>_globalRF2_noGapFill` |

Two suffixes per script, on purpose: `SAMPLES_SUFFIX` reads the **existing**
canonical `04_` sample assets (stage 04 is never re-run) and
`RUN_SUFFIX = '_yearlyRF2'` writes the new ones. `BANDS` (the exact 10, order
matters) and `numberOfTrees: 500` are held identical to the canonical run so the
comparison isolates the training change.

The assets still exist in GEE under `projects/ee-licanemartinez/assets/Retamap/`.

### Why there is no `07y` gap-fill, and what the second output is for

Stage 07 is a **temporal** filter: `sumOnes == N-1 → 1` and `sumOnes == 1 → 0`.
With `N = 2` both thresholds collapse onto `1`, so every pixel where 2017 and
2025 disagree gets *swapped* — the change signal is destroyed. Gap-fill has no
meaningful two-year analogue, so this variant skips it.

That would have made a comparison against the published map conflate two changes
(per-year RF2 *and* no gap-fill), so `08y` also runs the **same**
`applyPatchFilter` function over the canonical `06_` raw predictions, producing a
matched `globalRF2_noGapFill` baseline. The real A/B is
`globalRF2_noGapFill` vs `yearlyRF2_noGapFill`.

## What the validation said

Measured with the independent 09v reference (1130 human labels, 919 points, same
points and same labels for every arm — a **paired** comparison), using the
map-agnostic comparison tooling that stayed in the active pipeline:
`../../validation/09v_valMapLabels.js` + `../../validation/09v_valCompare.R` / `.Rmd`.

Figures below are from the **2026-08-11** re-run: reconciliation rules
`MIN_CONF = 5` + `OR_THRESHOLD = 2`, n = 826 / 813 evaluated point-years, and the
**corrected stratum weights** (measured at 10 m in UTM19S — see the note at the
end of this section). Re-running `09v_valCompare.R` after changing those knobs
will shift the figures; the direction and the verdict have held across every rule
version tried so far.

Both arms are compared **without gap-fill**, which is the paired baseline the
script uses (`08y_..._globalRF2_noGapFill`) — the yearly variant cannot go
through script 07, so comparing it against the published map would confound the
two changes.

| | 2017 global → yearly | 2025 global → yearly |
|---|---|---|
| **Retama omission** (crude) | 12.8% → **20.6%** | 16.0% → **20.7%** |
| **Retama omission** (area-weighted) | 54.0% → **58.1%** | 46.7% → **58.7%** |
| Retama commission (crude) | 10.2% → 7.4% *(better)* | 9.5% → 10.9% |
| Overall accuracy (unweighted) | 92.3% → 90.8% | 92.7% → 91.1% |
| McNemar (paired) | p = 0.097 | **p = 0.016** |

Two reasons it was rejected:

1. **It trades good omission for marginal commission.** The flow table shows the
   mechanism: of the 40 points it changes in 2017, **35 are retama it removes** —
   and only 12 of those removals were correct. It adds 23 false negatives to
   remove 12 false positives. Omission was already the map's weakest metric.
   2025 is the same story: 16 removals in S1, only 2 correct.
2. **It breaks the zero far-background commission property.** In 2025 the yearly
   map puts a false positive in the far-background (S3) stratum; the published
   global map has none in 420 point-years. Because S3 carries 92.5% of the ROI
   area, that single point drives the area-weighted retama commission to **91.2%
   [75.4, 107]**. Numerically that rests on **one** sampled point and the
   interval is enormous — but structurally it is the more serious of the two
   findings, since the far-commission-is-zero property is what the whole
   stratified design leans on.

> **Weights corrected 2026-08-11.** Every area-weighted figure on this page was
> recomputed after the stratum areas were found to be measured at 100 m over a
> presence-only raster, which inflated S1 by 4.4× (see `../../tools/13_areaAudit.js`
> and `docs/area_discrepancy_note.md`). The correction raised the weight of S2
> relative to S1 by ~5×, so *all* weighted omission figures — for the published
> map as much as for this variant — are substantially higher than previously
> reported. **The verdict is unaffected**: it rests on crude omission, the
> McNemar test and the flow table, none of which use the weights.
>
> One thing the re-run did surface: the *global* map without gap-fill also
> carries a heavy-stratum false positive in 2025 (weighted commission 20.3%,
> against 9.1% for the published map). That is a point in favour of the
> gap-fill, not of the yearly variant, and it is worth a closer look if the
> gap-fill question is reopened.

### Why this was, in hindsight, predictable

The training labels come from `03_RF2_stable_categories`, which is a **temporal**
construct: a pixel is "stable retama" because of how it behaves across nine
years. The pooled model sees each point nine times, with nine different spectra
and the same label, which forces it to learn a **year-invariant** signature. The
per-year model sees each point once and can fit that year's noise (a shadow, a
mosaic artefact). Pooling across years is not merely "more data" here — it *is*
the regularization mechanism of the method.

### Side finding worth keeping

`canon` (with gap-fill) vs `globalNoGF` barely differ **in how many points they
get right** (2017: 765 vs 762 of 826; 2025: 755 vs 754 of 813). Counting hits,
gap-fill looks like it contributes almost nothing.

Weighting by area tells a different story: in 2025 the no-gap-fill map carries a
false positive in a heavy stratum, so its weighted retama commission is 20.3%
against 9.1% for the published map. **A change that moves one point can move a
weighted metric by ten points** — which is the whole reason the design is
stratified, and a reminder not to judge a map version by its hit count. Worth
following up on its own.

## Caveats on the verdict

- Only **two years** were tested. A year with a poor mosaic might favour a
  year-specific model; neither of these two did.
- At the time of the test **48.5% of Lican's validation load was still unloaded**
  (201/390), concentrated in stratum S1 — which is exactly where the difference
  between arms lives. The deficit would have to move a long way to reverse an
  8-point omission gap, but it is not a closed question.
- The validation sample was stratified by the **canonical** map, so it has little
  power to measure commission that a variant introduces outside S1. That is
  precisely why the single far-background false positive matters more than its
  raw count suggests.

## If this is ever revisited

Don't re-run this as-is. The mechanism above says a per-year model is strictly
less regularized for the same label set, so the promising directions are the
intermediate ones:

- pooled training with **year as a feature**, or a year-interaction term;
- an **ensemble** of the global and per-year models;
- per-year models with the *pooled* model as a prior / warm start.

And first: finish loading the validation, since S1 is where the signal is.

To measure any of them, add an entry to `MAPS` in `../../validation/09v_valMapLabels.js` and
one to `CMP_VARIANTS` in `../../validation/09v_valCompare.R`; the comparison report
re-renders with no text changes.

## Running these scripts again

They are exploration-only and are **not** part of the daily workflow, but they
are still in the GEE repo, so they appear in the Code Editor under
`Retamap/explorations/yearlyRF2/`. Order: `05y` → `06y` → `08y`, waiting for each
stage's assets to finish ingesting before running the next.
