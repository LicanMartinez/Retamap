# Downstream R analysis

This folder holds the R analysis built on top of the GEE classification outputs
(`08_RF2_prediction_YYYY<SUFFIX>`): occupation/coverage trends, environmental and
landscape drivers, and the management-prioritization workflow described in the
paper's Methods §2.3–2.4.

## Scope

- **Occupation and coverage analysis** — annual area occupied by *C. scoparius*
  by subarea type (urban, National Parks, roads, rivers), fitted with a
  Generalized Additive Model (`mgcv`).
- **Drivers of presence and coverage** — spatial generalized linear models
  (`spmodel`) relating *C. scoparius* presence/coverage in a 1,000 × 1,000 m grid
  to environmental (temperature, precipitation, elevation), corridor (roads,
  rivers, lake shoreline), and landscape (urban area, National Park area,
  historic fires) predictors. Residual spatial autocorrelation checked with
  variograms (`gstat`) and Moran's I (`spdep`).
- **Coverage change model** — spatial model of the change in *C. scoparius*
  coverage between the first and last three years of the series, related to
  change in urban area.
- **Management index** — landscape invasion metrics (`landscapemetrics`: number
  of invaded patches, total invaded area, largest patch index) per grid cell,
  classified into a four-level management index (No invasion / Eradication /
  Containment / Restoration).

## Inputs

Scripts in this folder consume the final annual classification rasters produced
by the canonical GEE pipeline (`08_RF2_prediction_YYYY<SUFFIX>`, see the
repository [README](../README.md)) together with external layers (roads,
rivers, urban cover, National Parks, fire history, climate).

## Maintainer

This part of the repository is maintained by Sofía (paper lead author).
