# =============================================================================
# 09v_valCompare.R
# -----------------------------------------------------------------------------
# Compare SEVERAL map versions against the SAME validation points and the SAME
# human reference labels, by running 09v_valAnalysis_compute.R once per version.
#
# Consumes:
#   - validation/09v_valMapLabels.csv   (from 09v_valMapLabels.js, via Drive)
#         pxid + one column m<year>_<tag> per map version
#   - validation/09v_master_key.csv     (strata, validators, frozen map labels)
#   - validation/sheets_filled/<val>_filled.csv   (reference labels)
#
# Produces (in validation/analysis/):
#   - one subfolder per variant (<tag>/) with that variant's full CSV/figure set
#   - 09v_cmp_metrics_crude.csv       side-by-side crude metrics + Wilson CIs
#   - 09v_cmp_metrics_olofsson.csv    side-by-side area-weighted metrics
#   - 09v_cmp_area.csv                side-by-side retama area estimates
#   - 09v_cmp_mcnemar.csv             paired test, baseline vs test arm
#   - 09v_cmp_flow.csv                which points changed label, and were they right
#   - 09v_cmp_control.csv             control-column verification
#
# WHY THIS IS A VALID COMPARISON
# The point set, the strata, the area weights Wi and the human labels are all
# frozen properties of the sampling DESIGN. Evaluating a different map on them
# stays unbiased, and because every arm is scored on the SAME points with the
# SAME references, the comparison is PAIRED — which is what makes McNemar the
# right significance test rather than comparing two independent CIs.
#
# The catch, stated in the report: power is uneven. Retama an arm adds OUTSIDE
# the canonical S1 is only sampled through S2 (400 pts) and S3 (119 pts covering
# 93.5% of the ROI), so a new arm's commission in those areas carries very wide
# intervals.
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

# ---- 0. CONFIG --------------------------------------------------------------
PROJECT_DIR <- "D:/Lican/uni/Investigacion/Colaboraciones_y_ayudas/Retamap"
GEE_DIR     <- file.path(PROJECT_DIR, "gee_scripts")
WORK_DIR    <- file.path(GEE_DIR, "validation")
OUT_DIR     <- file.path(WORK_DIR, "analysis")
COMPUTE_R   <- file.path(GEE_DIR, "09v_valAnalysis_compute.R")
MASTER_KEY  <- file.path(WORK_DIR, "09v_master_key.csv")
# Allow a caller (or a smoke test) to point at a different labels file.
if (!exists("LABELS_CSV", inherits = FALSE)) {
  LABELS_CSV <- file.path(WORK_DIR, "09v_valMapLabels.csv")
}

# One entry per map version. `cols` are the columns of LABELS_CSV holding that
# version's map label for each year. Adding the future 9-year yearly run with
# gap-fill = one more entry here plus one more entry in 09v_valMapLabels.js.
CMP_VARIANTS <- list(
  list(tag = "canon",
       label = "Global RF2 + gap-fill (mapa publicado)",
       cols  = c("m2017_canon", "m2025_canon")),
  list(tag = "globalNoGF",
       label = "Global RF2, sin gap-fill",
       cols  = c("m2017_globalNoGF", "m2025_globalNoGF")),
  list(tag = "yearlyNoGF",
       label = "Yearly RF2, sin gap-fill",
       cols  = c("m2017_yearlyNoGF", "m2025_yearlyNoGF"))
)

# The A/B pair the experiment is actually about: same everything except RF2
# training. `canon` rides along as context (it shows the gap-fill effect).
BASELINE_TAG <- "globalNoGF"
TEST_TAG     <- "yearlyNoGF"
# Control column: re-extracted from the very asset that produced the frozen
# m2017/m2025 of the master key, so the two must agree row by row.
CONTROL_TAG  <- "canon"

CMP_YEARS <- c(2017, 2025)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. LOAD THE MAP LABELS -------------------------------------------------
if (!file.exists(LABELS_CSV)) {
  stop("Missing ", LABELS_CSV, "\n",
       "Run 09v_valMapLabels.js in the GEE Code Editor and download the CSV ",
       "from Drive (", "Retamap/Retamap_GEE_Exports", ") into validation/.")
}
lab <- read.csv(LABELS_CSV, stringsAsFactors = FALSE, colClasses = "character")
stopifnot("pxid" %in% names(lab), !anyDuplicated(lab$pxid))

# Blank (pixel masked out by the patch filter, or outside the map) => 0, the
# same convention the frozen master-key columns already use.
zero_blank <- function(x) ifelse(is.na(x) | x == "", 0L, as.integer(x))

need <- unlist(lapply(CMP_VARIANTS, `[[`, "cols"))
miss <- setdiff(need, names(lab))
if (length(miss)) stop("LABELS_CSV lacks column(s): ", paste(miss, collapse = ", "))

mk_cmp <- read.csv(MASTER_KEY, stringsAsFactors = FALSE, colClasses = "character")
stopifnot(setequal(mk_cmp$pxid, lab$pxid))

# ---- 2. CONTROL-COLUMN VERIFICATION (blocking) ------------------------------
# If the re-extracted canonical labels do not reproduce the frozen ones, the
# sampling geometry is misaligned and NO metric below is interpretable.
ctl <- CMP_VARIANTS[[which(vapply(CMP_VARIANTS, function(v) v$tag, "") == CONTROL_TAG)]]
i_mk <- match(mk_cmp$pxid, lab$pxid)
control_df <- do.call(rbind, lapply(seq_along(CMP_YEARS), function(k) {
  yr  <- CMP_YEARS[k]
  old <- zero_blank(mk_cmp[[paste0("m", yr)]])
  new <- zero_blank(lab[[ctl$cols[k]]][i_mk])
  data.frame(year = yr, n = length(old), n_mismatch = sum(old != new),
             n_retama_masterkey = sum(old == 1), n_retama_reextracted = sum(new == 1),
             row.names = NULL)
}))
print(control_df)
if (any(control_df$n_mismatch > 0)) {
  stop("Control column does not reproduce the master key (",
       paste(control_df$n_mismatch, collapse = "/"), " mismatches). ",
       "Grid misalignment in 09v_valMapLabels.js — resolve before reading metrics.")
}
message("Control OK: re-extracted canonical labels reproduce the master key exactly.")

# ---- 3. RUN THE COMPUTE ENGINE ONCE PER VARIANT -----------------------------
# Fresh environment per run so nothing leaks between arms. Each run also writes
# its own full output set to analysis/<tag>/.
run_variant <- function(v) {
  e <- new.env(parent = globalenv())
  assign("VARIANT_TAG", v$tag, envir = e)
  assign("MAP_LABELS", data.frame(
    pxid  = lab$pxid,
    m2017 = zero_blank(lab[[v$cols[1]]]),
    m2025 = zero_blank(lab[[v$cols[2]]]),
    stringsAsFactors = FALSE), envir = e)
  message("--- computing variant: ", v$tag, " (", v$label, ")")
  sys.source(COMPUTE_R, envir = e)
  e
}
runs <- setNames(lapply(CMP_VARIANTS, run_variant),
                 vapply(CMP_VARIANTS, function(v) v$tag, ""))

var_label <- setNames(vapply(CMP_VARIANTS, function(v) v$label, ""),
                      vapply(CMP_VARIANTS, function(v) v$tag,   ""))
TAGS <- names(runs)

# The reference reconciliation is map-independent, so every arm must end up
# scoring exactly the same rows. If it doesn't, a join went wrong upstream.
key_of <- function(e) paste(e$usable$pxid, e$usable$year)
stopifnot(all(vapply(runs, function(e) identical(key_of(e), key_of(runs[[1]])), TRUE)))
message("Paired check OK: all arms scored on the same ",
        nrow(runs[[1]]$usable), " (pxid, year) rows.")

# ---- 4. SIDE-BY-SIDE METRIC TABLES ------------------------------------------
CRUDE_METRICS <- c("OA",
                   "UA_retama", "CE_retama", "PA_retama", "OE_retama",
                   "UA_noret",  "CE_noret",  "PA_noret",  "OE_noret",
                   "CE_noret_near", "CE_noret_far", "OE_noret_near", "OE_noret_far")

cmp_crude <- do.call(rbind, lapply(CMP_YEARS, function(yr) {
  do.call(rbind, lapply(CRUDE_METRICS, function(mm) {
    do.call(rbind, lapply(TAGS, function(tg) {
      ci <- runs[[tg]]$metrics_ci_df
      r  <- ci[ci$year == yr & ci$metric == mm, ]
      if (!nrow(r)) return(NULL)
      data.frame(year = yr, metric = mm, variant = tg, label = var_label[[tg]],
                 k = r$k, n = r$n, est = r$est, lo = r$lo, hi = r$hi,
                 row.names = NULL)
    }))
  }))
}))

OLO_METRICS <- c("OA", "UA_retama", "CE_retama", "PA_retama", "OE_retama",
                 "UA_noret", "CE_noret", "PA_noret", "OE_noret")

cmp_olof <- NULL
cmp_area <- NULL
if (!is.null(runs[[1]]$olofsson_df)) {
  cmp_olof <- do.call(rbind, lapply(CMP_YEARS, function(yr) {
    do.call(rbind, lapply(OLO_METRICS, function(mm) {
      do.call(rbind, lapply(TAGS, function(tg) {
        d <- runs[[tg]]$olofsson_df
        r <- d[d$year == yr, ]
        data.frame(year = yr, metric = mm, variant = tg, label = var_label[[tg]],
                   est = r[[mm]], se = r[[paste0(mm, "_se")]],
                   lo  = r[[paste0(mm, "_lo")]], hi = r[[paste0(mm, "_hi")]],
                   row.names = NULL)
      }))
    }))
  }))
  cmp_area <- do.call(rbind, lapply(CMP_YEARS, function(yr) {
    do.call(rbind, lapply(TAGS, function(tg) {
      d <- runs[[tg]]$olofsson_area_df
      r <- d[d$year == yr, ]
      data.frame(year = yr, variant = tg, label = var_label[[tg]],
                 area_prop = r$area_prop_retama, se = r$area_se,
                 lo = r$area_lo, hi = r$area_hi, row.names = NULL)
    }))
  }))
}

# ---- 5. PAIRED COMPARISON: baseline arm vs test arm -------------------------
# Same points, same references, only the map differs => McNemar on per-point
# correctness. Reported as UNWEIGHTED: it treats every sampled point equally,
# ignoring the stratification. The intervals that describe the ROI are the
# area-weighted ones in section 4.
paired <- runs[[BASELINE_TAG]]$usable[, c("pxid", "year", "stratum", "zone", "ref")]
paired$map_base <- runs[[BASELINE_TAG]]$usable$map
paired$map_test <- runs[[TEST_TAG]]$usable$map
paired$ok_base  <- paired$map_base == paired$ref
paired$ok_test  <- paired$map_test == paired$ref

mcnemar_df <- do.call(rbind, lapply(CMP_YEARS, function(yr) {
  p <- paired[paired$year == yr, ]
  b <- sum( p$ok_base & !p$ok_test)   # baseline right, test wrong
  c_<- sum(!p$ok_base &  p$ok_test)   # baseline wrong, test right
  tt <- tryCatch(
    stats::mcnemar.test(matrix(c(sum(p$ok_base & p$ok_test), b,
                                 c_, sum(!p$ok_base & !p$ok_test)), 2, 2)),
    error = function(e) NULL)
  data.frame(year = yr, n = nrow(p),
             n_both_correct = sum(p$ok_base & p$ok_test),
             n_both_wrong   = sum(!p$ok_base & !p$ok_test),
             n_only_base_correct = b, n_only_test_correct = c_,
             acc_base = mean(p$ok_base), acc_test = mean(p$ok_test),
             delta_acc = mean(p$ok_test) - mean(p$ok_base),
             p_value = if (is.null(tt)) NA_real_ else unname(tt$p.value),
             row.names = NULL)
}))

# Which points changed label, in which direction, and was the change right?
flow_df <- do.call(rbind, lapply(CMP_YEARS, function(yr) {
  p <- paired[paired$year == yr & paired$map_base != paired$map_test, ]
  if (!nrow(p)) return(data.frame(year = yr, direction = character(),
                                  stratum = integer(), n = integer(),
                                  n_test_correct = integer(), row.names = NULL))
  p$direction <- ifelse(p$map_base == 0 & p$map_test == 1,
                        "0->1 (test adds retama)", "1->0 (test removes retama)")
  ag <- aggregate(data.frame(n = rep(1L, nrow(p)), n_test_correct = as.integer(p$ok_test)),
                  by = list(direction = p$direction, stratum = p$stratum), FUN = sum)
  data.frame(year = yr, ag[order(ag$direction, ag$stratum), ], row.names = NULL)
}))

# ---- 6. EXPORT --------------------------------------------------------------
w <- function(df, f) if (!is.null(df)) write.csv(df, file.path(OUT_DIR, f), row.names = FALSE)
w(control_df,  "09v_cmp_control.csv")
w(cmp_crude,   "09v_cmp_metrics_crude.csv")
w(cmp_olof,    "09v_cmp_metrics_olofsson.csv")
w(cmp_area,    "09v_cmp_area.csv")
w(mcnemar_df,  "09v_cmp_mcnemar.csv")
w(flow_df,     "09v_cmp_flow.csv")

cat("\ncompare done. Outputs in:", OUT_DIR, "\n")
print(mcnemar_df)
