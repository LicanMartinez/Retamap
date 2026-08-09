# =============================================================================
# 09v_valAnalysis_compute.R
# -----------------------------------------------------------------------------
# Compute engine for the independent-validation report (family 09v_).
# Consumes:
#   - validation/09v_master_key.csv                      (map labels + stratum +
#                                                        authoritative validators)
#   - validation/sheets_filled/<val>_filled.csv          (reference labels)
# Produces (in validation/analysis/):
#   - 09v_confusion_<year>.csv       overall 2x2 confusion (map vs reference)
#   - 09v_perStratum_<year>.csv      per-stratum (S1/S2/S3) confusion cells
#   - 09v_metrics_unweighted.csv     OA / commission / omission (sample-based)
#   - 09v_metrics_olofsson.csv       area-weighted (only if Wi provided)
#   - 09v_kappa.csv                  inter-observer agreement on paired points
#   - 09v_pairConfidence_scatter.png/.pdf   the two-observer confidence figure
#   - 09v_metrics_digest.md          human-readable dump of ALL numbers
#
# This script is standalone (base R + ggplot2, no pandoc) AND is sourced by
# 09v_valAnalysis.Rmd. Running it alone == "pass 1" (export tables/figures/digest);
# the Rmd then embeds the same objects and adds narrative == the final report.
#
# Reference reconciliation per (pxid, year), driven by master-key `validators`
# (source of truth after the 2026-07-27 rebalancing):
#   - single validator  -> that label
#   - pair (2 validators) -> consensus if they agree; DISCARDED that year if they
#                            disagree (reported as %)
#   - confidence == 0 (no determinable) -> that label dropped for that year
# Map label per (pxid, year) = m2017 / m2025 (blank = 0 no-retama, 1 = retama),
# from the master key by default or from MAP_LABELS when evaluating a variant.
# Strata: 1 = S1 mapped-retama, 2 = S2 near background, 3 = S3 far background.
# near/far breakdown of the no-retama class: near = S1 + S2 (both inside the
# ~1 km near-retama envelope by construction); far = S3.
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

# ---- 0. CONFIG --------------------------------------------------------------
PROJECT_DIR <- "D:/Lican/uni/Investigacion/Colaboraciones_y_ayudas/Retamap"
WORK_DIR    <- file.path(PROJECT_DIR, "gee_scripts", "validation")
SHEETS_DIR  <- file.path(WORK_DIR, "sheets_filled")
MASTER_KEY  <- file.path(WORK_DIR, "09v_master_key.csv")
YEARS       <- c(2017, 2025)
VALIDATORS  <- c("lican", "sofi", "jaime")
SEED        <- 42

# --- MAP VARIANT hooks (defaults reproduce the canonical run exactly) --------
# Set these BEFORE sourcing this file to evaluate an ALTERNATIVE map on the same
# points and the same human reference labels (see 09v_valCompare.R):
#   VARIANT_TAG : subfolder of analysis/ for this run's outputs ("" = canonical)
#   MAP_LABELS  : data.frame(pxid, m2017, m2025) replacing the master key's map
#                 columns. Everything downstream derives from those two columns,
#                 so this single override is the whole variant mechanism.
# The reference labels, strata and Wi are properties of the sampling DESIGN and
# stay frozen — only the map under evaluation changes.
# inherits = FALSE: only an explicit assignment in the environment this script is
# evaluated in counts, so a stray global of the same name can never leak into a
# run (09v_valCompare.R sources this repeatedly into fresh environments).
if (!exists("VARIANT_TAG", inherits = FALSE)) VARIANT_TAG <- ""
if (!exists("MAP_LABELS",  inherits = FALSE)) MAP_LABELS  <- NULL

OUT_DIR <- if (nzchar(VARIANT_TAG)) {
  file.path(WORK_DIR, "analysis", VARIANT_TAG)
} else {
  file.path(WORK_DIR, "analysis")
}

# Stratum AREA WEIGHTS Wi (proportion of ROI area in each stratum), printed by
# 09v_valGenerator.js section 5 ("Wi ... suma~1") in the GEE console.
# Fill all three to enable the Olofsson area-weighted estimator.
# Leave any as NA -> only the crude (sample-based) metrics are produced.
# Loaded 2026-08-07 from the GEE console. Stratum areas (m^2) behind them:
#   S1 =    55,178,037   S2 =   591,685,230   S3 = 9,273,053,977   (sum ~ 9.92e9)
# S3 is 93.5% of the ROI, so the area-weighted estimator is dominated by the far
# background — which is exactly why the crude metrics overstate the error.
Wi <- c(S1 = 0.005562348528027356,
        S2 = 0.059646186135871575,
        S3 = 0.934791465336101)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. LOAD & VALIDATE -----------------------------------------------------
mk <- read.csv(MASTER_KEY, stringsAsFactors = FALSE, colClasses = "character")
stopifnot(all(c("pxid", "stratum", "m2017", "m2025", "validators") %in% names(mk)))
mk$stratum <- as.integer(mk$stratum)
zero_blank <- function(x) ifelse(is.na(x) | x == "", 0L, as.integer(x))
mk$m2017 <- zero_blank(mk$m2017)
mk$m2025 <- zero_blank(mk$m2025)

# Alternative map: swap in the variant's labels, keyed by pxid. Must cover every
# pxid — a partial join would silently shrink the evaluated sample and make the
# variant look different for the wrong reason.
if (!is.null(MAP_LABELS)) {
  stopifnot(all(c("pxid", "m2017", "m2025") %in% names(MAP_LABELS)))
  stopifnot(!anyDuplicated(MAP_LABELS$pxid))
  miss <- setdiff(mk$pxid, MAP_LABELS$pxid)
  if (length(miss)) stop("MAP_LABELS missing ", length(miss), " pxid, e.g.: ",
                         paste(utils::head(miss, 5), collapse = ", "))
  i <- match(mk$pxid, MAP_LABELS$pxid)
  mk$m2017 <- zero_blank(MAP_LABELS$m2017[i])
  mk$m2025 <- zero_blank(MAP_LABELS$m2025[i])
}

mk_val   <- setNames(strsplit(mk$validators, ";"), mk$pxid)      # pxid -> chr vec
mk_strat <- setNames(mk$stratum, mk$pxid)
mapcol   <- list("2017" = setNames(mk$m2017, mk$pxid),
                 "2025" = setNames(mk$m2025, mk$pxid))

sheets <- do.call(rbind, lapply(VALIDATORS, function(v) {
  f <- file.path(SHEETS_DIR, paste0(v, "_filled.csv"))
  if (!file.exists(f)) { warning("missing sheet: ", f); return(NULL) }
  read.csv(f, stringsAsFactors = FALSE, colClasses = "character")
}))
stopifnot(all(c("validator", "pxid", "class_2017", "conf_2017",
                "class_2025", "conf_2025") %in% names(sheets)))

# domain checks (loud failure on any transcription slip)
bad_pxid <- setdiff(sheets$pxid, mk$pxid)
if (length(bad_pxid)) stop("pxid not in master key: ", paste(bad_pxid, collapse = ", "))
ok_class <- function(x) all(x %in% c("Retama", "No_retama", "", NA))
stopifnot(ok_class(sheets$class_2017), ok_class(sheets$class_2025))
sheets$conf_2017 <- suppressWarnings(as.integer(sheets$conf_2017))
sheets$conf_2025 <- suppressWarnings(as.integer(sheets$conf_2025))
chk_conf <- function(x) all(is.na(x) | x %in% 0:10)
stopifnot(chk_conf(sheets$conf_2017), chk_conf(sheets$conf_2025))

# ---- 2. LONG TABLE (pxid, validator, year, class01, conf) -------------------
to_long <- function(df, yr) {
  cl <- df[[paste0("class_", yr)]]
  cf <- df[[paste0("conf_",  yr)]]
  data.frame(pxid = df$pxid, validator = df$validator, year = yr,
             class01 = ifelse(cl == "Retama", 1L,
                       ifelse(cl == "No_retama", 0L, NA_integer_)),
             conf = cf, stringsAsFactors = FALSE)
}
long <- rbind(to_long(sheets, 2017), to_long(sheets, 2025))
long <- long[!is.na(long$class01), ]                 # keep only filled-class rows

# ---- 3. REFERENCE RECONCILIATION per (pxid, year) ---------------------------
# Uses ONLY the validators listed in the master key for that pxid.
ana_rows <- list()
for (px in mk$pxid) {
  assigned <- mk_val[[px]]
  for (yr in YEARS) {
    sub <- long[long$pxid == px & long$year == yr &
                long$validator %in% assigned, ]
    sub <- sub[!is.na(sub$conf) & sub$conf > 0, ]     # drop conf==0 (no determinable)
    n_obs <- nrow(sub)
    if (n_obs == 0) next
    labs <- unique(sub$class01)
    status <- if (length(labs) > 1) "disagree"
              else if (n_obs >= 2)  "agree" else "single"
    ref <- if (length(labs) == 1) labs else NA_integer_
    ana_rows[[length(ana_rows) + 1]] <- data.frame(
      pxid = px, year = yr, stratum = mk_strat[[px]],
      map = mapcol[[as.character(yr)]][[px]],
      ref = ref, status = status, n_obs = n_obs,
      n_assigned = length(assigned),
      # confidence of the lone observer (NA for 2-observer points): feeds the
      # single-validator marginal band of the confidence-lattice figures
      conf_single = if (n_obs == 1) sub$conf[1] else NA_integer_,
      stringsAsFactors = FALSE)
  }
}
ana <- do.call(rbind, ana_rows)

# usable rows for the confusion matrix: a defined reference (agree|single)
usable <- ana[ana$status %in% c("agree", "single") & !is.na(ana$ref), ]

# Zone for the no-retama near/far split. Derived from the STRATUM alone: S1 and
# S2 are both inside the ~1 km near-retama envelope by construction, S3 is the
# far background. Deliberately independent of the map being evaluated, so an
# alternative map (09v_valCompare.R) is split the same way as the canonical one.
# Identical to the previous stratum+map rule for every reported quantity: CE is
# computed only among map == 0 points, where S1 on-year (map == 1) drops out on
# its own, and OE is now counted per zone instead of being assumed.
usable$zone <- ifelse(usable$stratum == 3, "far", "near")

# ---- 4. CRUDE (sample-based) METRICS per year -------------------------------
conf2 <- function(map, ref) c(
  TP = sum(map == 1 & ref == 1), FP = sum(map == 1 & ref == 0),
  FN = sum(map == 0 & ref == 1), TN = sum(map == 0 & ref == 0))
safe <- function(a, b) if (b == 0) NA_real_ else unname(a / b)
metrics_from_conf <- function(cc) {
  TP <- cc["TP"]; FP <- cc["FP"]; FN <- cc["FN"]; TN <- cc["TN"]; N <- sum(cc)
  c(N = N, OA = safe(TP + TN, N),
    UA_retama = safe(TP, TP + FP), PA_retama = safe(TP, TP + FN),
    CE_retama = safe(FP, TP + FP), OE_retama = safe(FN, TP + FN),
    UA_noret  = safe(TN, TN + FN), PA_noret  = safe(TN, TN + FP),
    CE_noret  = safe(FN, FN + TN), OE_noret  = safe(FP, FP + TN))
}

# --- binomial CIs for the crude (sample-based) ratios ------------------------
# Every crude metric is "k successes out of n trials", so a Wilson score interval
# is the right tool here: unlike the normal approximation it stays inside [0,1]
# and still returns a usable upper bound when k == 0 (e.g. CE_far = 0/n, where a
# Wald interval would collapse to the uninformative [0, 0]).
# NOTE these are UNWEIGHTED intervals: they describe the precision of the sample
# as drawn, ignoring the stratification. The design-relevant intervals are the
# area-weighted ones in section 5 (Olofsson) — see section 5b for the comparison
# against the SE the sampling design was sized for.
Z95 <- qnorm(0.975)
wilson_ci <- function(k, n, z = Z95) {
  if (is.na(k) || is.na(n) || n <= 0)
    return(c(est = NA_real_, se = NA_real_, lo = NA_real_, hi = NA_real_))
  p   <- k / n
  den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den
  hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(est = p, se = sqrt(p * (1 - p) / n),
    lo = max(0, ctr - hw), hi = min(1, ctr + hw))
}
# the (k, n) pair behind each crude metric, given a confusion vector cc
crude_kn <- function(cc) {
  TP <- unname(cc["TP"]); FP <- unname(cc["FP"])
  FN <- unname(cc["FN"]); TN <- unname(cc["TN"])
  list(OA        = c(TP + TN, TP + FP + FN + TN),
       UA_retama = c(TP, TP + FP), CE_retama = c(FP, TP + FP),
       PA_retama = c(TP, TP + FN), OE_retama = c(FN, TP + FN),
       UA_noret  = c(TN, TN + FN), CE_noret  = c(FN, TN + FN),
       PA_noret  = c(TN, TN + FP), OE_noret  = c(FP, TN + FP))
}
ci_table <- function(kn, yr) do.call(rbind, lapply(names(kn), function(mm) {
  w <- wilson_ci(kn[[mm]][1], kn[[mm]][2])
  data.frame(year = yr, metric = mm, k = kn[[mm]][1], n = kn[[mm]][2],
             est = unname(w["est"]), se = unname(w["se"]),
             lo = unname(w["lo"]), hi = unname(w["hi"]),
             method = "Wilson 95%", row.names = NULL)
}))

conf_by_year   <- list()
metrics_uw     <- list()
metrics_ci     <- list()
perstratum     <- list()
nearfar        <- list()
qc_by_year     <- list()
for (yr in YEARS) {
  uy <- usable[usable$year == yr, ]
  cc <- conf2(uy$map, uy$ref)
  conf_by_year[[as.character(yr)]] <- cc
  m  <- metrics_from_conf(cc)
  # near/far breakdown of the no-retama class, counted per zone.
  # CE_noret = FN / (FN + TN) among map-no-retama points of the zone.
  # OE_noret = FP / (FP + TN) among reference-no-retama points of the zone.
  # FP is counted, not assumed: for the canonical map every FP necessarily sits
  # in S1 (S1 IS "map == 1 in 2017 and/or 2025", so no FP can fall in S2/S3) and
  # far therefore comes out at exactly 0 — but that identity is a property of
  # the strata-defining map only. An alternative map can put retama on S2/S3
  # points, and hard-coding OE_far = 0 would misattribute those errors to near.
  zcell <- function(z) {
    s  <- uy[uy$zone == z, ]
    FN <- sum(s$map == 0 & s$ref == 1); TN <- sum(s$map == 0 & s$ref == 0)
    FP <- sum(s$map == 1 & s$ref == 0)
    c(FN = FN, TN = TN, FP = FP, n = FN + TN,
      CE_noret = safe(FN, FN + TN), OE_noret = safe(FP, FP + TN))
  }
  near <- zcell("near"); far <- zcell("far")
  oe_near <- unname(near["OE_noret"]); oe_far <- unname(far["OE_noret"])
  # Wilson CIs for every crude ratio, including the near/far breakdown
  kn <- c(crude_kn(cc), list(
    CE_noret_near = c(unname(near["FN"]), unname(near["FN"] + near["TN"])),
    CE_noret_far  = c(unname(far["FN"]),  unname(far["FN"]  + far["TN"])),
    OE_noret_near = c(unname(near["FP"]), unname(near["FP"] + near["TN"])),
    OE_noret_far  = c(unname(far["FP"]),  unname(far["FP"]  + far["TN"]))))
  cit <- ci_table(kn, yr)
  metrics_ci[[as.character(yr)]] <- cit
  # same numbers in wide form (<metric>_lo / <metric>_hi) for inline use in the Rmd
  wide <- as.data.frame(setNames(
    as.list(as.vector(t(as.matrix(cit[, c("lo", "hi")])))),
    as.vector(t(outer(cit$metric, c("_lo", "_hi"), paste0)))))
  metrics_uw[[as.character(yr)]] <- data.frame(year = yr, t(m),
    CE_noret_near = unname(near["CE_noret"]), CE_noret_far = unname(far["CE_noret"]),
    OE_noret_near = unname(oe_near),          OE_noret_far = unname(oe_far),
    n_near = unname(near["n"]), n_far = unname(far["n"]), wide, row.names = NULL)
  nearfar[[as.character(yr)]] <- rbind(
    data.frame(year = yr, zone = "near", t(near)),
    data.frame(year = yr, zone = "far",  t(far)))
  # per-stratum confusion
  ps <- do.call(rbind, lapply(sort(unique(uy$stratum)), function(st) {
    s <- uy[uy$stratum == st, ]; c2 <- conf2(s$map, s$ref)
    data.frame(year = yr, stratum = st, n = nrow(s), t(c2),
               OA = safe(c2["TP"] + c2["TN"], sum(c2)), row.names = NULL)
  }))
  perstratum[[as.character(yr)]] <- ps
  # QC counts for this year
  ay <- ana[ana$year == yr, ]
  qc_by_year[[as.character(yr)]] <- data.frame(year = yr,
    n_usable = nrow(uy),
    n_disagree = sum(ay$status == "disagree"),
    n_pairs_eval = sum(ay$n_assigned == 2),
    # single-validator points feeding `usable`, split by WHY they're single:
    # by design (only 1 validator was ever assigned, e.g. most of S2/S3) vs.
    # a still-incomplete pair (assigned 2, only 1 has loaded data so far)
    n_single_design  = sum(ay$status == "single" & ay$n_assigned == 1),
    n_single_partial = sum(ay$status == "single" & ay$n_assigned == 2),
    n_conf0_excluded = NA_integer_)  # filled below from raw long
}
metrics_uw_df <- do.call(rbind, metrics_uw)
metrics_ci_df <- do.call(rbind, metrics_ci)
nearfar_df    <- do.call(rbind, nearfar)
qc_df         <- do.call(rbind, qc_by_year)

# conf==0 exclusions per year (labels present but dropped)
for (i in seq_len(nrow(qc_df))) {
  yr <- qc_df$year[i]
  qc_df$n_conf0_excluded[i] <- sum(long$year == yr & !is.na(long$conf) & long$conf == 0)
}

# ---- 5. OLOFSSON area-weighted estimator (only if Wi supplied) --------------
olofsson_df <- NULL
olofsson_area_df <- NULL
if (all(!is.na(Wi))) {
  z  <- qnorm(0.975)
  Wv <- c("1" = Wi[["S1"]], "2" = Wi[["S2"]], "3" = Wi[["S3"]])
  ol <- list(); ar <- list()
  for (yr in YEARS) {
    uy <- usable[usable$year == yr, ]
    strata <- c("1", "2", "3")
    # per-stratum cell fractions and n_h
    p <- matrix(0, 2, 2, dimnames = list(map = c("1", "0"), ref = c("1", "0")))
    nh <- setNames(integer(length(strata)), strata)
    ph_correct <- setNames(numeric(length(strata)), strata)   # y_h = fraction correct
    ph_refret  <- setNames(numeric(length(strata)), strata)   # fraction ref==retama
    # store per-stratum fractions for ratio variances
    fr <- list()
    for (st in strata) {
      s <- uy[uy$stratum == as.integer(st), ]
      nh[st] <- nrow(s)
      if (nrow(s) == 0) next
      cc <- conf2(s$map, s$ref) / nrow(s)   # fractions n_hij/n_h
      W  <- Wv[[st]]
      p["1", "1"] <- p["1", "1"] + W * cc["TP"]
      p["1", "0"] <- p["1", "0"] + W * cc["FP"]
      p["0", "1"] <- p["0", "1"] + W * cc["FN"]
      p["0", "0"] <- p["0", "0"] + W * cc["TN"]
      ph_correct[st] <- (cc["TP"] + cc["TN"])
      ph_refret[st]  <- (cc["TP"] + cc["FN"])
      fr[[st]] <- cc
    }
    OA <- p["1","1"] + p["0","0"]
    UA_ret <- safe(p["1","1"], p["1","1"] + p["1","0"])
    PA_ret <- safe(p["1","1"], p["1","1"] + p["0","1"])
    UA_nor <- safe(p["0","0"], p["0","0"] + p["0","1"])
    PA_nor <- safe(p["0","0"], p["0","0"] + p["1","0"])
    # OA variance (Olofsson eq 5)
    vOA <- sum(sapply(strata, function(st) if (nh[st] > 1)
      Wv[[st]]^2 * ph_correct[st] * (1 - ph_correct[st]) / (nh[st] - 1) else 0))
    # ratio variance (delta method) for a user/producer accuracy R = X/Y
    ratio_var <- function(xf, yf) {   # xf,yf: named-by-stratum fractions (numerator/denominator indicators)
      X <- sum(sapply(strata, function(st) Wv[[st]] * xf[[st]]))
      Y <- sum(sapply(strata, function(st) Wv[[st]] * yf[[st]]))
      if (is.na(Y) || Y == 0) return(c(R = NA_real_, V = NA_real_))
      R <- X / Y
      V <- 0
      for (st in strata) {
        # nh[[st]] (not nh[st]): single-bracket keeps the stratum name, which
        # then rides along into V and turns c(R=, V=) into c(R=, V.3=) — the
        # same name-propagation trap as the safe() bug fixed on 2026-07-27.
        nhs <- nh[[st]]
        if (nhs < 2) next
        xh <- xf[[st]]; yh <- yf[[st]]
        vx <- xh * (1 - xh) / (nhs - 1)
        vy <- yh * (1 - yh) / (nhs - 1)
        cxy <- (xh - xh * yh) / (nhs - 1)      # x subset of y => E[xy]=E[x]
        V <- V + Wv[[st]]^2 * (vx + R^2 * vy - 2 * R * cxy)
      }
      c(R = unname(R), V = unname(V / Y^2))
    }
    getf <- function(st, cell) if (is.null(fr[[st]])) 0 else unname(fr[[st]][cell])
    # numerator / denominator indicator fractions for each of the four ratios
    sumf <- function(...) { cells <- c(...); setNames(lapply(strata, function(st)
      sum(vapply(cells, function(cl) getf(st, cl), 0))), strata) }
    v_UAret <- ratio_var(sumf("TP"), sumf("TP", "FP"))
    v_PAret <- ratio_var(sumf("TP"), sumf("TP", "FN"))
    v_UAnor <- ratio_var(sumf("TN"), sumf("TN", "FN"))
    v_PAnor <- ratio_var(sumf("TN"), sumf("TN", "FP"))
    seOA <- sqrt(max(vOA, 0))
    # SE / CI of a ratio; the CE/OE complement (1 - R) keeps the same SE and has
    # its bounds mirrored, so it is derived rather than re-estimated.
    se_of <- function(v) if (is.na(v[["V"]])) NA_real_ else sqrt(max(v[["V"]], 0))
    ci    <- function(est, v) { s <- se_of(v)
      if (is.na(s) || is.na(est)) c(NA_real_, NA_real_) else est + c(-1, 1) * z * s }
    cUAret <- ci(UA_ret, v_UAret); cPAret <- ci(PA_ret, v_PAret)
    cUAnor <- ci(UA_nor, v_UAnor); cPAnor <- ci(PA_nor, v_PAnor)
    ol[[as.character(yr)]] <- data.frame(
      year = yr,
      OA = OA, OA_se = seOA, OA_lo = OA - z*seOA, OA_hi = OA + z*seOA,
      UA_retama = UA_ret, UA_retama_se = se_of(v_UAret), UA_retama_lo = cUAret[1], UA_retama_hi = cUAret[2],
      PA_retama = PA_ret, PA_retama_se = se_of(v_PAret), PA_retama_lo = cPAret[1], PA_retama_hi = cPAret[2],
      UA_noret  = UA_nor, UA_noret_se  = se_of(v_UAnor), UA_noret_lo  = cUAnor[1], UA_noret_hi  = cUAnor[2],
      PA_noret  = PA_nor, PA_noret_se  = se_of(v_PAnor), PA_noret_lo  = cPAnor[1], PA_noret_hi  = cPAnor[2],
      # commission / omission = 1 - UA / 1 - PA  (same SE, mirrored bounds)
      CE_retama = 1 - UA_ret, CE_retama_se = se_of(v_UAret), CE_retama_lo = 1 - cUAret[2], CE_retama_hi = 1 - cUAret[1],
      OE_retama = 1 - PA_ret, OE_retama_se = se_of(v_PAret), OE_retama_lo = 1 - cPAret[2], OE_retama_hi = 1 - cPAret[1],
      CE_noret  = 1 - UA_nor, CE_noret_se  = se_of(v_UAnor), CE_noret_lo  = 1 - cUAnor[2], CE_noret_hi  = 1 - cUAnor[1],
      OE_noret  = 1 - PA_nor, OE_noret_se  = se_of(v_PAnor), OE_noret_lo  = 1 - cPAnor[2], OE_noret_hi  = 1 - cPAnor[1],
      row.names = NULL)
    # area of retama (proportion) + CI (stratified proportion estimator)
    Aret <- sum(sapply(strata, function(st) Wv[[st]] * ph_refret[st]))
    vA <- sum(sapply(strata, function(st) if (nh[st] > 1)
      Wv[[st]]^2 * ph_refret[st] * (1 - ph_refret[st]) / (nh[st] - 1) else 0))
    ar[[as.character(yr)]] <- data.frame(year = yr,
      area_prop_retama = Aret, area_se = sqrt(max(vA, 0)),
      area_lo = Aret - z*sqrt(vA), area_hi = Aret + z*sqrt(vA),
      row.names = NULL)
  }
  olofsson_df      <- do.call(rbind, ol)
  olofsson_area_df <- do.call(rbind, ar)
}

# ---- 5b. DESIGN CHECK: achieved precision vs. the sample-size design ---------
# 09v_valGenerator.js sized each stratum with the Olofsson rule
#     n_i = U_i (1 - U_i) / SE_i^2      (floored at FLOOR_i = 100)
# using EXP_UA = {.80, .80, .95} and TARGET_SE = {.02, .02, .02}, which is why
# 400 / 400 / 119 points were drawn. This block asks whether the data actually
# delivers that precision, separating the two ways it can fail:
#   (a) LABELS INCOMPLETE — fewer points evaluated than drawn, so n_h < n_drawn.
#       `se_at_n_drawn` is the precision once every assigned label is loaded.
#   (b) ACCURACY DIFFERS FROM THE ASSUMPTION — if the realized p_h is further
#       from 0 or 1 than EXP_UA, p(1-p) is larger and the target SE needs more
#       points than the design bought, however complete the loading is.
# p_h = the stratum's proportion of correctly classified points (map == ref):
# the same quantity that enters the OA variance, and the realized analogue of
# the U_i the design assumed. Needs no Wi, so it runs even without area weights.
DESIGN <- data.frame(
  stratum   = 1:3,
  label     = c("S1 mapped-retama", "S2 near background", "S3 far background"),
  exp_UA    = c(0.80, 0.80, 0.95),   # EXP_UA    in 09v_valGenerator.js sec. 0
  target_SE = c(0.02, 0.02, 0.02),   # TARGET_SE in 09v_valGenerator.js sec. 0
  floor_n   = c(100L, 100L, 100L),   # FLOOR     in 09v_valGenerator.js sec. 0
  n_drawn   = as.integer(table(factor(mk$stratum, levels = 1:3))),
  stringsAsFactors = FALSE)
DESIGN$n_formula <- ceiling(DESIGN$exp_UA * (1 - DESIGN$exp_UA) / DESIGN$target_SE^2)
DESIGN$n_design  <- pmax(DESIGN$n_formula, DESIGN$floor_n)

design_rows <- list()
for (yr in YEARS) {
  uy <- usable[usable$year == yr, ]
  for (i in seq_len(nrow(DESIGN))) {
    st <- DESIGN$stratum[i]
    s  <- uy[uy$stratum == st, ]
    nl <- nrow(s)
    p  <- if (nl > 0) mean(s$map == s$ref) else NA_real_
    se_now   <- if (nl > 0) sqrt(p * (1 - p) / nl) else NA_real_
    se_full  <- if (nl > 0) sqrt(p * (1 - p) / DESIGN$n_drawn[i]) else NA_real_
    # points the realized p_h would require to hit the target SE
    n_req    <- if (nl > 0) ceiling(p * (1 - p) / DESIGN$target_SE[i]^2) else NA_real_
    design_rows[[length(design_rows) + 1]] <- data.frame(
      year = yr, stratum = st, label = DESIGN$label[i],
      exp_UA = DESIGN$exp_UA[i], obs_p_correct = p,
      target_SE = DESIGN$target_SE[i],
      n_design = DESIGN$n_design[i], n_drawn = DESIGN$n_drawn[i], n_loaded = nl,
      pct_loaded = if (DESIGN$n_drawn[i] > 0) nl / DESIGN$n_drawn[i] else NA_real_,
      se_achieved = se_now, se_at_n_drawn = se_full,
      se_ratio_vs_target = if (is.na(se_now)) NA_real_ else se_now / DESIGN$target_SE[i],
      n_required_at_obs_p = n_req,
      meets_target_now  = !is.na(se_now)  && se_now  <= DESIGN$target_SE[i],
      meets_target_full = !is.na(se_full) && se_full <= DESIGN$target_SE[i],
      row.names = NULL)
  }
}
design_df <- do.call(rbind, design_rows)

# OA-level design check (needs Wi: SE(OA) is an area-weighted combination).
# se_OA_design = the value 09v_valGenerator.js prints as "SE(OA) achieved given
# nFinal" — the design's own promise, computed from EXP_UA at the drawn n.
design_oa_df <- NULL
if (all(!is.na(Wi))) {
  Wv <- c(Wi[["S1"]], Wi[["S2"]], Wi[["S3"]])
  se_oa_design <- sqrt(sum(Wv^2 * DESIGN$exp_UA * (1 - DESIGN$exp_UA) / DESIGN$n_drawn))
  design_oa_df <- do.call(rbind, lapply(YEARS, function(yr) {
    d  <- design_df[design_df$year == yr, ]
    ob <- olofsson_df[olofsson_df$year == yr, ]
    # SE(OA) the realized p_h would give at full loading
    se_full <- sqrt(sum(Wv^2 * d$obs_p_correct * (1 - d$obs_p_correct) / d$n_drawn,
                        na.rm = TRUE))
    data.frame(year = yr,
      se_OA_design = se_oa_design,          # promised by the design
      se_OA_achieved = ob$OA_se,            # from the data as loaded now
      se_OA_at_full_loading = se_full,      # projected once all labels are in
      ratio_achieved_vs_design = ob$OA_se / se_oa_design,
      row.names = NULL)
  }))
}

# ---- 6. INTER-OBSERVER (paired points only) ---------------------------------
pair_pxids <- names(mk_val)[vapply(mk_val, length, 1L) == 2]
kappa_rows <- list()
pair_scatter <- list()
for (yr in YEARS) {
  a_lab <- c(); b_lab <- c()
  for (px in pair_pxids) {
    vs <- sort(mk_val[[px]])                    # alphabetical: A=vs[1], B=vs[2]
    s  <- long[long$pxid == px & long$year == yr & long$validator %in% vs &
               !is.na(long$conf) & long$conf > 0, ]
    if (length(unique(s$validator)) < 2) next   # need both observers this year
    A <- s[s$validator == vs[1], ][1, ]
    B <- s[s$validator == vs[2], ][1, ]
    a_lab <- c(a_lab, A$class01); b_lab <- c(b_lab, B$class01)
    vote <- if (A$class01 == 1 && B$class01 == 1) "both_retama"
            else if (A$class01 == 0 && B$class01 == 0) "both_noret" else "split"
    pair_scatter[[length(pair_scatter) + 1]] <- data.frame(
      pxid = px, year = yr, vA = vs[1], vB = vs[2],
      x = A$conf, y = B$conf, vote = vote,
      map = mapcol[[as.character(yr)]][[px]], stringsAsFactors = FALSE)
  }
  n <- length(a_lab)
  po <- if (n) mean(a_lab == b_lab) else NA
  pe <- if (n) {
    pa1 <- mean(a_lab == 1); pb1 <- mean(b_lab == 1)
    pa1 * pb1 + (1 - pa1) * (1 - pb1)
  } else NA
  kappa <- if (!is.na(pe) && pe < 1) (po - pe) / (1 - pe) else NA
  kappa_rows[[as.character(yr)]] <- data.frame(year = yr, n_pairs = n,
    pct_agreement = po, kappa = kappa,
    both_retama = sum(a_lab == 1 & b_lab == 1),
    both_noret  = sum(a_lab == 0 & b_lab == 0),
    split       = sum(a_lab != b_lab))
}
kappa_df <- do.call(rbind, kappa_rows)
pair_df  <- if (length(pair_scatter)) do.call(rbind, pair_scatter) else
            data.frame(pxid=character(), year=integer(), vA=character(), vB=character(),
                       x=integer(), y=integer(), vote=character(), map=integer())

# ---- 7. PAIRED-CONFIDENCE SCATTER (shared jitter for point + segment) -------
set.seed(SEED)
pair_df$jx <- pair_df$x + runif(nrow(pair_df), -0.18, 0.18)
pair_df$jy <- pair_df$y + runif(nrow(pair_df), -0.18, 0.18)
pair_df$vote     <- factor(pair_df$vote, levels = c("both_retama", "split", "both_noret"))
pair_df$map_lab  <- factor(ifelse(pair_df$map == 1, "map: Retama", "map: No_retama"),
                           levels = c("map: Retama", "map: No_retama"))
pair_df$year_f   <- factor(pair_df$year, levels = YEARS)

# segments connecting 2017 -> 2025 of the same pxid (uses the jittered coords)
w17 <- pair_df[pair_df$year == YEARS[1], c("pxid", "jx", "jy")]
w25 <- pair_df[pair_df$year == YEARS[2], c("pxid", "jx", "jy")]
seg <- merge(w17, w25, by = "pxid", suffixes = c("17", "25"))

fill_cols <- c(both_retama = "#E69F00", split = "#F0E442", both_noret = "#009E73")
brd_cols  <- c("map: Retama" = "#B8860B", "map: No_retama" = "#006400")

p_scatter <- ggplot() +
  { if (nrow(seg))
      geom_segment(data = seg,
        aes(x = jx17, y = jy17, xend = jx25, yend = jy25),
        color = "grey70", alpha = 0.5, linewidth = 0.4,
        arrow = arrow(length = unit(0.16, "cm"), type = "closed")) } +
  geom_point(data = pair_df,
    aes(x = jx, y = jy, fill = vote, color = map_lab),
    shape = 21, size = 3.6, stroke = 1.4, alpha = 0.95) +
  scale_fill_manual(values = fill_cols,
    labels = c(both_retama = "ambos: Retama", split = "1 vs 1 (split)",
               both_noret = "ambos: No_retama"),
    name = "Referencia:\nvotos de los\ndos evaluadores") +
  scale_color_manual(values = brd_cols, name = "Mapa\n(lo evaluado)") +
  scale_x_continuous(limits = c(-0.3, 10.3), breaks = 0:10) +
  scale_y_continuous(limits = c(-0.3, 10.3), breaks = 0:10) +
  labs(
    title = "Confianza de los dos evaluadores en puntos con doble validacion",
    subtitle = paste0("Cada punto = un (pxid, anio); la flecha conecta ", YEARS[1],
                      " -> ", YEARS[2], " del mismo pixel."),
    caption = "Ejes: confianza del evaluador A / B por orden alfabetico dentro de cada par.",
    x = "Confianza evaluador A (0-10)", y = "Confianza evaluador B (0-10)") +
  coord_equal() + theme_bw(base_size = 12) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "09v_pairConfidence_scatter.png"), p_scatter,
       width = 9, height = 6.5, dpi = 150)
ggsave(file.path(OUT_DIR, "09v_pairConfidence_scatter.pdf"), p_scatter,
       width = 9, height = 6.5)

# ---- 7b. GRID VERSIONS OF THE PAIRED-CONFIDENCE FIGURE (per year) -----------
# Terminology (applies to the whole script, not just from here): the REFERENCE
# ("valor real") is always the validators' photo-interpretation; the MAP is
# always the product being evaluated. Every metric upstream (conf2/metrics_from_conf,
# section 4-5) already computes error this way (FP/FN defined against `ref`); this
# section just needs to keep repeating it in labels since it's visual, not tabular.
#
# Both figures live on the same 10x10 discrete lattice (confidence of evaluator
# A x confidence of evaluator B) PLUS a marginal band for single-validator points,
# which have only one confidence value and so can't sit on the 2-D lattice:
#   Figure 1 (bubble): grid = % of that year's PAIRED points on each (A,B)
#                      combination; band ("1 val.", y=Y_BAND) = % of that year's
#                      SINGLE-validator points at that one confidence. Each
#                      region is normalized against its own n (declared in the
#                      subtitle/legend) since single points outnumber pairs ~4x.
#   Figure 2 (vote grid): each node (grid AND band) holds a 3x2 block of cells,
#                      rows = reference vote (both Retama / 1 vs 1 / both
#                      No_retama for pairs; for the band, a lone Retama/No_retama
#                      label reuses the top/bottom row and the middle "1 vs 1"
#                      row is structurally empty -- impossible with 1 observer),
#                      cols = map class (Retama / No_retama). Number = % WITHIN
#                      that block (6 cells sum to 100%), number size = absolute
#                      count, colour = blue agreement / red disagreement / grey
#                      no-consensus (band never has grey).
CONF_LEVELS <- 1:10          # conf 0 ("no determinable") is dropped upstream
Y_BAND      <- -0.45         # y-position of the single-validator marginal band

VOTE_LEVELS <- c("both_retama", "split", "both_noret")   # block rows, top->bottom
VOTE_DY     <- c(both_retama = 0.31, split = 0, both_noret = -0.31)
MAP_DX      <- c("1" = -0.155, "0" = 0.155)              # map Retama left, No_retama right
CELL_W      <- 0.30
CELL_H      <- 0.30
STATUS_COLS <- c(correct = "#2166AC", incorrect = "#B2182B", split = "grey55")
STATUS_LABS <- c(correct   = "mapa = referencia",
                 incorrect = "mapa != referencia",
                 split     = "1 vs 1 (sin consenso)")
Y_BREAKS <- c(Y_BAND, CONF_LEVELS)
Y_LABS   <- c("1 val.", as.character(CONF_LEVELS))
Y_LIMS   <- c(Y_BAND - 0.62, 10.6)

bubble_plots <- list(); votegrid_plots <- list()
bubble_tabs  <- list(); votegrid_tabs  <- list()

for (yr in YEARS) {
  d    <- pair_df[pair_df$year == yr, ]
  ntot <- nrow(d)
  # single-validator points usable this year: only 1 label, so 1 confidence.
  # reuses the vote/status machinery below by mapping their lone class onto
  # the top ("both_retama") / bottom ("both_noret") row of the same 3-row scheme.
  s <- usable[usable$year == yr & usable$status == "single", ]
  s$vote <- ifelse(s$ref == 1, "both_retama", "both_noret")
  ns <- nrow(s)

  ## --- Figure 1: bubble lattice (grid = pairs, band = singles) ------------
  bub <- expand.grid(A = CONF_LEVELS, B = CONF_LEVELS)
  bub$n   <- mapply(function(a, b) sum(d$x == a & d$y == b), bub$A, bub$B)
  bub$pct <- if (ntot) 100 * bub$n / ntot else NA_real_
  bub$region <- "pair"

  band1 <- data.frame(A = CONF_LEVELS, B = Y_BAND)
  band1$n   <- vapply(CONF_LEVELS, function(cv) sum(s$conf_single == cv, na.rm = TRUE), integer(1))
  band1$pct <- if (ns) 100 * band1$n / ns else NA_real_
  band1$region <- "single"

  bub_all <- rbind(bub, band1)
  bub_all$year <- yr
  bubble_tabs[[as.character(yr)]] <- bub_all

  bub_nz <- bub_all[bub_all$n > 0, ]
  bubble_plots[[as.character(yr)]] <- ggplot(bub_nz, aes(x = A, y = B, size = pct)) +
    geom_hline(yintercept = 0.35, color = "grey75", linewidth = 0.3) +
    geom_point(shape = 21, fill = "#4292C6", color = "grey25",
               alpha = 0.85, stroke = 0.4) +
    scale_size_area(max_size = 14, name = "% dentro de\nsu grupo\n(pares o\nsingles)") +
    scale_x_continuous(breaks = CONF_LEVELS, limits = c(0.4, 10.6)) +
    scale_y_continuous(breaks = Y_BREAKS, labels = Y_LABS, limits = Y_LIMS) +
    labs(title = paste0("Confianza conjunta de los evaluadores - ", yr),
         subtitle = paste0("Grilla (n=", ntot, " pares) vs banda '1 val.' (n=", ns,
                           " puntos de 1 validador); % dentro de cada region"),
         caption = "Cada region (grilla / banda) esta normalizada contra su propio total, no un total comun.",
         x = "Confianza evaluador A (1-10)",
         y = "Confianza evaluador B (1-10)  /  banda: confianza del unico validador") +
    coord_equal() + theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  ## --- Figure 2: 6-cell vote/map block per lattice node (grid + band) -----
  cells <- expand.grid(A = CONF_LEVELS, B = CONF_LEVELS, vote = VOTE_LEVELS,
                       map = c(1, 0), stringsAsFactors = FALSE)
  cells$n <- mapply(function(a, b, v, m)
      sum(d$x == a & d$y == b & as.character(d$vote) == v & d$map == m),
      cells$A, cells$B, cells$vote, cells$map)
  blk <- aggregate(n ~ A + B, data = cells, FUN = sum)
  names(blk)[names(blk) == "n"] <- "n_block"
  cells <- merge(cells, blk, by = c("A", "B"))
  cells$region <- "pair"

  cells_band <- expand.grid(A = CONF_LEVELS, B = Y_BAND, vote = VOTE_LEVELS,
                            map = c(1, 0), stringsAsFactors = FALSE)
  cells_band$n <- mapply(function(cv, v, m)
      sum(s$conf_single == cv & s$vote == v & s$map == m, na.rm = TRUE),
      cells_band$A, cells_band$vote, cells_band$map)
  blk_b <- aggregate(n ~ A, data = cells_band, FUN = sum)
  names(blk_b)[names(blk_b) == "n"] <- "n_block"
  cells_band <- merge(cells_band[, setdiff(names(cells_band), "n_block")], blk_b, by = "A")
  cells_band$region <- "single"

  cells_full <- rbind(cells, cells_band)
  cells_full$pct <- ifelse(cells_full$n_block > 0, 100 * cells_full$n / cells_full$n_block, NA_real_)
  cells_full$status <- ifelse(cells_full$vote == "split", "split",
                         ifelse((cells_full$vote == "both_retama") == (cells_full$map == 1),
                                "correct", "incorrect"))
  cells_full$cx <- cells_full$A + MAP_DX[as.character(cells_full$map)]
  cells_full$cy <- cells_full$B + VOTE_DY[cells_full$vote]
  cells_full$year <- yr
  votegrid_tabs[[as.character(yr)]] <- cells_full[, c("year", "region", "A", "B", "vote", "map",
                                                       "n", "n_block", "pct", "status")]

  frame_cells <- cells_full[cells_full$n_block > 0, ]  # draw all 6 cells of live blocks
  lab_cells   <- cells_full[cells_full$n > 0, ]         # label only non-empty cells
  lab_cells$lab <- ifelse(lab_cells$pct < 1, "<1", as.character(round(lab_cells$pct)))
  blocks <- unique(cells_full[cells_full$n_block > 0, c("A", "B")])

  votegrid_plots[[as.character(yr)]] <- ggplot() +
    geom_hline(yintercept = 0.35, color = "grey85", linewidth = 0.3) +
    geom_tile(data = blocks, aes(x = A, y = B),
              width = 0.64, height = 0.96, fill = "grey97", color = NA) +
    geom_tile(data = frame_cells, aes(x = cx, y = cy),
              width = CELL_W, height = CELL_H,
              fill = NA, color = "grey80", linewidth = 0.18) +
    geom_text(data = lab_cells,
              aes(x = cx, y = cy, label = lab, size = n, color = status),
              fontface = "bold") +
    scale_size(range = c(1.7, 4.6), name = "n de puntos\nen la celda") +
    scale_color_manual(values = STATUS_COLS, labels = STATUS_LABS,
                       name = "Acuerdo mapa /\nreferencia") +
    scale_x_continuous(breaks = CONF_LEVELS, limits = c(0.4, 10.6)) +
    scale_y_continuous(breaks = Y_BREAKS, labels = Y_LABS, limits = Y_LIMS) +
    labs(title = paste0("Voto de la referencia vs clase del mapa, por confianza - ", yr),
         subtitle = paste0("Grilla: n = ", ntot, " puntos par-validados. Banda '1 val.': n = ", ns,
                           " puntos con 1 solo validador (fila '1 vs 1' vacia por construccion)."),
         caption = paste("Bloque: columnas = mapa (izq. Retama / der. No_retama);",
                         "filas = referencia (arriba Retama, medio 1 vs 1 solo en la grilla,",
                         "abajo No_retama). Cada bloque (grilla o banda) suma 100% por si mismo."),
         x = "Confianza evaluador A (1-10)",
         y = "Confianza evaluador B (1-10)  /  banda: confianza del unico validador") +
    guides(size  = guide_legend(order = 1,
             override.aes = list(label = "12", color = "grey30")),
           color = guide_legend(order = 2,
             override.aes = list(label = "50", size = 4.2))) +
    coord_equal() + theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.2, color = "grey92"))

  ggsave(file.path(OUT_DIR, sprintf("09v_confLattice_bubble_%d.png", yr)),
         bubble_plots[[as.character(yr)]], width = 9, height = 7.3, dpi = 150)
  ggsave(file.path(OUT_DIR, sprintf("09v_confLattice_voteGrid_%d.png", yr)),
         votegrid_plots[[as.character(yr)]], width = 11, height = 10.5, dpi = 150)
  ggsave(file.path(OUT_DIR, sprintf("09v_confLattice_voteGrid_%d.pdf", yr)),
         votegrid_plots[[as.character(yr)]], width = 11, height = 10.5)
}
bubble_df   <- do.call(rbind, bubble_tabs)
votegrid_df <- do.call(rbind, votegrid_tabs)

# ---- 8. EXPORT TABLES -------------------------------------------------------
for (yr in YEARS) {
  cc <- conf_by_year[[as.character(yr)]]
  cm <- data.frame(
    "map\\ref" = c("map_Retama", "map_No_retama"),
    ref_Retama    = c(cc["TP"], cc["FN"]),
    ref_No_retama = c(cc["FP"], cc["TN"]), check.names = FALSE)
  write.csv(cm, file.path(OUT_DIR, sprintf("09v_confusion_%d.csv", yr)), row.names = FALSE)
  write.csv(perstratum[[as.character(yr)]],
            file.path(OUT_DIR, sprintf("09v_perStratum_%d.csv", yr)), row.names = FALSE)
}
write.csv(metrics_uw_df, file.path(OUT_DIR, "09v_metrics_unweighted.csv"), row.names = FALSE)
write.csv(metrics_ci_df, file.path(OUT_DIR, "09v_metrics_unweighted_ci.csv"), row.names = FALSE)
write.csv(design_df,     file.path(OUT_DIR, "09v_designCheck_perStratum.csv"), row.names = FALSE)
if (!is.null(design_oa_df))
  write.csv(design_oa_df, file.path(OUT_DIR, "09v_designCheck_OA.csv"), row.names = FALSE)
write.csv(nearfar_df,    file.path(OUT_DIR, "09v_nearfar_noretama.csv"),   row.names = FALSE)
write.csv(kappa_df,      file.path(OUT_DIR, "09v_kappa.csv"),              row.names = FALSE)
write.csv(qc_df,         file.path(OUT_DIR, "09v_qc.csv"),                 row.names = FALSE)
write.csv(bubble_df,     file.path(OUT_DIR, "09v_confLattice_bubble_counts.csv"),   row.names = FALSE)
write.csv(votegrid_df,   file.path(OUT_DIR, "09v_confLattice_voteGrid_counts.csv"), row.names = FALSE)
if (!is.null(olofsson_df)) {
  write.csv(olofsson_df,      file.path(OUT_DIR, "09v_metrics_olofsson.csv"),  row.names = FALSE)
  write.csv(olofsson_area_df, file.path(OUT_DIR, "09v_area_olofsson.csv"),     row.names = FALSE)
}

# ---- 9. HUMAN-READABLE DIGEST (for pass-2 narrative) ------------------------
fmt <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(as.numeric(x), format = "f", digits = d))
dg <- c()
dg <- c(dg, "# 09v validation — metrics digest", "",
        sprintf("Generated: %s", format(Sys.time())),
        sprintf("Filled rows loaded: %s",
                paste(sprintf("%s=%d", VALIDATORS,
                      sapply(VALIDATORS, function(v) sum(sheets$validator == v))), collapse = ", ")),
        sprintf("Master-key points: %d  (S1=%d, S2=%d, S3=%d)",
                nrow(mk), sum(mk$stratum==1), sum(mk$stratum==2), sum(mk$stratum==3)),
        sprintf("Paired (2-observer) pxids: %d", length(pair_pxids)),
        sprintf("Wi (area weights): %s",
                if (all(!is.na(Wi))) paste(sprintf("%s=%.4f", names(Wi), Wi), collapse=", ")
                else "NOT SET -> Olofsson block skipped"),
        "")
for (yr in YEARS) {
  cc <- conf_by_year[[as.character(yr)]]
  m  <- metrics_uw_df[metrics_uw_df$year == yr, ]
  q  <- qc_df[qc_df$year == yr, ]
  dg <- c(dg, sprintf("## Year %d", yr),
    sprintf("Confusion (map x ref): TP=%d FP=%d FN=%d TN=%d  (N usable=%d)",
            cc["TP"], cc["FP"], cc["FN"], cc["TN"], sum(cc)),
    sprintf("QC: disagreements discarded=%d | conf0 labels excluded=%d",
            q$n_disagree, q$n_conf0_excluded),
    "Crude (sample-based, Wilson 95% CI, unweighted):",
    "")
  ciy <- metrics_ci_df[metrics_ci_df$year == yr, ]
  for (i in seq_len(nrow(ciy))) dg <- c(dg, sprintf(
    "  %-14s = %s [%s, %s]  (SE=%s; k=%d/n=%d)",
    ciy$metric[i], fmt(ciy$est[i]), fmt(ciy$lo[i]), fmt(ciy$hi[i]),
    fmt(ciy$se[i], 4), ciy$k[i], ciy$n[i]))
  dg <- c(dg, "")
  if (!is.null(olofsson_df)) {
    o <- olofsson_df[olofsson_df$year == yr, ]; a <- olofsson_area_df[olofsson_area_df$year == yr, ]
    dg <- c(dg, "Olofsson (area-weighted, 95% CI, SE from the stratified variance):",
      sprintf("  OA         = %s [%s, %s]  (SE=%s)", fmt(o$OA), fmt(o$OA_lo), fmt(o$OA_hi), fmt(o$OA_se, 4)),
      sprintf("  UA_retama  = %s [%s, %s]  (SE=%s)   CE_retama = %s [%s, %s]",
              fmt(o$UA_retama), fmt(o$UA_retama_lo), fmt(o$UA_retama_hi), fmt(o$UA_retama_se, 4),
              fmt(o$CE_retama), fmt(o$CE_retama_lo), fmt(o$CE_retama_hi)),
      sprintf("  PA_retama  = %s [%s, %s]  (SE=%s)   OE_retama = %s [%s, %s]",
              fmt(o$PA_retama), fmt(o$PA_retama_lo), fmt(o$PA_retama_hi), fmt(o$PA_retama_se, 4),
              fmt(o$OE_retama), fmt(o$OE_retama_lo), fmt(o$OE_retama_hi)),
      sprintf("  UA_noret   = %s [%s, %s]  (SE=%s)   CE_noret  = %s [%s, %s]",
              fmt(o$UA_noret), fmt(o$UA_noret_lo), fmt(o$UA_noret_hi), fmt(o$UA_noret_se, 4),
              fmt(o$CE_noret), fmt(o$CE_noret_lo), fmt(o$CE_noret_hi)),
      sprintf("  PA_noret   = %s [%s, %s]  (SE=%s)   OE_noret  = %s [%s, %s]",
              fmt(o$PA_noret), fmt(o$PA_noret_lo), fmt(o$PA_noret_hi), fmt(o$PA_noret_se, 4),
              fmt(o$OE_noret), fmt(o$OE_noret_lo), fmt(o$OE_noret_hi)),
      sprintf("  Area(retama) proportion = %s [%s, %s]  (SE=%s)",
              fmt(a$area_prop_retama), fmt(a$area_lo), fmt(a$area_hi), fmt(a$area_se, 4)),
      "")
  }
}

# --- design check: achieved precision vs. what the sampling design promised ---
dg <- c(dg, "## Design check — achieved SE vs. Olofsson target SE",
  sprintf("Design (09v_valGenerator.js): EXP_UA = %s | TARGET_SE = %s | FLOOR = %s",
          paste(DESIGN$exp_UA, collapse = "/"), paste(DESIGN$target_SE, collapse = "/"),
          paste(DESIGN$floor_n, collapse = "/")),
  sprintf("n: formula = %s -> design = %s -> actually drawn = %s",
          paste(DESIGN$n_formula, collapse = "/"), paste(DESIGN$n_design, collapse = "/"),
          paste(DESIGN$n_drawn, collapse = "/")),
  "")
for (yr in YEARS) {
  dg <- c(dg, sprintf("### %d (per stratum)", yr))
  d <- design_df[design_df$year == yr, ]
  for (i in seq_len(nrow(d))) dg <- c(dg, sprintf(
    "  %-19s n %d/%d loaded (%s%%) | p_correct=%s (design assumed %s) | SE now=%s (%sx target %s) | SE at full n=%s | n needed at observed p=%s | target met now/full: %s/%s",
    d$label[i], d$n_loaded[i], d$n_drawn[i], fmt(100 * d$pct_loaded[i], 0),
    fmt(d$obs_p_correct[i]), fmt(d$exp_UA[i], 2), fmt(d$se_achieved[i], 4),
    fmt(d$se_ratio_vs_target[i], 2), fmt(d$target_SE[i], 2),
    fmt(d$se_at_n_drawn[i], 4), d$n_required_at_obs_p[i],
    ifelse(d$meets_target_now[i], "YES", "no"), ifelse(d$meets_target_full[i], "YES", "no")))
  if (!is.null(design_oa_df)) {
    o <- design_oa_df[design_oa_df$year == yr, ]
    dg <- c(dg, sprintf(
      "  OA: SE design=%s | SE achieved=%s (%sx design) | SE at full loading=%s",
      fmt(o$se_OA_design, 4), fmt(o$se_OA_achieved, 4),
      fmt(o$ratio_achieved_vs_design, 2), fmt(o$se_OA_at_full_loading, 4)))
  } else {
    dg <- c(dg, "  OA: SE(OA) comparison needs Wi (area weights) — not set, skipped.")
  }
  dg <- c(dg, "")
}
for (yr in YEARS) {
  k <- kappa_df[kappa_df$year == yr, ]
  dg <- c(dg, sprintf("## Inter-observer %d: n_pairs=%d  agreement=%s  kappa=%s  (both_retama=%d, split=%d, both_noret=%d)",
    yr, k$n_pairs, fmt(k$pct_agreement), fmt(k$kappa), k$both_retama, k$split, k$both_noret))
}
dg <- c(dg, "", "## Confidence-lattice figures (figures 1 & 2): grid (pairs) + marginal band (singles)")
for (yr in YEARS) {
  b <- bubble_tabs[[as.character(yr)]]; v <- votegrid_tabs[[as.character(yr)]]
  bp <- b[b$region == "pair", ]; bs <- b[b$region == "single", ]
  vp <- v[v$region == "pair", ]; vs <- v[v$region == "single", ]
  ntot <- sum(bp$n); nsin <- sum(bs$n)
  top <- bp[order(-bp$n), ][1:5, ]
  st  <- tapply(vp$n, vp$status, sum)
  sts <- tapply(vs$n, vs$status, sum)
  qc  <- qc_df[qc_df$year == yr, ]
  dg <- c(dg,
    sprintf("### %d", yr),
    sprintf("  GRID (paired points): n=%d, %d of the 100 (A,B) combinations occupied",
            ntot, sum(bp$n > 0)),
    sprintf("  top combinations (A,B): %s",
            paste(sprintf("(%d,%d) n=%d [%s%%]", top$A, top$B, top$n,
                          fmt(top$pct, 1)), collapse = ", ")),
    sprintf("  both-conf>=8 share = %s%% | both-conf<=5 share = %s%%",
            fmt(100 * sum(bp$n[bp$A >= 8 & bp$B >= 8]) / ntot, 1),
            fmt(100 * sum(bp$n[bp$A <= 5 & bp$B <= 5]) / ntot, 1)),
    sprintf("  grid cell status totals: correct=%d (%s%%), incorrect=%d (%s%%), split=%d (%s%%)",
            st[["correct"]],   fmt(100 * st[["correct"]]   / ntot, 1),
            st[["incorrect"]], fmt(100 * st[["incorrect"]] / ntot, 1),
            st[["split"]],     fmt(100 * st[["split"]]     / ntot, 1)),
    sprintf("  mean confidence of A/B where mapa!=referencia: %s / %s ; where coincide: %s / %s",
            fmt(mean(pair_df$x[pair_df$year == yr &
                     ((pair_df$vote == "both_retama") != (pair_df$map == 1)) &
                      pair_df$vote != "split"]), 2),
            fmt(mean(pair_df$y[pair_df$year == yr &
                     ((pair_df$vote == "both_retama") != (pair_df$map == 1)) &
                      pair_df$vote != "split"]), 2),
            fmt(mean(pair_df$x[pair_df$year == yr &
                     ((pair_df$vote == "both_retama") == (pair_df$map == 1)) &
                      pair_df$vote != "split"]), 2),
            fmt(mean(pair_df$y[pair_df$year == yr &
                     ((pair_df$vote == "both_retama") == (pair_df$map == 1)) &
                      pair_df$vote != "split"]), 2)),
    sprintf("  BAND (single-validator points): n=%d  (n_single_design=%d by-design + n_single_partial=%d pair-still-incomplete)",
            nsin, qc$n_single_design, qc$n_single_partial),
    sprintf("  band cell status totals: correct=%d (%s%%), incorrect=%d (%s%%)",
            ifelse(is.na(sts[["correct"]]), 0, sts[["correct"]]),
            fmt(100 * ifelse(is.na(sts[["correct"]]), 0, sts[["correct"]]) / nsin, 1),
            ifelse(is.na(sts[["incorrect"]]), 0, sts[["incorrect"]]),
            fmt(100 * ifelse(is.na(sts[["incorrect"]]), 0, sts[["incorrect"]]) / nsin, 1)),
    sprintf("  check: n_pairs(%d) = n_agree(%d, correct+incorrect above) + n_disagree(%d, qc, discarded from usable): %s",
            ntot, st[["correct"]] + st[["incorrect"]], qc$n_disagree,
            ifelse(st[["correct"]] + st[["incorrect"]] + qc$n_disagree == ntot, "OK", "MISMATCH")),
    sprintf("  check: n_agree(%d) + n_single_design(%d) + n_single_partial(%d) = %d, n_usable(%d): %s",
            st[["correct"]] + st[["incorrect"]], qc$n_single_design, qc$n_single_partial,
            st[["correct"]] + st[["incorrect"]] + qc$n_single_design + qc$n_single_partial, qc$n_usable,
            ifelse(st[["correct"]] + st[["incorrect"]] + qc$n_single_design + qc$n_single_partial == qc$n_usable,
                   "OK", "MISMATCH")),
    "")
}
writeLines(dg, file.path(OUT_DIR, "09v_metrics_digest.md"))

cat("compute done. Outputs in:", OUT_DIR, "\n")
cat(readLines(file.path(OUT_DIR, "09v_metrics_digest.md")), sep = "\n")
