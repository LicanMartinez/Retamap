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
# Map label per (pxid, year) = m2017 / m2025 (blank = 0 no-retama, 1 = retama).
# Strata: 1 = S1 mapped-retama, 2 = S2 near background, 3 = S3 far background.
# near/far breakdown of the no-retama class: near = S2 (+ S1 in its off-year,
# i.e. stratum 1 pixels that are map-no-retama that year); far = S3.
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

# ---- 0. CONFIG --------------------------------------------------------------
PROJECT_DIR <- "D:/Lican/uni/Investigacion/Colaboraciones_y_ayudas/Retamap"
WORK_DIR    <- file.path(PROJECT_DIR, "gee_scripts", "validation")
SHEETS_DIR  <- file.path(WORK_DIR, "sheets_filled")
OUT_DIR     <- file.path(WORK_DIR, "analysis")
MASTER_KEY  <- file.path(WORK_DIR, "09v_master_key.csv")
YEARS       <- c(2017, 2025)
VALIDATORS  <- c("lican", "sofi", "jaime")
SEED        <- 42

# Stratum AREA WEIGHTS Wi (proportion of ROI area in each stratum), printed by
# 09v_valGenerator.js section 5 ("Wi ... suma~1") in the GEE console.
# Fill all three to enable the Olofsson area-weighted estimator.
# Leave any as NA -> only the crude (sample-based) metrics are produced.
Wi <- c(S1 = NA_real_, S2 = NA_real_, S3 = NA_real_)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. LOAD & VALIDATE -----------------------------------------------------
mk <- read.csv(MASTER_KEY, stringsAsFactors = FALSE, colClasses = "character")
stopifnot(all(c("pxid", "stratum", "m2017", "m2025", "validators") %in% names(mk)))
mk$stratum <- as.integer(mk$stratum)
zero_blank <- function(x) ifelse(is.na(x) | x == "", 0L, as.integer(x))
mk$m2017 <- zero_blank(mk$m2017)
mk$m2025 <- zero_blank(mk$m2025)
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

# zone for the no-retama near/far split
usable$zone <- ifelse(usable$stratum == 3, "far",
               ifelse(usable$stratum == 2, "near",
               ifelse(usable$stratum == 1 & usable$map == 0, "near", NA)))  # S1 on-year (map=1) => NA

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

conf_by_year   <- list()
metrics_uw     <- list()
perstratum     <- list()
nearfar        <- list()
qc_by_year     <- list()
for (yr in YEARS) {
  uy <- usable[usable$year == yr, ]
  cc <- conf2(uy$map, uy$ref)
  conf_by_year[[as.character(yr)]] <- cc
  m  <- metrics_from_conf(cc)
  # near/far commission of no-retama (among map-no-retama points, by zone)
  nz <- uy[uy$map == 0, ]
  cef <- function(z) {
    s <- nz[nz$zone == z & !is.na(nz$zone), ]
    FN <- sum(s$ref == 1); TN <- sum(s$ref == 0)
    c(FN = FN, TN = TN, n = FN + TN, CE_noret = safe(FN, FN + TN))
  }
  near <- cef("near"); far <- cef("far")
  # omission of no-retama by zone: FP live in S1 on-year (near); far has no FP
  fp_all <- sum(uy$map == 1 & uy$ref == 0)
  oe_near <- safe(fp_all, fp_all + near["TN"])   # all FP are near (S1)
  oe_far  <- safe(0,        0 + far["TN"])       # = 0 by stratum construction
  metrics_uw[[as.character(yr)]] <- data.frame(year = yr, t(m),
    CE_noret_near = unname(near["CE_noret"]), CE_noret_far = unname(far["CE_noret"]),
    OE_noret_near = unname(oe_near),          OE_noret_far = unname(oe_far),
    n_near = unname(near["n"]), n_far = unname(far["n"]), row.names = NULL)
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
      if (Y == 0) return(c(R = NA, V = NA))
      R <- X / Y
      V <- 0
      for (st in strata) {
        if (nh[st] < 2) next
        xh <- xf[[st]]; yh <- yf[[st]]
        vx <- xh * (1 - xh) / (nh[st] - 1)
        vy <- yh * (1 - yh) / (nh[st] - 1)
        cxy <- (xh - xh * yh) / (nh[st] - 1)   # x subset of y => E[xy]=E[x]
        V <- V + Wv[[st]]^2 * (vx + R^2 * vy - 2 * R * cxy)
      }
      c(R = R, V = V / Y^2)
    }
    getf <- function(st, cell) if (is.null(fr[[st]])) 0 else unname(fr[[st]][cell])
    xf_UAret <- setNames(lapply(strata, getf, "TP"), strata)
    yf_UAret <- setNames(lapply(strata, function(st) getf(st,"TP") + getf(st,"FP")), strata)
    xf_PAret <- setNames(lapply(strata, getf, "TP"), strata)
    yf_PAret <- setNames(lapply(strata, function(st) getf(st,"TP") + getf(st,"FN")), strata)
    v_UAret <- ratio_var(xf_UAret, yf_UAret)
    v_PAret <- ratio_var(xf_PAret, yf_PAret)
    ci <- function(est, v) if (is.na(v)) c(NA, NA) else est + c(-1, 1) * z * sqrt(max(v, 0))
    ol[[as.character(yr)]] <- data.frame(
      year = yr,
      OA = OA, OA_lo = OA - z*sqrt(vOA), OA_hi = OA + z*sqrt(vOA),
      UA_retama = UA_ret, UA_retama_lo = ci(UA_ret, v_UAret["V"])[1], UA_retama_hi = ci(UA_ret, v_UAret["V"])[2],
      PA_retama = PA_ret, PA_retama_lo = ci(PA_ret, v_PAret["V"])[1], PA_retama_hi = ci(PA_ret, v_PAret["V"])[2],
      CE_retama = 1 - UA_ret, OE_retama = 1 - PA_ret,
      UA_noret = UA_nor, PA_noret = PA_nor,
      CE_noret = 1 - UA_nor, OE_noret = 1 - PA_nor, row.names = NULL)
    # area of retama (proportion) + CI (stratified proportion estimator)
    Aret <- sum(sapply(strata, function(st) Wv[[st]] * ph_refret[st]))
    vA <- sum(sapply(strata, function(st) if (nh[st] > 1)
      Wv[[st]]^2 * ph_refret[st] * (1 - ph_refret[st]) / (nh[st] - 1) else 0))
    ar[[as.character(yr)]] <- data.frame(year = yr,
      area_prop_retama = Aret, area_lo = Aret - z*sqrt(vA), area_hi = Aret + z*sqrt(vA),
      row.names = NULL)
  }
  olofsson_df      <- do.call(rbind, ol)
  olofsson_area_df <- do.call(rbind, ar)
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
    "Crude (sample-based):",
    sprintf("  OA = %s", fmt(m$OA)),
    sprintf("  Retama:    commission(CE)=%s  omission(OE)=%s  (UA=%s PA=%s)",
            fmt(m$CE_retama), fmt(m$OE_retama), fmt(m$UA_retama), fmt(m$PA_retama)),
    sprintf("  No_retama: commission(CE)=%s  omission(OE)=%s  (UA=%s PA=%s)",
            fmt(m$CE_noret), fmt(m$OE_noret), fmt(m$UA_noret), fmt(m$PA_noret)),
    sprintf("  No_retama near/far: CE_near=%s (n=%d)  CE_far=%s (n=%d) | OE_near=%s  OE_far=%s",
            fmt(m$CE_noret_near), m$n_near, fmt(m$CE_noret_far), m$n_far,
            fmt(m$OE_noret_near), fmt(m$OE_noret_far)),
    "")
  if (!is.null(olofsson_df)) {
    o <- olofsson_df[olofsson_df$year == yr, ]; a <- olofsson_area_df[olofsson_area_df$year == yr, ]
    dg <- c(dg, "Olofsson (area-weighted):",
      sprintf("  OA = %s [%s, %s]", fmt(o$OA), fmt(o$OA_lo), fmt(o$OA_hi)),
      sprintf("  Retama:    CE=%s OE=%s  UA=%s [%s,%s]  PA=%s [%s,%s]",
              fmt(o$CE_retama), fmt(o$OE_retama), fmt(o$UA_retama), fmt(o$UA_retama_lo),
              fmt(o$UA_retama_hi), fmt(o$PA_retama), fmt(o$PA_retama_lo), fmt(o$PA_retama_hi)),
      sprintf("  Area(retama) proportion = %s [%s, %s]", fmt(a$area_prop_retama), fmt(a$area_lo), fmt(a$area_hi)),
      "")
  }
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
