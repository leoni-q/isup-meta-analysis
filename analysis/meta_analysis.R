###############################################################################
# Evolution of ISUP Grade Group 1 Detection Across Prostate Cancer Screening
# Strategies: meta-analysis per protocol.md (section 8, Statistical Analysis)
#
# Models per protocol:
#   - Risk ratios: inverse-variance, DerSimonian-Laird tau, Hartung-Knapp CI,
#     0.5 continuity correction for zero cells, random effects.
#   - Proportions: generalized linear mixed models (logit), random effects.
#   - Heterogeneity: Cochran's Q (p < 0.05 significant), plus tau^2 / I^2.
#   - Leave-one-out sensitivity analysis on every pooled estimate.
#   - Publication bias: NOT assessed (no analysis reaches k >= 10, the
#     protocol-specified threshold).
###############################################################################

suppressPackageStartupMessages({
  library(meta)
  library(metafor)
})

settings.meta(digits = 3, digits.pval = 4)

BASE <- "/Users/lukaszkujawski/Code/ania/meta_analysis/analysis"
OUT  <- file.path(BASE, "output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

LOG <- file.path(OUT, "results.txt")
cat("", file = LOG)
say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  cat(txt, "\n", file = LOG, sep = "", append = TRUE)
}
rule <- function(title = NULL) {
  say(strrep("=", 78))
  if (!is.null(title)) { say(title); say(strrep("=", 78)) }
}

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
arm <- read.csv(file.path(BASE, "arm_level_data.csv"), stringsAsFactors = FALSE)

STRAT_LEVELS <- c("No screening / usual care", "PSA + DRE", "PSA",
                  "PSA -> MRI", "Biomarker -> MRI", "MRI-first")
arm$strategy <- factor(arm$strategy, levels = STRAT_LEVELS)
stopifnot(!any(is.na(arm$strategy)))

# internal consistency: gg1 + cs <= pca in every arm
stopifnot(all(arm$gg1 + arm$cs <= arm$pca))

# Wide (trial-level) frame for the randomized comparisons.
# RRs always use ITT (randomized / invited) denominators - the randomized
# comparison must not be conditioned on being screened.
i <- arm[arm$arm == "I", ]
c_ <- arm[arm$arm == "C", ]
i <- i[order(i$trial), ]; c_ <- c_[order(c_$trial), ]
stopifnot(identical(i$trial, c_$trial))

tri <- data.frame(
  trial   = i$trial,
  year    = i$year,
  family  = i$family,
  label   = paste0(i$trial, " ", i$year),
  n.e = i$n_itt,  n.c = c_$n_itt,
  gg1.e = i$gg1,  gg1.c = c_$gg1,
  cs.e  = i$cs,   cs.c  = c_$cs,
  pca.e = i$pca,  pca.c = c_$pca,
  bx.e  = i$bx,   bx.c  = c_$bx,
  stringsAsFactors = FALSE
)
tri$fam_lab <- ifelse(
  tri$family == "A",
  "A. Screening invitation vs no screening / usual care",
  "B. Diagnostic / biopsy strategy (both arms screened)"
)
tri <- tri[order(tri$family, tri$year, tri$trial), ]

# ---------------------------------------------------------------------------
# Descriptives
# ---------------------------------------------------------------------------
rule("TABLE 1. Extracted data, internal consistency check")
desc <- arm[, c("trial", "arm", "strategy", "n_itt", "n_screened",
                "gg1", "cs", "pca", "bx")]
desc$gg1_pct_screened <- round(100 * desc$gg1 / desc$n_screened, 3)
desc$cs_pct_screened  <- round(100 * desc$cs  / desc$n_screened, 3)
desc$gg1_of_pca_pct   <- round(100 * desc$gg1 / desc$pca, 1)
desc$bx_pct_screened  <- round(100 * desc$bx  / desc$n_screened, 2)
print(desc, row.names = FALSE)
capture.output(print(desc, row.names = FALSE), file = LOG, append = TRUE)
write.csv(desc, file.path(OUT, "table1_descriptives.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
het <- function(m) {
  # GLMM fits return Q / df.Q / pval.Q of length 2: [1] Cochran's Q (Wald-type),
  # [2] the likelihood-ratio test of heterogeneity. Protocol specifies Cochran's Q.
  Q <- m$Q[1]; dfQ <- m$df.Q[1]; pQ <- m$pval.Q[1]
  q <- sprintf("Cochran's Q = %.2f, df = %s, p = %s", Q,
               ifelse(is.na(dfQ), "NA", as.character(round(dfQ))),
               ifelse(is.na(pQ), "NA", format.pval(pQ, digits = 4, eps = 1e-4)))
  i2 <- if (is.na(m$I2[1])) "I2 = NA" else sprintf("I2 = %.1f%%", 100 * m$I2[1])
  t2 <- if (is.na(m$tau2[1])) "tau2 = NA" else sprintf("tau2 = %.4f", m$tau2[1])
  sig <- if (!is.na(pQ) && pQ < 0.05) "  ** significant heterogeneity **" else ""
  extra <- if (length(m$pval.Q) > 1 && !is.na(m$pval.Q[2]))
    sprintf("; LRT of heterogeneity p = %s",
            format.pval(m$pval.Q[2], digits = 4, eps = 1e-4)) else ""
  paste0("  Heterogeneity: ", q, "; ", i2, "; ", t2, extra, sig)
}
first <- function(x) if (length(x) == 0) NA else x[1]

# Random-effects RR: inverse variance + DL tau + Hartung-Knapp CI, incr 0.5
fit_rr <- function(ee, ne, ec, nc, data, title, subgroup = TRUE) {
  d <- data[!is.na(data[[ee]]) & !is.na(data[[ec]]), ]
  args <- list(event.e = d[[ee]], n.e = d[[ne]], event.c = d[[ec]], n.c = d[[nc]],
               studlab = d$label, sm = "RR", method = "Inverse",
               method.tau = "DL", method.random.ci = "HK",
               incr = 0.5, method.incr = "only0",
               random = TRUE, common = FALSE, title = title)
  if (subgroup && length(unique(d$fam_lab)) > 1) {
    args$subgroup <- d$fam_lab
    args$subgroup.name <- "Contrast family"
    args$test.subgroup <- TRUE
  }
  do.call(metabin, args)
}

report_rr <- function(m, title) {
  rule(title)
  say(sprintf("k = %d trials | model: random effects, inverse variance, DL tau, Hartung-Knapp CI",
              m$k))
  say(sprintf("  Pooled RR = %.3f (95%% CI %.3f to %.3f), p = %s",
              exp(m$TE.random), exp(m$lower.random), exp(m$upper.random),
              format.pval(m$pval.random, digits = 4, eps = 1e-4)))
  say(het(m))
  if (!is.null(m$subgroup)) {
    say("  Subgroup estimates (random effects):")
    for (j in seq_along(m$subgroup.levels)) {
      say(sprintf("    %-55s k=%d  RR = %.3f (%.3f to %.3f)",
                  m$subgroup.levels[j], m$k.w[j],
                  exp(m$TE.random.w[j]), exp(m$lower.random.w[j]),
                  exp(m$upper.random.w[j])))
    }
    if (!is.na(m$pval.Q.b.random))
      say(sprintf("    Test for subgroup differences: Q = %.2f, df = %d, p = %s",
                  m$Q.b.random, m$df.Q.b,
                  format.pval(m$pval.Q.b.random, digits = 4, eps = 1e-4)))
  }
  say("  Per-trial risk ratios:")
  for (j in seq_len(m$k.all)) {
    say(sprintf("    %-22s %6d/%-7d vs %6d/%-7d   RR = %.3f (%.3f to %.3f)",
                m$studlab[j], m$event.e[j], m$n.e[j], m$event.c[j], m$n.c[j],
                exp(m$TE[j]), exp(m$lower[j]), exp(m$upper[j])))
  }
  invisible(m)
}

loo_rr <- function(m, title) {
  if (m$k < 3) { say(sprintf("  Leave-one-out skipped for %s (k = %d < 3).", title, m$k)); return(invisible(NULL)) }
  s <- metainf(m, pooled = "random")
  say(sprintf("  Leave-one-out (random effects, %s):", title))
  keep <- seq_len(m$k)
  for (j in keep) {
    say(sprintf("    omitting %-22s RR = %.3f (%.3f to %.3f)  I2 = %s",
                s$studlab[j], exp(s$TE[j]), exp(s$lower[j]), exp(s$upper[j]),
                ifelse(is.na(s$I2[j]), "NA", sprintf("%.1f%%", 100 * s$I2[j]))))
  }
  invisible(s)
}

plot_rr <- function(m, file, label.left, label.right, at = c(0.2, 0.5, 1, 2, 5, 10)) {
  nsub <- if (is.null(m$subgroup)) 0 else length(m$subgroup.levels)
  h <- 0.8 + 0.28 * (m$k.all + 3 * nsub + 5)
  for (dev in c("pdf", "png")) {
    f <- file.path(OUT, paste0(file, ".", dev))
    if (dev == "pdf") pdf(f, width = 11.5, height = h)
    else png(f, width = 11.5, height = h, units = "in", res = 200)
    forest(m, layout = "RevMan5", common = FALSE, random = TRUE,
           label.e = "Intervention", label.c = "Comparator",
           label.left = label.left, label.right = label.right,
           xlim = range(at), at = at,
           print.subgroup.name = FALSE, col.subgroup = "black",
           print.I2 = TRUE, print.tau2 = TRUE, print.Q = TRUE, print.pval.Q = TRUE,
           digits = 2, col.diamond = "#2166ac", col.diamond.lines = "black",
           col.square = "#4393c3", col.square.lines = "black")
    grid::grid.text(m$title, x = 0.5, y = 0.99, just = "top",
                    gp = grid::gpar(cex = 1.05, fontface = "bold"))
    dev.off()
  }
}

# Random-effects pooled proportion via GLMM (logit); k=1 groups reported exactly
fit_prop <- function(ev, nn, data, title, subgroup = TRUE) {
  d <- data[!is.na(data[[ev]]) & !is.na(data[[nn]]), ]
  args <- list(event = d[[ev]], n = d[[nn]],
               studlab = paste0(d$trial, " ", d$year, " (", d$arm, ")"),
               sm = "PLOGIT", method = "GLMM", method.tau = "ML",
               random = TRUE, common = FALSE, title = title)
  if (subgroup) {
    args$subgroup <- droplevels(d$strategy)
    args$subgroup.name <- "Screening strategy"
  }
  suppressWarnings(do.call(metaprop, args))
}

report_prop <- function(m, title) {
  rule(title)
  say(sprintf("k = %d arms | model: random effects GLMM (logit-normal), back-transformed", m$k))
  say(sprintf("  Overall pooled proportion = %.4f%% (95%% CI %.4f%% to %.4f%%)",
              100 * plogis(m$TE.random), 100 * plogis(m$lower.random),
              100 * plogis(m$upper.random)))
  say(het(m))
  if (!is.null(m$subgroup)) {
    say("  By screening strategy (random effects):")
    for (j in seq_along(m$subgroup.levels)) {
      if (m$k.w[j] == 1) {
        idx <- which(as.character(m$subgroup) == m$subgroup.levels[j])
        ci <- binom.test(m$event[idx], m$n[idx])$conf.int
        say(sprintf("    %-28s k=1  %.4f%% (exact 95%% CI %.4f%% to %.4f%%)  [%s]",
                    m$subgroup.levels[j], 100 * m$event[idx] / m$n[idx],
                    100 * ci[1], 100 * ci[2], m$studlab[idx]))
      } else {
        say(sprintf("    %-28s k=%d  %.4f%% (95%% CI %.4f%% to %.4f%%)  Q p = %s",
                    m$subgroup.levels[j], m$k.w[j],
                    100 * plogis(m$TE.random.w[j]),
                    100 * plogis(m$lower.random.w[j]),
                    100 * plogis(m$upper.random.w[j]),
                    ifelse(is.na(first(m$pval.Q.w[j])), "NA",
                           format.pval(first(m$pval.Q.w[j]), digits = 3, eps = 1e-4))))
      }
    }
  }
  say("  Per-arm observed proportions:")
  for (j in seq_len(m$k.all)) {
    say(sprintf("    %-26s %-28s %7d/%-7d = %.4f%%",
                m$studlab[j], as.character(m$subgroup)[j],
                m$event[j], m$n[j], 100 * m$event[j] / m$n[j]))
  }
  invisible(m)
}

loo_prop <- function(m, title) {
  if (m$k < 3) return(invisible(NULL))
  s <- suppressWarnings(metainf(m, pooled = "random"))
  say(sprintf("  Leave-one-out (random effects, %s):", title))
  for (j in seq_len(m$k)) {
    say(sprintf("    omitting %-26s %.4f%% (%.4f%% to %.4f%%)",
                s$studlab[j], 100 * plogis(s$TE[j]),
                100 * plogis(s$lower[j]), 100 * plogis(s$upper[j])))
  }
  invisible(s)
}

plot_prop <- function(m, file, xlab, smlab = "Pooled proportion (%)") {
  nsub <- if (is.null(m$subgroup)) 0 else length(m$subgroup.levels)
  h <- 0.8 + 0.24 * (m$k.all + 3 * nsub + 5)
  for (dev in c("pdf", "png")) {
    f <- file.path(OUT, paste0(file, ".", dev))
    if (dev == "pdf") pdf(f, width = 11, height = h)
    else png(f, width = 11, height = h, units = "in", res = 200)
    forest(m, layout = "meta", common = FALSE, random = TRUE,
           xlab = xlab, pscale = 100, smlab = smlab,
           leftcols = c("studlab", "event", "n", "effect", "ci"),
           leftlabs = c("Trial (arm)", "Events", "N", "%", "95% CI"),
           rightcols = FALSE,
           print.subgroup.name = FALSE, col.subgroup = "black",
           print.I2 = TRUE, print.Q = TRUE, print.pval.Q = TRUE,
           digits = 2, col.diamond = "#b2182b", col.diamond.lines = "black",
           col.square = "#d6604d", col.square.lines = "black")
    grid::grid.text(m$title, x = 0.5, y = 0.99, just = "top",
                    gp = grid::gpar(cex = 1.0, fontface = "bold"))
    dev.off()
  }
}

collect <- list()
add_row <- function(analysis, scale, group, k, est, lo, hi, q, dfq, pq, i2) {
  collect[[length(collect) + 1]] <<- data.frame(
    analysis = analysis, scale = scale, group = group, k = k,
    estimate = est, ci_lower = lo, ci_upper = hi,
    Q = q, df_Q = dfq, p_Q = pq, I2_pct = i2, stringsAsFactors = FALSE)
}
harvest_rr <- function(m, name) {
  add_row(name, "RR", "Overall (not interpreted across families)", m$k,
          exp(m$TE.random), exp(m$lower.random), exp(m$upper.random),
          m$Q, m$df.Q, m$pval.Q, 100 * m$I2)
  if (!is.null(m$subgroup))
    for (j in seq_along(m$subgroup.levels))
      add_row(name, "RR", m$subgroup.levels[j], m$k.w[j],
              exp(m$TE.random.w[j]), exp(m$lower.random.w[j]), exp(m$upper.random.w[j]),
              m$Q.w[j], m$df.Q.w[j], m$pval.Q.w[j], 100 * m$I2.w[j])
}
harvest_prop <- function(m, name) {
  add_row(name, "Proportion (%)", "All arms pooled", m$k,
          100 * plogis(m$TE.random), 100 * plogis(m$lower.random),
          100 * plogis(m$upper.random), first(m$Q), first(m$df.Q),
          first(m$pval.Q), 100 * first(m$I2))
  if (!is.null(m$subgroup))
    for (j in seq_along(m$subgroup.levels)) {
      if (m$k.w[j] == 1) {
        idx <- which(as.character(m$subgroup) == m$subgroup.levels[j])
        ci <- binom.test(m$event[idx], m$n[idx])$conf.int
        add_row(name, "Proportion (%)", m$subgroup.levels[j], 1,
                100 * m$event[idx] / m$n[idx], 100 * ci[1], 100 * ci[2],
                NA, NA, NA, NA)
      } else {
        add_row(name, "Proportion (%)", m$subgroup.levels[j], m$k.w[j],
                100 * plogis(m$TE.random.w[j]), 100 * plogis(m$lower.random.w[j]),
                100 * plogis(m$upper.random.w[j]),
                first(m$Q.w[j]), first(m$df.Q.w[j]), first(m$pval.Q.w[j]),
                100 * first(m$I2.w[j]))
      }
    }
}

###############################################################################
# PART 1. COMPARATIVE ANALYSES - RISK RATIOS (ITT denominators)
###############################################################################
rule("PART 1. COMPARATIVE META-ANALYSES (RISK RATIOS, ITT DENOMINATORS)")
say("Denominators are randomized / invited N throughout, preserving the")
say("randomized comparison. Subgrouped by contrast family; the overall pooled")
say("RR across both families is computed but NOT interpreted (protocol s.8/s.9).")

m_gg1 <- fit_rr("gg1.e", "n.e", "gg1.c", "n.c", tri,
                "Primary outcome: ISUP grade group 1 detection per man randomized/invited")
report_rr(m_gg1, "1.1 PRIMARY OUTCOME - ISUP GG1 detection per man randomized/invited (RR)")
loo_rr(m_gg1, "GG1 detection RR")
harvest_rr(m_gg1, "GG1 detection per man (RR)")
plot_rr(m_gg1, "fig1_rr_gg1_detection", "Favours intervention", "Favours comparator")

m_cs <- fit_rr("cs.e", "n.e", "cs.c", "n.c", tri,
               "Secondary outcome: ISUP grade group >=2 detection per man randomized/invited")
report_rr(m_cs, "1.2 SECONDARY - ISUP >=2 (csPCa) detection per man randomized/invited (RR)")
loo_rr(m_cs, "csPCa detection RR")
harvest_rr(m_cs, "csPCa detection per man (RR)")
plot_rr(m_cs, "fig2_rr_cspca_detection", "Fewer csPCa", "More csPCa",
        at = c(0.5, 1, 2, 5, 10, 20, 50, 100))

m_any <- fit_rr("pca.e", "n.e", "pca.c", "n.c", tri,
                "Any prostate cancer detection per man randomized/invited")
report_rr(m_any, "1.3 CONTEXT - Any prostate cancer detection per man randomized/invited (RR)")
loo_rr(m_any, "any PCa detection RR")
harvest_rr(m_any, "Any PCa detection per man (RR)")
plot_rr(m_any, "fig3_rr_any_pca", "Fewer PCa", "More PCa",
        at = c(0.5, 1, 2, 5, 10, 25, 50))

m_mix <- fit_rr("gg1.e", "pca.e", "gg1.c", "pca.c", tri,
                "Grade mix: ISUP GG1 as a proportion of all detected prostate cancers")
report_rr(m_mix, "1.4 PRIMARY (grade mix) - GG1 as a share of all detected PCa (RR)")
loo_rr(m_mix, "GG1 share of detected PCa RR")
harvest_rr(m_mix, "GG1 share of detected PCa (RR)")
plot_rr(m_mix, "fig4_rr_gg1_share_of_pca",
        "Lower GG1 share", "Higher GG1 share",
        at = c(0.2, 0.35, 0.5, 0.7, 1, 1.5, 2))

m_bx <- fit_rr("bx.e", "n.e", "bx.c", "n.c", tri,
               "Secondary outcome: prostate biopsies performed per man randomized")
report_rr(m_bx, "1.5 SECONDARY - Biopsies performed per man randomized (RR)")
say("  NOTE: restricted to trials reporting biopsies in BOTH arms.")
say("  ERSPC, CAP and ProScreen report intervention-arm biopsies only; PLCO reports neither.")
loo_rr(m_bx, "biopsy RR")
harvest_rr(m_bx, "Biopsies per man (RR)")
plot_rr(m_bx, "fig5_rr_biopsies", "Fewer biopsies", "More biopsies",
        at = c(0.25, 0.5, 1, 2, 4))

###############################################################################
# PART 2. SINGLE-ARM POOLED PROPORTIONS BY SCREENING STRATEGY
###############################################################################
rule("PART 2. POOLED PROPORTIONS BY SCREENING STRATEGY (GLMM)")
say("Arms are classified by the triage test that triggers work-up. Denominators")
say("are the number of men actually screened where reported; control arms of")
say("invitation trials were not screened and retain the invited/randomized N.")
say("")
say("IMPORTANT CAVEAT: in ERSPC, PLCO and CAP the numerators are cancers")
say("accrued across the whole randomized arm over 9-10 years of follow-up,")
say("whereas the per-screened denominator counts only men actually tested.")
say("These per-screened rates are therefore upward-biased for those three")
say("trials. Part 3 repeats the analysis on ITT denominators.")
say("")
say("Arms per strategy group:")
for (s in levels(droplevels(arm$strategy)))
  say(sprintf("  %-28s %s", s,
              paste(paste0(arm$trial[arm$strategy == s], "-", arm$arm[arm$strategy == s]),
                    collapse = ", ")))

p_gg1 <- fit_prop("gg1", "n_screened", arm,
                  "ISUP GG1 detection rate per man screened, by screening strategy")
report_prop(p_gg1, "2.1 PRIMARY - GG1 detection rate per man screened, by strategy")
loo_prop(p_gg1, "GG1 rate per screened man")
harvest_prop(p_gg1, "GG1 rate per man screened (proportion)")
plot_prop(p_gg1, "fig6_prop_gg1_per_screened", "ISUP GG1 detected per 100 men screened (%)")

p_cs <- fit_prop("cs", "n_screened", arm,
                 "ISUP >=2 detection rate per man screened, by screening strategy")
report_prop(p_cs, "2.2 SECONDARY - csPCa detection rate per man screened, by strategy")
loo_prop(p_cs, "csPCa rate per screened man")
harvest_prop(p_cs, "csPCa rate per man screened (proportion)")
plot_prop(p_cs, "fig7_prop_cspca_per_screened", "ISUP >=2 detected per 100 men screened (%)")

p_mix <- fit_prop("gg1", "pca", arm,
                  "ISUP GG1 as a proportion of all detected prostate cancers, by strategy")
report_prop(p_mix, "2.3 PRIMARY (grade mix) - GG1 share of all detected PCa, by strategy")
loo_prop(p_mix, "GG1 share of detected PCa")
harvest_prop(p_mix, "GG1 share of detected PCa (proportion)")
plot_prop(p_mix, "fig8_prop_gg1_share_of_pca", "ISUP GG1 as % of all detected prostate cancers")

p_bx <- fit_prop("bx", "n_screened", arm,
                 "Prostate biopsies per man screened, by screening strategy")
report_prop(p_bx, "2.4 SECONDARY - Biopsies per man screened, by strategy")
say("  Arms with biopsy counts NR are excluded (ERSPC-C, PLCO both arms, CAP-C, ProScreen-C).")
loo_prop(p_bx, "biopsies per screened man")
harvest_prop(p_bx, "Biopsies per man screened (proportion)")
plot_prop(p_bx, "fig9_prop_biopsies_per_screened", "Biopsies per 100 men screened (%)")

###############################################################################
# PART 3. SENSITIVITY ANALYSES
###############################################################################
rule("PART 3. SENSITIVITY ANALYSES")

# 3.1 proportions on ITT denominators
say("")
say("--- 3.1 Proportions on ITT (randomized/invited) denominators ---")
s_gg1 <- fit_prop("gg1", "n_itt", arm, "GG1 per man randomized/invited, by strategy")
report_prop(s_gg1, "3.1a GG1 detection rate per man randomized/invited, by strategy")
harvest_prop(s_gg1, "SENS: GG1 rate per man randomized (proportion)")
plot_prop(s_gg1, "figS1_prop_gg1_per_randomized", "ISUP GG1 detected per 100 men randomized/invited (%)")

s_cs <- fit_prop("cs", "n_itt", arm, "csPCa per man randomized/invited, by strategy")
report_prop(s_cs, "3.1b csPCa detection rate per man randomized/invited, by strategy")
harvest_prop(s_cs, "SENS: csPCa rate per man randomized (proportion)")
plot_prop(s_cs, "figS2_prop_cspca_per_randomized", "ISUP >=2 detected per 100 men randomized/invited (%)")

# 3.2 exclude PROBASE (comparator is delayed screening; control csPCa=0 inferred)
say("")
say("--- 3.2 RRs excluding PROBASE (delayed-screening comparator; control csPCa = 0 inferred) ---")
noprob <- tri[tri$trial != "PROBASE", ]
for (nm in list(c("gg1.e", "gg1.c", "GG1 detection"), c("cs.e", "cs.c", "csPCa detection"),
                c("pca.e", "pca.c", "any PCa detection"))) {
  m <- fit_rr(nm[1], "n.e", nm[2], "n.c", noprob, paste0(nm[3], " (PROBASE excluded)"))
  report_rr(m, paste0("3.2 ", nm[3], " RR - PROBASE excluded"))
  harvest_rr(m, paste0("SENS no PROBASE: ", nm[3], " (RR)"))
}
m <- fit_rr("gg1.e", "pca.e", "gg1.c", "pca.c", noprob, "GG1 share of PCa (PROBASE excluded)")
report_rr(m, "3.2 GG1 share of detected PCa RR - PROBASE excluded")
harvest_rr(m, "SENS no PROBASE: GG1 share of detected PCa (RR)")

# 3.3 exclude ERSPC (Gleason counts back-calculated from percentages)
say("")
say("--- 3.3 RRs excluding ERSPC (Gleason counts back-calculated from reported percentages) ---")
noerspc <- tri[tri$trial != "ERSPC", ]
for (nm in list(c("gg1.e", "gg1.c", "GG1 detection"), c("cs.e", "cs.c", "csPCa detection"))) {
  m <- fit_rr(nm[1], "n.e", nm[2], "n.c", noerspc, paste0(nm[3], " (ERSPC excluded)"))
  report_rr(m, paste0("3.3 ", nm[3], " RR - ERSPC excluded"))
  harvest_rr(m, paste0("SENS no ERSPC: ", nm[3], " (RR)"))
}

# 3.4 MVP and PROSA on randomized rather than analysed denominators
say("")
say("--- 3.4 RRs with MVP and PROSA on randomized (not analysed) denominators ---")
tri2 <- tri
tri2$n.e[tri2$trial == "MVP"]   <- 259; tri2$n.c[tri2$trial == "MVP"]   <- 266
tri2$n.e[tri2$trial == "PROSA"] <- 374; tri2$n.c[tri2$trial == "PROSA"] <- 385
for (nm in list(c("gg1.e", "gg1.c", "GG1 detection"), c("cs.e", "cs.c", "csPCa detection"))) {
  m <- fit_rr(nm[1], "n.e", nm[2], "n.c", tri2, paste0(nm[3], " (randomized denominators)"))
  report_rr(m, paste0("3.4 ", nm[3], " RR - MVP/PROSA randomized denominators"))
  harvest_rr(m, paste0("SENS randomized denom: ", nm[3], " (RR)"))
}

# 3.5 continuity-correction sensitivity for the single zero cell (PROBASE csPCa)
say("")
say("--- 3.5 Continuity-correction sensitivity, csPCa RR (single zero cell: PROBASE control) ---")
for (inc in c(0.1, 0.5, 1)) {
  m <- suppressWarnings(metabin(tri$cs.e, tri$n.e, tri$cs.c, tri$n.c, studlab = tri$label,
               sm = "RR", method = "Inverse", method.tau = "DL",
               method.random.ci = "HK", incr = inc, method.incr = "only0",
               random = TRUE, common = FALSE, subgroup = tri$fam_lab))
  say(sprintf("  incr = %.1f : overall RR = %.3f (%.3f to %.3f); family A RR = %.3f (%.3f to %.3f)",
              inc, exp(m$TE.random), exp(m$lower.random), exp(m$upper.random),
              exp(m$TE.random.w[1]), exp(m$lower.random.w[1]), exp(m$upper.random.w[1])))
}

###############################################################################
# PART 4. PUBLICATION BIAS
###############################################################################
rule("PART 4. PUBLICATION BIAS")
say("The protocol restricts publication-bias assessment (Peters' regression and")
say("funnel plots) to analyses containing 10 or more studies. The largest")
say(sprintf("comparative analysis here contains k = %d trials, and every strategy", nrow(tri)))
say("subgroup contains fewer still. Publication bias was therefore NOT assessed,")
say("and small-study effects cannot be excluded.")

###############################################################################
# Export summary table
###############################################################################
res <- do.call(rbind, collect)
res$estimate <- round(res$estimate, 4)
res$ci_lower <- round(res$ci_lower, 4)
res$ci_upper <- round(res$ci_upper, 4)
res$Q <- round(res$Q, 3); res$p_Q <- signif(res$p_Q, 4); res$I2_pct <- round(res$I2_pct, 1)
write.csv(res, file.path(OUT, "table2_pooled_estimates.csv"), row.names = FALSE)

rule("SESSION INFO")
say(sprintf("R %s.%s | meta %s | metafor %s",
            R.version$major, R.version$minor,
            packageVersion("meta"), packageVersion("metafor")))
say(sprintf("Outputs written to %s", OUT))
