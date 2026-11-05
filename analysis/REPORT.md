# Meta-analysis technical report

**Evolution of ISUP Grade Group 1 Detection Across Prostate Cancer Screening Strategies**

Analysis: [meta_analysis.R](meta_analysis.R) · Numeric log:
[output/results.txt](output/results.txt) · Source:
[primary-outcome table](../Prostate%20Cancer%20Screening%20-%20Table%203.%20Primary%20outcome%20data%20used%20in%20the%20main%20meta-analysis.csv)

R 4.6.0; `meta` 8.5.0; `metafor` 5.0.1.

## 1. Data and analysis units

The extraction contains nine randomized trials and 18 arms. ProScreen (two arms)
is excluded from every analysis, leaving eight trials and 16 arms. Comparative
models use paired randomized arms; proportion models use 14 study-strategy units.

| Cohort | Included trials | Use |
| --- | --- | --- |
| Screening invitation vs no organised screening | ERSPC, PLCO, CAP | Pooled comparative RRs |
| Diagnostic/biopsy strategy; both arms screened | STHLM3-MRI, MVP, GÖTEBORG-2, PROSA | Pooled comparative RRs |
| Timing of screening initiation | PROBASE | Standalone estimate |

GÖTEBORG-2 arms assigned to the same broad `PSA → MRI` category are combined.
Different strategies within a trial remain separate. PROBASE contributes only
its immediate PSA arm to strategy proportions. Its delayed comparator included
6,537 DRE assessments and two cancers, but no coherent PSA-screened denominator
at this endpoint.

CAP randomized 573 practices. No ICC or cluster-adjusted grade-specific variance
was available, so its confidence intervals and weights are unadjusted; all
invitation-family outcomes are also analyzed without CAP.

### Grade-yield inputs

`NR` denotes no coherent reported denominator.

| Trial | GG1 / biopsied, I | GG1 / biopsied, C | GG≥2 / biopsied, I | GG≥2 / biopsied, C |
| --- | ---: | ---: | ---: | ---: |
| ERSPC | 4,325 / 17,543 (24.7%) | NR | 1,665 / 17,543 (9.5%) | NR |
| PLCO | NR | NR | NR | NR |
| CAP | NR | NR | NR | NR |
| PROBASE | 15 / 120 (12.5%) | 2 / 37 (5.4%) | 33 / 120 (27.5%) | 0 / 37 (0.0%) |
| STHLM3-MRI | 34 / 311 (10.9%) | 73 / 438 (16.7%) | 192 / 331 (58.0%) | 106 / 438 (24.2%) |
| MVP | 4 / 24 (16.7%) | 4 / 28 (14.3%) | 11 / 24 (45.8%) | 4 / 28 (14.3%) |
| GÖTEBORG-2 | 68 / 254 (26.8%) | 159 / 588 (27.0%) | 117 / 254 (46.1%) | 139 / 588 (23.6%) |
| PROSA | 9 / 40 (22.5%) | 4 / 20 (20.0%) | 17 / 40 (42.5%) | 7 / 20 (35.0%) |

| Trial | GG≥2 / all detected cancers, I | GG≥2 / all detected cancers, C |
| --- | ---: | ---: |
| ERSPC | 1,665 / 5,990 (27.8%) | 1,947 / 4,307 (45.2%) |
| PLCO | 1,104 / 3,452 (32.0%) | 1,120 / 2,974 (37.7%) |
| CAP | 6,649 / 12,013 (55.3%) | 8,019 / 12,958 (61.9%) |
| PROBASE | 33 / 48 (68.8%) | 0 / 2 (0.0%) |
| STHLM3-MRI | 192 / 226 (85.0%) | 106 / 179 (59.2%) |
| MVP | 11 / 15 (73.3%) | 4 / 8 (50.0%) |
| GÖTEBORG-2 | 117 / 185 (63.2%) | 139 / 298 (46.6%) |
| PROSA | 17 / 26 (65.4%) | 7 / 11 (63.6%) |

CAP is excluded from grade-per-biopsy analyses because its 5,848 biopsies concern
the first screening round, while grade counts are cumulative over 15 years.
STHLM3-MRI reports intervention GG≥2 as `192/331`, although its general biopsy
count is 311; the outcome-specific denominator is retained.

## 2. Analysis specification

- Detection-rate RRs use randomized/invited denominators, preserving the
  randomized comparison.
- Detection-rate proportions are reported on both denominators: per man screened
  (protocol s.9) and per man randomized/invited.
- Cancer-grade composition uses all detected prostate cancers.
- Biopsy yield uses outcome-specific biopsy denominators.
- Biopsy burden is reported both comparatively (RR of biopsies per
  randomized/invited participant, trials reporting both arms) and per
  screened/assessed participant.
- Invitation and diagnostic-strategy cohorts are fitted separately.
- Strategy proportions are fitted separately and are not head-to-head comparisons.

Per-screened detection rates carry a denominator caveat: in ERSPC, PLCO and CAP
the numerators accrue across the whole randomized arm over long follow-up, while
the per-screened denominator counts only men actually tested, so those rates are
upward-biased. The per-randomized/invited rates are the internally consistent
ones. Both are reported so the gap is visible rather than assumed.

Models:

- Risk ratios: inverse-variance random effects, DerSimonian–Laird variance,
  Hartung–Knapp confidence intervals, and a 0.5 all-cell correction for zero cells.
- Proportions: random-effects logit-normal GLMMs; exact binomial intervals for
  single-study groups.
- Single-trial contrasts are not pooled: they carry the trial's own estimate with
  a normal-approximation interval, and no heterogeneity statistics.
- Prediction intervals are reported for models with at least three studies.
- Leave-one-out analysis omits each contribution once, and runs only where at
  least three studies contribute; omitting one of two would leave a single trial
  rather than a pooled estimate.

## 3. Comparative risk ratios

### 3.1 GG1 detection per randomized/invited participant

| Cohort | k | RR (95% CI) | I² | 95% prediction interval |
| --- | ---: | ---: | ---: | ---: |
| Invitation vs usual care | 3 | 1.573 (0.726–3.408) | 99.3% | 0.335–7.378 |
| Diagnostic/biopsy strategy | 4 | 0.563 (0.152–2.081) | 74.6% | 0.093–3.397 |
| PROBASE timing contrast | 1 | 7.487 (1.712–32.736) | — | — |

<!-- figures: fig1_rr_gg1_detection -->

### 3.2 GG1 among all detected prostate cancers

| Cohort | k | RR (95% CI) | I² | 95% prediction interval |
| --- | ---: | ---: | ---: | ---: |
| Invitation vs usual care | 3 | 1.223 (0.952–1.572) | 96.7% | 0.750–1.997 |
| Diagnostic/biopsy strategy | 4 | 0.567 (0.302–1.064) | 68.9% | 0.152–2.113 |
| PROBASE timing contrast | 1 | 0.380 (0.198–0.729) | — | — |

<!-- figures: fig2_rr_gg1_share_of_pca -->

### 3.3 GG≥2 detection per randomized/invited participant

| Cohort | k | RR (95% CI) | I² | 95% prediction interval |
| --- | ---: | ---: | ---: | ---: |
| Invitation vs usual care | 3 | 0.993 (0.886–1.113) | 64.6% | 0.805–1.225 |
| Diagnostic/biopsy strategy | 4 | 1.292 (0.571–2.925) | 73.8% | 0.373–4.475 |
| PROBASE timing contrast | 1 | 66.885 (4.098–1091.529) | — | — |

The PROBASE comparator has zero inferred GG≥2 events; continuity-correction
sensitivity is reported below.

<!-- figures: fig3_rr_gg2plus_detection -->

### 3.4 GG≥2 among all detected prostate cancers

| Trial | Intervention strategy: GG≥2 / all cancers | Comparator strategy: GG≥2 / all cancers | Trial RR |
| --- | ---: | ---: | ---: |
| STHLM3-MRI | Biomarker → MRI: 192/226 (85.0%) | PSA: 106/179 (59.2%) | 1.435 |
| MVP | MRI-first: 11/15 (73.3%) | PSA: 4/8 (50.0%) | 1.467 |
| GÖTEBORG-2 | MRI-targeted only: 117/185 (63.2%) | Systematic ± targeted: 139/298 (46.6%) | 1.356 |
| PROSA | MRI-first: 17/26 (65.4%) | PSA-triggered MRI: 7/11 (63.6%) | 1.027 |

For STHLM3-MRI, 85.0% / 59.2% gives RR 1.435. The pooled RR is 1.388,
equivalent to a 38.8% higher GG≥2 share among detected cancers.

| Cohort | k | RR (95% CI) | I² | 95% prediction interval |
| --- | ---: | ---: | ---: | ---: |
| Invitation vs usual care | 3 | 0.776 (0.469–1.286) | 98.8% | 0.273–2.206 |
| Diagnostic/biopsy strategy | 4 | 1.388 (1.233–1.563) | 0.0% | 1.179–1.634 |
| PROBASE timing contrast | 1 | 4.102 (0.324–51.881) | — | — |

<!-- figures: fig4_rr_gg2plus_share_of_pca -->

### 3.5 GG1 among participants biopsied

| Cohort | k | RR (95% CI) | I² | 95% prediction interval |
| --- | ---: | ---: | ---: | ---: |
| Diagnostic/biopsy strategy | 3 | 1.002 (0.888–1.130) | 0.0% | 0.601–1.669 |
| PROBASE timing contrast | 1 | 2.313 (0.554–9.650) | — | — |

STHLM3-MRI was excluded from the diagnostic/biopsy-strategy analysis.

<!-- figures: fig5_rr_gg1_per_biopsy -->

### 3.6 GG≥2 among participants biopsied

| Cohort | k | RR (95% CI) | I² | 95% prediction interval |
| --- | ---: | ---: | ---: | ---: |
| Diagnostic/biopsy strategy | 3 | 1.875 (0.861–4.085) | 25.3% | 0.646–5.442 |
| PROBASE timing contrast | 1 | 21.041 (1.320–335.297) | — | — |

STHLM3-MRI was excluded from the diagnostic/biopsy-strategy analysis.
No invitation-family biopsy-yield RR is estimable.

<!-- figures: fig6_rr_gg2plus_per_biopsy -->

### 3.7 Biopsies performed per randomized/invited participant

Protocol secondary endpoint (s.6): the number of biopsies performed. This is the
randomized comparison of biopsy burden, distinct from the yield outcomes above.

| Cohort | k | RR (95% CI) | I² | 95% prediction interval |
| --- | ---: | ---: | ---: | ---: |
| Diagnostic/biopsy strategy | 4 | 0.713 (0.236–2.152) | 92.0% | 0.175–2.900 |
| PROBASE timing contrast | 1 | 3.238 (2.241–4.678) | — | — |

The pooled estimate is dominated by opposing trial effects: GÖTEBORG-2 (0.43) and
STHLM3-MRI (0.48) roughly halve biopsies, while PROSA increases them (2.06) and
MVP is neutral (0.88). With I² = 92.0% the pooled value should not be read as a
single biopsy-burden effect. No invitation-family estimate is possible: ERSPC
and CAP report intervention-arm biopsies only, and PLCO reports none.

<!-- figures: fig7_rr_biopsies_per_participant -->

## 4. Proportions by screening strategy

These are separate descriptive models, not randomized head-to-head comparisons.
`k` is the number of study-strategy units.

### 4.1 GG1 per screened participant

Protocol-specified denominator (s.9): detection per number of screened men.

| Strategy | k | Proportion (95% CI) | I² |
| --- | ---: | ---: | ---: |
| No screening / usual care | 3 | 2.702% (1.634–4.438) | 99.9% |
| PSA + DRE | 1 | 6.962% (6.688–7.244) | — |
| PSA | 5 | 1.422% (0.317–6.140) | 99.4% |
| PSA → MRI | 2 | 1.707% (1.502–1.939) | 0.8% |
| Biomarker → MRI | 1 | 0.447% (0.310–0.624) | — |
| MRI-first | 2 | 2.114% (1.231–3.606) | 0.0% |

<!-- figures: fig8_prop_gg1_per_screened -->

### 4.2 GG1 per randomized/invited participant

| Strategy | k | Proportion (95% CI) | I² |
| --- | ---: | ---: | ---: |
| No screening / usual care | 3 | 2.702% (1.634–4.438) | 99.9% |
| PSA + DRE | 1 | 5.918% (5.684–6.159) | — |
| PSA | 5 | 1.112% (0.287–4.216) | 99.8% |
| PSA → MRI | 2 | 1.706% (1.501–1.939) | 2.9% |
| Biomarker → MRI | 1 | 0.447% (0.310–0.624) | — |
| MRI-first | 2 | 2.054% (1.196–3.504) | 0.0% |

The two denominators separate where screening uptake was incomplete. For PSA +
DRE (PLCO) the rate is 6.962% per screened man against 5.918% per randomized
man. Control arms of invitation trials were never screened, so their two
denominators coincide, and the No screening / usual care rows in 4.1 and 4.2 are
identical. The Biomarker → MRI rows also match, because STHLM3-MRI, its only
contributor, reports the same screened and randomized denominator.

<!-- figures: fig9_prop_gg1_per_randomized -->

### 4.3 GG1 among all detected prostate cancers

| Strategy | k | Proportion (95% CI) | I² |
| --- | ---: | ---: | ---: |
| No screening / usual care | 3 | 46.730% (30.240–63.967) | 99.9% |
| PSA + DRE | 1 | 65.730% (64.120–67.314) | — |
| PSA | 5 | 45.788% (31.125–61.218) | 99.8% |
| PSA → MRI | 2 | 46.761% (42.398–51.175) | 0.0% |
| Biomarker → MRI | 1 | 15.044% (10.649–20.383) | — |
| MRI-first | 2 | 31.707% (19.387–47.266) | 0.0% |

<!-- figures: fig10_prop_gg1_share_of_pca -->

### 4.4 GG≥2 per screened participant

| Strategy | k | Proportion (95% CI) | I² |
| --- | ---: | ---: | ---: |
| No screening / usual care | 3 | 2.858% (2.245–3.631) | 99.6% |
| PSA + DRE | 1 | 3.387% (3.194–3.590) | — |
| PSA | 5 | 1.629% (0.478–5.398) | 99.8% |
| PSA → MRI | 2 | 1.943% (1.724–2.190) | 0.0% |
| Biomarker → MRI | 1 | 2.523% (2.183–2.901) | — |
| MRI-first | 2 | 4.553% (3.162–6.515) | 0.0% |

<!-- figures: fig11_prop_gg2plus_per_screened -->

### 4.5 GG≥2 per randomized/invited participant

| Strategy | k | Proportion (95% CI) | I² |
| --- | ---: | ---: | ---: |
| No screening / usual care | 3 | 2.858% (2.245–3.631) | 99.6% |
| PSA + DRE | 1 | 2.879% (2.714–3.052) | — |
| PSA | 5 | 1.271% (0.456–3.491) | 99.3% |
| PSA → MRI | 2 | 1.943% (1.723–2.189) | 0.0% |
| Biomarker → MRI | 1 | 2.523% (2.183–2.901) | — |
| MRI-first | 2 | 4.423% (3.071–6.332) | 0.0% |

<!-- figures: fig12_prop_gg2plus_per_randomized -->

### 4.6 GG≥2 among all detected prostate cancers

| Strategy | k | Proportion (95% CI) | I² |
| --- | ---: | ---: | ---: |
| No screening / usual care | 3 | 48.262% (36.887–59.820) | 99.7% |
| PSA + DRE | 1 | 31.982% (30.427–33.567) | — |
| PSA | 5 | 51.494% (37.320–65.432) | 99.7% |
| PSA → MRI | 2 | 53.239% (48.825–57.602) | 0.0% |
| Biomarker → MRI | 1 | 84.956% (79.617–89.351) | — |
| MRI-first | 2 | 68.293% (52.734–80.613) | 0.0% |

<!-- figures: fig13_prop_gg2plus_share_of_pca -->

### 4.7 GG1 among participants biopsied

| Strategy | k | Proportion (95% CI) | I² |
| --- | ---: | ---: | ---: |
| PSA | 3 | 18.164% (11.539–27.415) | 81.0% |
| PSA → MRI | 2 | 26.798% (23.947–29.856) | 0.0% |
| MRI-first | 2 | 20.312% (12.177–31.910) | 0.0% |

STHLM3-MRI was excluded from this outcome-specific analysis, so no Biomarker →
MRI unit remains.

<!-- figures: fig14_prop_gg1_per_biopsy -->

### 4.8 GG≥2 among participants biopsied

| Strategy | k | Proportion (95% CI) | I² |
| --- | ---: | ---: | ---: |
| PSA | 3 | 15.471% (8.370–26.834) | 94.9% |
| PSA → MRI | 2 | 30.510% (27.526–33.668) | 0.0% |
| MRI-first | 2 | 43.750% (32.187–56.034) | 0.0% |

STHLM3-MRI was excluded from this outcome-specific analysis, so no Biomarker →
MRI unit remains.

<!-- figures: fig15_prop_gg2plus_per_biopsy -->

### 4.9 Biopsies per screened/assessed participant

| Strategy | k | Proportion (95% CI) | I² |
| --- | ---: | ---: | ---: |
| PSA | 5 | 6.811% (2.018–20.599) | 100.0% |
| PSA → MRI | 2 | 6.368% (5.969–6.792) | 0.0% |
| Biomarker → MRI | 1 | 4.087% (3.654–4.557) | — |
| MRI-first | 2 | 10.406% (8.229–13.079) | 0.0% |

This outcome uses screened/assessed rather than randomized/invited denominators.
For CAP, the calculation is 5,848/75,694 (7.7%), not 5,848/189,326 (3.1%).

<!-- figures: fig16_prop_biopsies_per_screened -->

## 5. Sensitivity analyses

| Analysis | Outcome | k | Estimate (95% CI) |
| --- | --- | ---: | ---: |
| Invitation cohort excluding CAP | GG1 / randomized | 2 | RR 1.687 (0.044–64.583) |
| Invitation cohort excluding CAP | GG1 / all cancers | 2 | RR 1.199 (0.360–3.996) |
| Invitation cohort excluding CAP | GG≥2 / randomized | 2 | RR 1.022 (0.696–1.501) |
| Invitation cohort excluding CAP | GG≥2 / all cancers | 2 | RR 0.722 (0.093–5.616) |
| Invitation cohort excluding ERSPC | GG1 / randomized | 2 | RR 1.320 (0.806–2.161) |
| Invitation cohort excluding ERSPC | GG1 / all cancers | 2 | RR 1.178 (0.439–3.162) |
| Invitation cohort excluding ERSPC | GG≥2 / randomized | 2 | RR 0.964 (0.863–1.077) |
| Invitation cohort excluding ERSPC | GG≥2 / all cancers | 2 | RR 0.880 (0.651–1.191) |
| MVP/PROSA assessed denominators | GG1 detection per assessed participant | 4 | RR 0.562 (0.152–2.083) |
| MVP/PROSA assessed denominators | GG≥2 detection per assessed participant | 4 | RR 1.291 (0.572–2.913) |

PROBASE standalone GG≥2 continuity-correction sensitivity:

| Added value | RR (95% CI) |
| ---: | ---: |
| 0.1 | 330.433 (0.666–164009.952) |
| 0.5 | 66.885 (4.098–1091.529) |
| 1.0 | 33.942 (4.647–247.929) |

## 6. Additional checks and limitations

The analysis generated 97 leave-one-out refits. Complete results are in
[output/table3_leave_one_out.csv](output/table3_leave_one_out.csv). Analyses with
fewer than three contributing studies have no leave-one-out rows by design.

Publication-bias analysis was not performed because the largest model contains
five studies, below the protocol threshold of ten. The threshold is evaluated at
run time, so a future search update that pushes any model to ten or more studies
triggers Peters' regression and funnel plots automatically.

### Protocol deviations

Recorded in full at the end of [output/results.txt](output/results.txt).

1. **Stratification of comparative analyses (protocol s.8).** The protocol
   specifies pooling risk ratios "stratified by the screening strategy". A risk
   ratio is a contrast *between* two strategies, so a single strategy label
   cannot describe one. Trials are grouped instead by contrast family
   (invitation, diagnostic/biopsy strategy, timing), and families are never
   pooled together. Single-arm proportions are stratified by screening strategy
   exactly as the protocol states.
2. **Continuity correction for proportions (protocol s.8).** The protocol states
   a 0.5 correction for zero cells in both analyses. It is applied to the
   risk-ratio models. The proportion models are logit-normal GLMMs fitted by
   likelihood; they use the exact binomial likelihood and need no continuity
   correction, and adding one would bias them. None is applied there.
3. **Curative-treatment outcome (protocol s.2 and s.9).** "Percentage of patients
   receiving curative treatment over the course of follow-up" is a registered
   secondary outcome. The included reports provide no extractable arm-level
   treatment counts, so it is not analyzed. This is an outcome-level deviation.
4. **Adverse events and procedural harms (protocol s.9).** Contextual data only,
   not formal endpoints; not extracted and not summarized.
5. **Exclusion of an extracted trial.** ProScreen was extracted (its rows remain
   in `arm_level_data.csv`) but is excluded from every analysis, including the
   descriptive tables.

RoB 2 and the PRISMA flow diagram are handled outside this pipeline.

## 7. Output files

- [output/table1_descriptives.csv](output/table1_descriptives.csv): arm-level data.
- [output/table1b_study_strategy_units.csv](output/table1b_study_strategy_units.csv):
  strategy-proportion inputs.
- [output/table2_pooled_estimates.csv](output/table2_pooled_estimates.csv): pooled
  estimates, heterogeneity, and prediction intervals.
- [output/table3_leave_one_out.csv](output/table3_leave_one_out.csv): leave-one-out
  results.
- `output/fig1_*`–`output/fig7_*`: comparative and timing RR forest plots.
- `output/fig8_*`–`output/fig16_*`: strategy-proportion forest plots.

Reproducibility checks: [validate_analysis.R](validate_analysis.R) verifies the
structure of the generated outputs, and [tests/test_calculations.R](tests/test_calculations.R)
unit-tests the calculation helpers against hand-computed values and independent
`metafor` refits.
