# Meta-analysis technical report

**Evolution of ISUP Grade Group 1 Detection Across Prostate Cancer Screening Strategies**

Analysis: [meta_analysis.R](meta_analysis.R) · Numeric log:
[output/results.txt](output/results.txt) · Source:
[primary-outcome table](../Prostate%20Cancer%20Screening%20-%20Table%203.%20Primary%20outcome%20data%20used%20in%20the%20main%20meta-analysis.csv)

R 4.6.0; `meta` 8.5.0; `metafor` 5.0.1.

## 1. Data and analysis units

The extraction contains nine randomized trials and 18 arms. Comparative models
use paired randomized arms; proportion models use 16 study-strategy units.

| Cohort | Included trials | Use |
| --- | --- | --- |
| Screening invitation vs no organised screening | ERSPC, PLCO, CAP, ProScreen | Pooled comparative RRs |
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
| ProScreen | 39 / 262 (14.9%) | NR | 172 / 262 (65.6%) | NR |
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
| ProScreen | 172 / 219 (78.5%) | 282 / 355 (79.4%) |
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

- Detection-rate RRs and proportions use randomized/invited denominators.
- Cancer-grade composition uses all detected prostate cancers.
- Biopsy yield uses outcome-specific biopsy denominators.
- Biopsy use is expressed per screened/assessed participant.
- Invitation and diagnostic-strategy cohorts are fitted separately.
- Strategy proportions are fitted separately and are not head-to-head comparisons.
- A comparative biopsy-yield RR requires usable denominators in both arms.

Models:

- Risk ratios: inverse-variance random effects, DerSimonian–Laird variance,
  Hartung–Knapp confidence intervals, and a 0.5 all-cell correction for zero cells.
- Proportions: random-effects logit-normal GLMMs; exact binomial intervals for
  single-study groups.
- Prediction intervals are reported for models with at least three studies.
- Leave-one-out analysis omits each eligible contribution once within its cohort
  or strategy.

## 3. Comparative risk ratios

### 3.1 GG1 detection per randomized/invited participant

| Cohort | k | RR (95% CI) | I² | 95% prediction interval |
| --- | ---: | ---: | ---: | ---: |
| Invitation vs usual care | 4 | 1.614 (1.049–2.483) | 99.0% | 0.533–4.887 |
| Diagnostic/biopsy strategy | 4 | 0.563 (0.152–2.081) | 74.6% | 0.093–3.397 |
| PROBASE timing contrast | 1 | 7.487 (1.712–32.736) | — | — |

<!-- figures: fig1_rr_gg1_detection -->

### 3.2 GG1 among all detected prostate cancers

| Cohort | k | RR (95% CI) | I² | 95% prediction interval |
| --- | ---: | ---: | ---: | ---: |
| Invitation vs usual care | 4 | 1.204 (1.003–1.444) | 95.2% | 0.840–1.725 |
| Diagnostic/biopsy strategy | 4 | 0.567 (0.302–1.064) | 68.9% | 0.152–2.113 |
| PROBASE timing contrast | 1 | 0.380 (0.198–0.729) | — | — |

<!-- figures: fig2_rr_gg1_share_of_pca -->

### 3.3 GG≥2 detection per randomized/invited participant

| Cohort | k | RR (95% CI) | I² | 95% prediction interval |
| --- | ---: | ---: | ---: | ---: |
| Invitation vs usual care | 4 | 1.117 (0.721–1.732) | 93.6% | 0.692–1.802 |
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
| Invitation vs usual care | 4 | 0.823 (0.594–1.141) | 98.3% | 0.409–1.659 |
| Diagnostic/biopsy strategy | 4 | 1.388 (1.233–1.563) | 0.0% | 1.179–1.634 |
| PROBASE timing contrast | 1 | 4.102 (0.324–51.881) | — | — |

<!-- figures: fig4_rr_gg2plus_share_of_pca -->

### 3.5 GG1 among participants biopsied

| Cohort | k | RR (95% CI) | I² | 95% prediction interval |
| --- | ---: | ---: | ---: | ---: |
| Diagnostic/biopsy strategy | 4 | 0.878 (0.601–1.281) | 16.2% | 0.512–1.504 |
| PROBASE timing contrast | 1 | 2.313 (0.554–9.650) | — | — |

<!-- figures: fig5_rr_gg1_per_biopsy -->

### 3.6 GG≥2 among participants biopsied

| Cohort | k | RR (95% CI) | I² | 95% prediction interval |
| --- | ---: | ---: | ---: | ---: |
| Diagnostic/biopsy strategy | 4 | 2.094 (1.421–3.087) | 44.5% | 1.175–3.732 |
| PROBASE timing contrast | 1 | 21.041 (1.320–335.297) | — | — |

No invitation-family biopsy-yield RR is estimable.

<!-- figures: fig6_rr_gg2plus_per_biopsy -->

## 4. Proportions by screening strategy

These are separate descriptive models, not randomized head-to-head comparisons.
`k` is the number of study-strategy units.

### 4.1 GG1 per randomized/invited participant

| Strategy | k | Proportion (95% CI) | I² |
| --- | ---: | ---: | ---: |
| No screening / usual care | 4 | 1.309% (0.355–4.706) | 99.8% |
| PSA + DRE | 1 | 5.918% (5.684–6.159) | — |
| PSA | 5 | 1.112% (0.287–4.217) | 99.8% |
| PSA → MRI | 2 | 1.706% (1.501–1.939) | 2.9% |
| Biomarker → MRI | 2 | 0.332% (0.226–0.487) | 82.2% |
| MRI-first | 2 | 2.054% (1.196–3.504) | 0.0% |

<!-- figures: fig7_prop_gg1_per_randomized -->

### 4.2 GG1 among all detected prostate cancers

| Strategy | k | Proportion (95% CI) | I² |
| --- | ---: | ---: | ---: |
| No screening / usual care | 4 | 38.549% (22.295–57.834) | 99.8% |
| PSA + DRE | 1 | 65.730% (64.120–67.314) | — |
| PSA | 5 | 45.788% (31.125–61.218) | 99.8% |
| PSA → MRI | 2 | 46.761% (42.398–51.175) | 0.0% |
| Biomarker → MRI | 2 | 16.405% (13.247–20.141) | 0.0% |
| MRI-first | 2 | 31.707% (19.387–47.266) | 0.0% |

<!-- figures: fig8_prop_gg1_share_of_pca -->

### 4.3 GG≥2 per randomized/invited participant

| Strategy | k | Proportion (95% CI) | I² |
| --- | ---: | ---: | ---: |
| No screening / usual care | 4 | 1.958% (0.999–3.803) | 99.8% |
| PSA + DRE | 1 | 2.879% (2.714–3.052) | — |
| PSA | 5 | 1.271% (0.456–3.491) | 99.3% |
| PSA → MRI | 2 | 1.943% (1.723–2.189) | 0.0% |
| Biomarker → MRI | 2 | 1.689% (0.966–2.936) | 98.3% |
| MRI-first | 2 | 4.423% (3.071–6.332) | 0.0% |

<!-- figures: fig9_prop_gg2plus_per_randomized -->

### 4.4 GG≥2 among all detected prostate cancers

| Strategy | k | Proportion (95% CI) | I² |
| --- | ---: | ---: | ---: |
| No screening / usual care | 4 | 56.901% (39.783–72.515) | 99.6% |
| PSA + DRE | 1 | 31.982% (30.427–33.567) | — |
| PSA | 5 | 51.494% (37.320–65.432) | 99.7% |
| PSA → MRI | 2 | 53.239% (48.825–57.602) | 0.0% |
| Biomarker → MRI | 2 | 81.858% (76.967–85.900) | 67.2% |
| MRI-first | 2 | 68.293% (52.734–80.613) | 0.0% |

<!-- figures: fig10_prop_gg2plus_share_of_pca -->

### 4.5 GG1 among participants biopsied

| Strategy | k | Proportion (95% CI) | I² |
| --- | ---: | ---: | ---: |
| PSA | 4 | 18.003% (13.108–24.217) | 87.9% |
| PSA → MRI | 2 | 26.798% (23.947–29.856) | 0.0% |
| Biomarker → MRI | 2 | 12.740% (10.250–15.728) | 49.6% |
| MRI-first | 2 | 20.313% (12.177–31.910) | 0.0% |

<!-- figures: fig11_prop_gg1_per_biopsy -->

### 4.6 GG≥2 among participants biopsied

| Strategy | k | Proportion (95% CI) | I² |
| --- | ---: | ---: | ---: |
| PSA | 4 | 17.665% (10.883–27.375) | 97.7% |
| PSA → MRI | 2 | 30.510% (27.526–33.668) | 0.0% |
| Biomarker → MRI | 2 | 61.613% (56.173–66.777) | 72.2% |
| MRI-first | 2 | 43.750% (32.187–56.034) | 0.0% |

<!-- figures: fig12_prop_gg2plus_per_biopsy -->

### 4.7 Biopsies per screened/assessed participant

| Strategy | k | Proportion (95% CI) | I² |
| --- | ---: | ---: | ---: |
| PSA | 5 | 6.811% (2.018–20.600) | 100.0% |
| PSA → MRI | 2 | 6.368% (5.969–6.792) | 0.0% |
| Biomarker → MRI | 2 | 3.724% (3.266–4.244) | 81.1% |
| MRI-first | 2 | 10.407% (8.229–13.079) | 0.0% |

This outcome uses screened/assessed rather than randomized/invited denominators.
For ProScreen, the calculation is 262/7,744 (3.38%), not 262/15,201 (1.7%).

<!-- figures: fig13_prop_biopsies_per_screened -->

## 5. Sensitivity analyses

| Analysis | Outcome | k | Estimate (95% CI) |
| --- | --- | ---: | ---: |
| Invitation cohort excluding CAP | GG1 / randomized | 3 | RR 1.718 (0.817–3.611) |
| Invitation cohort excluding CAP | GG1 / all cancers | 3 | RR 1.162 (0.820–1.648) |
| Invitation cohort excluding CAP | GG≥2 / randomized | 3 | RR 1.209 (0.534–2.735) |
| Invitation cohort excluding CAP | GG≥2 / all cancers | 3 | RR 0.801 (0.438–1.464) |
| Invitation cohort excluding ERSPC | GG1 / randomized | 3 | RR 1.336 (1.071–1.666) |
| Invitation cohort excluding ERSPC | GG1 / all cancers | 3 | RR 1.152 (0.862–1.541) |
| Invitation cohort excluding ERSPC | GG≥2 / randomized | 3 | RR 1.168 (0.489–2.786) |
| Invitation cohort excluding ERSPC | GG≥2 / all cancers | 3 | RR 0.902 (0.759–1.073) |
| MVP/PROSA assessed denominators | GG1 / assessed | 4 | RR 0.562 (0.152–2.083) |
| MVP/PROSA assessed denominators | GG≥2 / assessed | 4 | RR 1.291 (0.572–2.913) |

PROBASE standalone GG≥2 continuity-correction sensitivity:

| Added value | RR (95% CI) |
| ---: | ---: |
| 0.1 | 330.433 (0.666–164009.952) |
| 0.5 | 66.885 (4.098–1091.529) |
| 1.0 | 33.942 (4.647–247.929) |

## 6. Additional checks and limitations

The analysis generated 131 leave-one-out refits. Complete results are in
[output/table3_leave_one_out.csv](output/table3_leave_one_out.csv).

Publication-bias analysis was not performed because the largest model contains
five studies, below the protocol threshold of ten. Curative treatment,
procedural harms, adverse events, RoB 2, and PRISMA flow were not analyzed.

## 7. Output files

- [output/table1_descriptives.csv](output/table1_descriptives.csv): arm-level data.
- [output/table1b_study_strategy_units.csv](output/table1b_study_strategy_units.csv):
  strategy-proportion inputs.
- [output/table2_pooled_estimates.csv](output/table2_pooled_estimates.csv): pooled
  estimates, heterogeneity, and prediction intervals.
- [output/table3_leave_one_out.csv](output/table3_leave_one_out.csv): leave-one-out
  results.
- `output/fig1_*`–`output/fig6_*`: comparative and timing RR forest plots.
- `output/fig7_*`–`output/fig13_*`: strategy-proportion forest plots.
