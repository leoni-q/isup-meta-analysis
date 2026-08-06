# Meta-analysis results

**Evolution of ISUP Grade Group 1 Detection Across Prostate Cancer Screening Strategies**

Analysis performed per `protocol.md` §8 (Statistical Analysis) and §9 (Data Synthesis).
R 4.6.0, `meta` 8.5.0, `metafor` 5.0.1. Script: [meta_analysis.R](meta_analysis.R).
Full numeric output: [output/results.txt](output/results.txt).

---

## 1. Included data

Nine randomized trials, 18 arms, 782,328 men randomized/invited. All nine parse
internally consistently: ISUP GG1 + ISUP ≥2 ≤ any PCa in every arm, and the
reported "GG1 / any PCa" percentages reproduce exactly from the component counts.

The trials divide into two contrast families, which are **never pooled into a single
interpreted estimate** (protocol §8/§9; also required by the PROSA extraction note):


| Family | Trials                                                        | Contrast                                          |
| ------ | ------------------------------------------------------------- | ------------------------------------------------- |
| **A**  | ERSPC 2009, PLCO 2009, CAP 2018, PROBASE 2022, ProScreen 2024 | Screening invitation vs no screening / usual care |
| **B**  | STHLM3-MRI 2021, MVP 2022, GÖTEBORG-2 2024, PROSA 2026        | Diagnostic / biopsy strategy, both arms screened  |


Arms were classified by **the triage test that triggers work-up**:


| Strategy                  | Arms (k)                                                  |
| ------------------------- | --------------------------------------------------------- |
| No screening / usual care | ERSPC-C, PLCO-C, CAP-C, ProScreen-C (4)                   |
| PSA + DRE                 | PLCO-I (1)                                                |
| PSA                       | ERSPC-I, CAP-I, PROBASE-I, PROBASE-C, STHLM3-C, MVP-C (6) |
| PSA → MRI                 | GÖTEBORG-2-I, GÖTEBORG-2-C, PROSA-C (3)                   |
| Biomarker → MRI           | ProScreen-I, STHLM3-I (2)                                 |
| MRI-first                 | MVP-I, PROSA-I (2)                                        |


---



## 2. Primary outcome — ISUP GG1 detection



### 2.1 Comparative (risk ratios, ITT denominators) — Figure 1

Random effects, inverse variance, DerSimonian–Laird τ, Hartung–Knapp CI, 0.5 correction for zero cells.


| Family                              | k   | RR (95% CI)                     |
| ----------------------------------- | --- | ------------------------------- |
| **A. Invitation vs no screening**   | 5   | **1.76 (1.09 – 2.84)**          |
| **B. Diagnostic / biopsy strategy** | 4   | **0.56 (0.15 – 2.08)**          |
| Test for subgroup differences       |     | Q = 6.51, df = 1, **p = 0.011** |


Inviting men to screening increases GG1 detection by roughly three quarters
(all five family-A trials point the same way: RR 1.27 to 7.49). Within already-screened
populations, modern pathways point toward *reduced* GG1 detection, but the estimate
is imprecise and directionally inconsistent (STHLM3-MRI 0.31 and GÖTEBORG-2 0.43 vs
PROSA 2.32 and MVP 1.03) — the two large trials drive the reduction, the two small
ones oppose it. The families differ significantly.

Heterogeneity is extreme throughout: overall Q = 388.3, df = 8, p < 0.0001, I² = 97.9%.
Leave-one-out: family-level conclusions hold; no single trial drives the family-A result
(omitting ERSPC gives the largest shift, overall RR 1.17 → 1.02).

### 2.2 Grade mix — GG1 as a share of all detected cancers — Figure 4

This is the more comparable outcome, because it is internally normalized within each
trial and so is not distorted by differences in follow-up duration, screening rounds, or age.


| Family                              | k   | RR (95% CI)            | Excluding PROBASE      |
| ----------------------------------- | --- | ---------------------- | ---------------------- |
| A. Invitation vs no screening       | 5   | 1.08 (0.64 – 1.82)     | **1.21 (1.00 – 1.47)** |
| **B. Diagnostic / biopsy strategy** | 4   | **0.57 (0.30 – 1.06)** | 0.57 (0.30 – 1.06)     |
| Subgroup difference                 |     | p = 0.019              | p = 0.00025            |


Screening invitation shifts the detected-cancer mix *toward* low grade (ERSPC 1.32,
CAP 1.30, PLCO 1.09 — all significant). Modern diagnostic pathways shift it *away*
from low grade by about 43% relative, just short of conventional significance.

### 2.3 Pooled proportions by strategy (GLMM) — Figures 6 and 8

**GG1 as % of all detected prostate cancers** — the evolution signal:


| Strategy                    | k   | GG1 share of detected PCa (95% CI) |
| --------------------------- | --- | ---------------------------------- |
| PSA + DRE                   | 1   | 65.7% (64.1 – 67.3)                |
| PSA                         | 6   | 49.4% (35.1 – 63.9)                |
| PSA → MRI                   | 3   | 44.3% (33.9 – 55.2)                |
| No screening / usual care   | 4   | 39.8% (23.8 – 58.4)                |
| MRI-first                   | 2   | 31.7% (19.4 – 47.3)                |
| Biomarker → MRI             | 2   | 16.4% (13.2 – 20.1)                |
| Between-strategy difference |     | χ² = 325.0, df = 5, p < 0.0001     |


A clear gradient: roughly two thirds of screen-detected cancers were GG1 in the
PSA + DRE era, about half under PSA alone, and one sixth to one third under
biomarker- or MRI-triaged pathways. This is the quantitative backbone of the review's
thesis. Heterogeneity within the two MRI/biomarker groups is nil (I² = 0%, Q p = 0.43
and 0.60), so those two estimates are the most trustworthy in the table.

**GG1 detection rate per man screened** (denominators = men actually screened where reported):


| Strategy                  | k   | GG1 per 100 screened (95% CI) |
| ------------------------- | --- | ----------------------------- |
| PSA + DRE                 | 1   | 6.96% (6.69 – 7.24)           |
| MRI-first                 | 2   | 2.11% (1.23 – 3.61)           |
| PSA → MRI                 | 3   | 1.46% (0.87 – 2.45)           |
| No screening / usual care | 4   | 1.20% (0.33 – 4.30)           |
| PSA                       | 6   | 0.60% (0.08 – 4.14)           |
| Biomarker → MRI           | 2   | 0.48% (0.38 – 0.60)           |


> **These per-screened rates should not be read as strategy effects.** The PSA group
> has τ² = 5.90 and spans a 120-fold range (ERSPC-I 7.22% to PROBASE-I 0.06%), so its
> pooled estimate of 0.60% lies *below almost every arm contributing to it* — a
> mathematical artefact of logit-scale pooling under extreme heterogeneity, not a finding.
> The dominant drivers are follow-up duration (9–10 years in ERSPC/PLCO/CAP vs a single
> round elsewhere), number of screening rounds, and age at screening (45 in PROBASE),
> not the screening strategy. Only the two homogeneous groups (Biomarker → MRI, I² = 0%;
> MRI-first, I² = 0%) support direct interpretation.

---



## 3. Secondary outcomes



### 3.1 ISUP ≥2 (csPCa) detection — Figures 2 and 7


| Family                          | k   | RR (95% CI)        |
| ------------------------------- | --- | ------------------ |
| A. Invitation vs no screening   | 5   | 1.16 (0.76 – 1.77) |
| B. Diagnostic / biopsy strategy | 4   | 1.29 (0.57 – 2.93) |
| Subgroup difference             |     | p = 0.71           |


No detectable difference in clinically significant cancer detection in either family, and
no difference between families. Read alongside §2.1, this is the central asymmetry:
screening invitation raised GG1 detection by ~76% while csPCa detection was essentially
unchanged (RR 1.05, 0.99, 1.04 in ERSPC, PLCO, CAP respectively) — i.e. the excess
detection was almost entirely low grade.

Pooled csPCa per 100 men screened: MRI-first 4.55% (3.16 – 6.52), PSA + DRE 3.39%,
Biomarker → MRI 2.37% (2.14 – 2.62), PSA → MRI 1.94% (1.72 – 2.19), no screening 1.69%,
PSA 0.47%. The same denominator caveat as §2.3 applies to the PSA and no-screening rows.

### 3.2 Biopsies performed — Figures 5 and 9

Comparative analysis is limited to the 5 trials reporting biopsies in **both** arms
(ERSPC, CAP and ProScreen report the intervention arm only; PLCO reports neither):


| Family                          | k   | RR (95% CI)        |
| ------------------------------- | --- | ------------------ |
| A (PROBASE only)                | 1   | 3.24 (2.24 – 4.68) |
| B. Diagnostic / biopsy strategy | 4   | 0.71 (0.24 – 2.15) |
| Subgroup difference             |     | p = 0.00013        |


Within family B the direction splits sharply by design: the two triage-tightening trials
cut biopsies by more than half (STHLM3-MRI 0.48, GÖTEBORG-2 0.43) while the two
MRI-first trials increased them (PROSA 2.06, MVP 0.88). Pooled biopsies per 100 men
screened: MRI-first 10.41% (8.23 – 13.08), PSA → MRI 5.70%, PSA 3.80%,
Biomarker → MRI 3.72% (3.27 – 4.24).

### 3.3 Curative treatment — not analysable

The protocol lists "percentage of patients receiving curative treatment over follow-up"
as a secondary outcome. These data are not present in the extraction and this outcome
could not be analysed. Adverse-event / procedural-harm data (protocol §5, §9) are
likewise absent from the extraction.

---



## 4. Heterogeneity, sensitivity, publication bias

**Heterogeneity** was significant by Cochran's Q (p < 0.05) in every pooled analysis,
with I² from 68.9% to 99.9%. Per protocol §8 this triggered exploratory analysis, which
is what the family and strategy stratifications above provide.

**Leave-one-out** was run on every pooled estimate (see `output/results.txt`). No
conclusion depends on a single trial.

**Targeted sensitivity analyses:**


| Analysis                                                                                                            | Effect                                                                                                                                     |
| ------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| Exclude **PROBASE** (delayed-screening comparator; control csPCa = 0 was *inferred by the extractor, not reported*) | Family A GG1 RR 1.76 → 1.67 (1.12 – 2.49); family A grade-mix RR 1.08 → 1.21 (1.00 – 1.47), subgroup difference strengthens to p = 0.00025 |
| Exclude **ERSPC** (Gleason counts back-calculated from published percentages)                                       | Family A GG1 RR 1.76 → 1.50 (0.96 – 2.34) — loses significance                                                                             |
| **MVP / PROSA on randomized** rather than analysed denominators                                                     | No material change (GG1 RR 1.167 both ways)                                                                                                |
| Continuity correction **0.1 / 0.5 / 1.0** on the single zero cell                                                   | Family A csPCa RR 1.14 / 1.16 / 1.17 — robust                                                                                              |
| Proportions on **ITT rather than per-screened** denominators (Figs. S1, S2)                                         | Reorders the per-man rate tables; the grade-mix ranking in §2.3 is unaffected, as it uses detected cancers as denominator                  |


**Publication bias was not assessed.** The protocol restricts Peters' regression and
funnel plots to analyses with ≥10 studies; the largest here has k = 9. Small-study
effects cannot be excluded — relevant given that the two trials opposing the family-B
GG1 reduction (MVP, PROSA) are also the two smallest.

**Protocol deviation:** §8 specifies a 0.5 continuity correction for zero cells in both
the risk-ratio and proportion analyses. It was applied to the risk-ratio models only.
The GLMM proportion models use a logit-normal likelihood that accommodates zero counts
directly, and adding 0.5 would bias those estimates; this is standard practice.

---



## 5. Interpretation

1. **Invitation to screening increases low-grade detection without a matching increase
  in significant cancer** (GG1 RR 1.76 vs csPCa RR 1.16), and shifts the detected-cancer
   mix toward GG1 (RR 1.21 excluding PROBASE). This is overdiagnosis quantified.
2. **Modern triage reverses the grade mix.** GG1 fell from 65.7% of detected cancers
  under PSA + DRE and 49.4% under PSA to 31.7% (MRI-first) and 16.4%
   (biomarker → MRI), p < 0.0001 between strategies; within-trial RR 0.57 in family B.
   The two biomarker → MRI arms agree closely (I² = 0%).
3. **The cost side differs by strategy type.** Biomarker-triaged pathways reduce GG1
  detection *and* biopsies (STHLM3-MRI, GÖTEBORG-2: biopsy RR 0.48 and 0.43).
   MRI-first pathways achieve the lowest GG1 share but at the highest biopsy rate
   (10.4 per 100 screened) — though on only 615 men across two small, early-terminated
   or single-round trials.
4. **The main threat to these conclusions is between-trial comparability, not sampling
  error.** Trials differ by 9–10 years vs one round of follow-up, by age 45 vs 55–69,
   and in control-arm contamination (PLCO's control arm had 40–52% PSA uptake, which
   biases its RRs toward 1). Per-man detection rates are therefore not interpretable
   across trials; the grade-mix outcome is, and it carries the review's message.

---



## 6. Output files


| File                                                                     | Contents                                                                                        |
| ------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------- |
| [output/results.txt](output/results.txt)                                 | Complete numeric output, all models, all leave-one-out runs                                     |
| [output/table1_descriptives.csv](output/table1_descriptives.csv)         | Arm-level extracted data with derived rates                                                     |
| [output/table2_pooled_estimates.csv](output/table2_pooled_estimates.csv) | All pooled estimates with Q, df, p, I²                                                          |
| `output/fig1`–`fig5` (pdf + png)                                         | Risk-ratio forest plots: GG1, csPCa, any PCa, grade mix, biopsies                               |
| `output/fig6`–`fig9` (pdf + png)                                         | Proportion forest plots by strategy: GG1/screened, csPCa/screened, GG1 share, biopsies/screened |
| `output/figS1`–`figS2` (pdf + png)                                       | Sensitivity: proportions on ITT denominators                                                    |
| [arm_level_data.csv](arm_level_data.csv)                                 | Analysis dataset, with per-arm provenance and caveat flags                                      |


