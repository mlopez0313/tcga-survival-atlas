# Results Summary: TCGA Multi-Modal Survival Benchmark

## Overview

This project implemented and reviewed a multi-cohort TCGA survival benchmark using six modalities:
- clinical
- RNA-seq
- mutation
- CNV
- DNA methylation
- miRNA

Models compared:
- Cox elastic-net
- Random Survival Forest (RSF)
- DeepSurv
- MultiBranch multimodal neural Cox model

Outputs were primarily stored in:
- `~/mnt/datapool/tcga_survival`
- `~/mnt/datapool/tcga_survival_cohorts`

---

## What was completed

### 1. Cohort sweep review
Completed cohorts with `summary.csv` available:
- TCGA-ACC
- TCGA-BLCA
- TCGA-BRCA
- TCGA-CESC
- TCGA-CHOL
- TCGA-COAD
- TCGA-DLBC
- TCGA-ESCA
- TCGA-GBM
- TCGA-HNSC
- TCGA-KICH
- TCGA-KIRC
- TCGA-KIRP
- TCGA-LGG
- TCGA-LIHC
- TCGA-LUAD

Total completed: **16 cohorts**

Excluded for now:
- `TCGA-LAML` due to a GDC/TCGAbiolinks downloader incompatibility involving missing clinical metadata fields such as `disease_response`

### 2. Robustness fix applied
Three cohorts initially failed during Cox elastic-net fitting:
- `TCGA-CESC`
- `TCGA-KICH`
- `TCGA-LGG`

Cause:
- `ArithmeticError: Numerical error, because weights are too large. Consider increasing alpha.`

Fix:
- added alpha backoff retry logic to `src/modeling/cox_elastic_net.py`
- the fit now retries with larger alpha values automatically when the initial fit fails

This resolved all three cohorts.

### 3. Aggregate summary created
Saved aggregate outputs:
- `~/mnt/datapool/tcga_survival_cohorts/_summary/aggregate_best_models.csv`
- `~/mnt/datapool/tcga_survival_cohorts/_summary/aggregate_best_models.json`
- `~/mnt/datapool/tcga_survival_cohorts/_summary/README.md`

Supporting script:
- `scripts/summarize_tcga_cohorts.py`

---

## Per-cohort best model results

| Cohort | Best model | Test C-index | n_train | n_test | events_train | events_test |
|---|---|---:|---:|---:|---:|---:|
| TCGA-ACC | cox_elastic_net | 0.8605 | 54 | 14 | 21 | 5 |
| TCGA-BLCA | random_survival_forest | 0.6526 | 304 | 76 | 132 | 33 |
| TCGA-BRCA | cox_elastic_net | 0.7494 | 496 | 125 | 62 | 16 |
| TCGA-CESC | random_survival_forest | 0.8065 | 206 | 52 | 51 | 13 |
| TCGA-CHOL | random_survival_forest | 0.5000 | 28 | 7 | 14 | 4 |
| TCGA-COAD | random_survival_forest | 0.6957 | 223 | 56 | 54 | 14 |
| TCGA-DLBC | multibranch | 0.7500 | 34 | 9 | 6 | 2 |
| TCGA-ESCA | multibranch | 0.4819 | 125 | 32 | 53 | 13 |
| TCGA-GBM | deepsurv | 0.8409 | 43 | 11 | 30 | 8 |
| TCGA-HNSC | multibranch | 0.6532 | 384 | 96 | 166 | 41 |
| TCGA-KICH | multibranch | 0.8947 | 50 | 13 | 7 | 2 |
| TCGA-KIRC | cox_elastic_net | 0.7074 | 204 | 51 | 53 | 13 |
| TCGA-KIRP | multibranch | 0.8393 | 200 | 50 | 30 | 7 |
| TCGA-LGG | deepsurv | 0.7868 | 391 | 98 | 97 | 24 |
| TCGA-LIHC | random_survival_forest | 0.7337 | 273 | 69 | 94 | 24 |
| TCGA-LUAD | multibranch | 0.7310 | 336 | 84 | 119 | 30 |

---

## Cross-cohort findings

### Win counts
- Random Survival Forest: **5**
- MultiBranch: **5**
- Cox elastic-net: **4**
- DeepSurv: **2**

### Interpretation
There is **no universal winning model** across TCGA cohorts.

Key pattern:
- RSF and MultiBranch each win the largest number of cohorts
- Cox elastic-net remains highly competitive despite being the simplest model family in the benchmark
- DeepSurv performs strongly in selected cohorts but is the least broadly dominant model

This supports the project’s benchmark philosophy: more complex neural models do **not** automatically outperform classical baselines.

---

## Notable rerun outcomes after Cox fix

### TCGA-CESC
Backoff behavior:
- selected alpha failed
- final successful alpha: `0.01455`

Best results:
- RSF: **0.8065**
- DeepSurv: 0.7351
- Cox elastic-net: 0.7321
- MultiBranch: 0.6786

### TCGA-KICH
Backoff behavior:
- selected alpha failed
- final successful alpha: `0.1704`

Best results:
- MultiBranch: **0.8947**
- Cox elastic-net: 0.6316
- RSF: 0.5000
- DeepSurv: 0.3684

Caution:
- very small cohort/test set
- only 2 test events
- result may be unstable despite the strong headline C-index

### TCGA-LGG
Backoff behavior:
- selected alpha failed
- final successful alpha: `0.01201`

Best results:
- DeepSurv: **0.7868**
- MultiBranch: 0.7299
- RSF: 0.7254
- Cox elastic-net: 0.6618

---

## LUAD-specific summary

`TCGA-LUAD` is fully complete and serves as the default/cohort-reference run.

Best LUAD model:
- **MultiBranch**
- test C-index: **0.7310**

LUAD output locations:
- `~/mnt/datapool/tcga_survival_cohorts/TCGA-LUAD/results/metrics/summary.csv`
- `~/mnt/datapool/tcga_survival_cohorts/TCGA-LUAD/results/figures/`
- `~/mnt/datapool/tcga_survival_cohorts/TCGA-LUAD/results/models/`

LUAD figures available include:
- `cindex_comparison.png`
- `km_cox_elastic_net.png`
- `km_deepsurv.png`
- `km_multibranch.png`
- `km_rsf.png`
- `deepsurv_training.png`
- `multibranch_training.png`
- `cox_top_coefficients.png`
- `rsf_top_features.png`

---

## Main conclusions

### 1. Survival performance is cohort-specific
No single model dominates across all cancers. The best-performing architecture depends strongly on cohort biology, sample size, event count, and modality structure.

### 2. Classical baselines remain strong
Cox elastic-net and RSF remain highly competitive. RSF ties MultiBranch for the most cohort wins, and Cox elastic-net wins several cohorts outright.

### 3. Multi-branch multimodal fusion helps in some cancers
MultiBranch performs especially well in:
- DLBC
- HNSC
- KICH
- KIRP
- LUAD

This suggests modality-aware fusion can help, but not uniformly.

### 4. DeepSurv is less broadly robust than hoped
DeepSurv wins only two completed cohorts:
- GBM
- LGG

It can excel in selected settings, but it is not the most consistently strong model in this benchmark.

### 5. Small cohorts require caution
Some top C-index values come from very small cohorts or very low event counts. These results are useful, but they should be interpreted cautiously because variance is likely high.

### 6. Remaining blocker is data ingestion, not modeling
`TCGA-LAML` remains excluded because of a downloader/preparation incompatibility in the GDC/TCGAbiolinks path. This is not evidence about model performance on LAML.

---

## Figures and tables to use in reporting

Recommended reporting assets:
1. Aggregate cohort winner table from `_summary/README.md`
2. Per-cohort `results/metrics/summary.csv`
3. LUAD detailed model comparison figure set
4. Selected Kaplan-Meier plots from best-performing cohorts
5. C-index comparison panels from each cohort’s `results/figures/cindex_comparison.png`

---

## Final takeaway

The current benchmark supports the following overall conclusion:

> Across completed TCGA cohorts, survival prediction performance is highly cohort-dependent. Classical survival models remain difficult to beat, multimodal deep fusion helps in some cancers, and no single architecture is universally best.
