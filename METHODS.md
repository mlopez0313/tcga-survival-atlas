# Methods: TCGA Multi-Modal Survival Prediction

This document gives a deeper technical description of the models, losses, feature engineering, and design choices used in this repository. It has been revised against standard references on Cox regression, tied-event handling, random survival forests, and deep survival analysis, including Cox (1972), Breslow (1974), Efron (1977), Faraggi and Simon (1995), Ishwaran et al. (2008), Katzman et al. (2018), and recent deep-survival reviews [1–8].

---

## 1. Problem setup

We model overall survival from TCGA patient-level multi-modal data. For each patient \(i\), we observe:

- feature vector(s) \(x_i\)
- event or censoring time \(t_i\)
- event indicator \(\delta_i \in \{0,1\}\), where 1 means death observed and 0 means right-censored

The goal is to learn a risk function
\[
r(x_i) \in \mathbb{R}
\]
that orders patients so higher predicted risk corresponds to shorter survival.

This repository focuses on **relative risk prediction**, not direct calibration of full survival distributions. Models are compared primarily by the **concordance index (C-index)** and by **Kaplan-Meier separation of predicted risk groups** [1,8].

---

## 2. Survival analysis background

### 2.1 Hazard and survival functions

Let \(T\) denote the event time.

- Survival function:
\[
S(t) = P(T > t)
\]
- Density:
\[
f(t) = -\frac{d}{dt}S(t)
\]
- Hazard:
\[
h(t) = \lim_{\Delta t \to 0} \frac{P(t \le T < t+\Delta t \mid T \ge t)}{\Delta t}
\]

The hazard can be written as
\[
h(t) = \frac{f(t)}{S(t)}
\]
and therefore
\[
S(t) = \exp\left(-\int_0^t h(u)\,du\right).
\]

Equivalently, if \(H(t)=\int_0^t h(u)\,du\) is the cumulative hazard, then
\[
S(t)=\exp(-H(t)).
\]

In right-censored survival data, not all deaths are observed by the end of follow-up, so learning methods must use partial information correctly [1,8].

### 2.2 Cox proportional hazards model

The Cox model assumes
\[
h(t \mid x_i) = h_0(t) \exp(\eta_i)
\]
where:
- \(h_0(t)\) is an unspecified baseline hazard
- \(\eta_i\) is the linear predictor or risk score

For classical Cox regression,
\[
\eta_i = x_i^\top \beta.
\]

For neural survival models in this repo,
\[
\eta_i = f_\theta(x_i)
\]
where \(f_\theta\) is a neural network.

The proportional hazards assumption means covariates scale the hazard multiplicatively and their effect is time-invariant [1]. This is central not only to classical Cox regression but also to DeepSurv-style neural Cox models [6].

---

## 3. Partial likelihood and survival loss

### 3.1 Cox partial likelihood

Suppose event times are ordered. For each patient \(i\) with an observed event (\(\delta_i=1\)), define the risk set
\[
R_i = \{j : t_j \ge t_i\}
\]
as the set of patients still under observation just before time \(t_i\).

Under the Cox model, the probability that patient \(i\) is the one who fails at time \(t_i\), conditional on one event occurring in risk set \(R_i\), is
\[
P(i \text{ fails at } t_i \mid R_i) = \frac{\exp(\eta_i)}{\sum_{j \in R_i} \exp(\eta_j)}.
\]

Multiplying over observed events gives the partial likelihood:
\[
L(\eta) = \prod_{i: \delta_i = 1} \frac{\exp(\eta_i)}{\sum_{j \in R_i} \exp(\eta_j)}.
\]

Taking logs:
\[
\ell(\eta) = \sum_{i: \delta_i = 1} \left(\eta_i - \log \sum_{j \in R_i} \exp(\eta_j)\right).
\]

The optimization target is the negative partial log-likelihood:
\[
\mathcal{L}_{\text{Cox}} = -\ell(\eta).
\]

This is the core loss family used by:
- Cox elastic-net
- DeepSurv
- MultiBranch

By contrast, Random Survival Forest does **not** optimize Cox partial likelihood; it is a tree-ensemble survival method based on recursive partitioning and ensemble survival estimation [5].

### 3.2 Handling tied event times

When multiple patients share the same observed event time, the Cox partial likelihood needs a tie-handling convention [2,3].

Classically used options are:
- **Breslow** approximation [2]
- **Efron** approximation [3]
- **exact** partial likelihood

Breslow is computationally simple and widely used. Efron is generally more accurate when many ties are present, while the exact approach is most appropriate for heavily discrete event times but is more expensive [2,3].

In this repository:
- the **Cox elastic-net** model is implemented with `sksurv.linear_model.CoxnetSurvivalAnalysis`
- the **DeepSurv** and **MultiBranch** models use a custom Cox-style loss in `src/modeling/_torch_surv.py`

Because `CoxnetSurvivalAnalysis` is not the same estimator as `CoxPHSurvivalAnalysis`, it is not best to claim a uniform explicit Breslow implementation across all models. More accurately, the repository uses **practical Cox-style partial-likelihood implementations with approximate tie handling**. The neural models use a sorted risk-set loss without an explicit Efron correction.

Code-level note: the shared neural loss is documented in code as “Negative Cox partial log-likelihood with the Breslow approximation” and is implemented via descending-time sorting plus `torch.logcumsumexp`, normalized by the number of observed events.

### 3.3 Why optimize relative risk instead of survival time directly?

TCGA survival labels are censored and often heterogeneous across cohorts. Direct regression on observed survival time would incorrectly treat censored times as exact outcomes. Cox-style ranking avoids that problem by only requiring correct ordering of risk among comparable patients [1,6,8].

### 3.4 Survival-loss details in prose

At every observed death time, the model compares the patient who died against everyone still at risk at that moment. The model is rewarded when it assigns the deceased patient a higher risk score than the others who were still under observation.

So rather than predicting “this patient dies in 14.2 months,” the model learns statements like:
- patient A should be riskier than patient B
- patient C should be riskier than the patients still alive at C’s event time

This ranking-style learning is why Cox-based methods are natural for censored survival data. The neural models in this repo do not change that principle; they only change how the risk score is computed [1,6].

---

## 4. Evaluation metrics

### 4.1 Concordance index

The C-index measures whether higher predicted risk corresponds to earlier observed events among comparable pairs. This is the standard discrimination metric used in classical survival modeling and in modern deep-survival benchmarks, but it measures **ranking**, not calibration [8].

For comparable patient pairs \((i,j)\), concordance is:
- 1 if the patient with shorter survival has higher predicted risk
- 0.5 for a tie in risk
- 0 otherwise

The empirical C-index is the fraction of concordant comparable pairs.

Interpretation:
- 0.5: random ordering
- 1.0: perfect ordering
- below 0.5: systematically reversed ordering

Code-level note: this repository computes Harrell’s C-index with `sksurv.metrics.concordance_index_censored` in `src/evaluation/metrics.py`.

### 4.2 Kaplan-Meier separation by predicted risk

For visualization, patients are split into high- and low-risk groups by median predicted risk. Kaplan-Meier curves are then plotted for the two groups and compared with a log-rank test.

This does **not** train the model; it is only a downstream way to assess whether the learned risk score produces clinically meaningful separation.

Code-level note: the grouping threshold is the median risk by default (`q=0.5`), implemented in `risk_to_groups()` and used by `plot_km_by_risk()` in `src/evaluation/survival_curves.py`.

### 4.3 What is not currently emphasized

Modern survival benchmarking often also reports:
- time-dependent AUC
- integrated Brier score
- calibration of predicted survival curves

This repository currently emphasizes discrimination-oriented benchmarking rather than full survival-probability calibration.

---

## 5. Data modalities and feature engineering

The pipeline supports six main modalities:
- clinical
- RNA-seq
- somatic mutation
- copy-number variation (CNV)
- DNA methylation
- miRNA

All preprocessing is driven by `config.yaml`.

### 5.1 Clinical features

Clinical preprocessing aims to construct a compact patient-level covariate table and the survival target.

Inputs include candidate fields for:
- patient ID
- age
- sex
- stage
- smoking

Because TCGA/Xena exports are not perfectly standardized, the code searches over candidate column names and uses the first available match. This improves portability across cohorts and Xena schema variations.

Feature engineering considerations:
- **Age** is retained as a numeric covariate.
- **Sex**, **stage**, and **smoking-related variables** are typically encoded categorically.
- Rows missing the survival target are dropped if `drop_missing_target: true`.

Rationale:
- Clinical variables are low-dimensional and usually provide a strong baseline signal.
- They also guarantee the multi-branch model has at least one branch with broad cohort coverage.

### 5.2 RNA-seq

RNA-seq is high-dimensional, noisy, and often strongly right-skewed. In the TCGA/Xena setting used here, the repository expects a precomputed expression matrix such as FPKM rather than raw sequencing counts, so the preprocessing is best understood as a pragmatic feature-preparation step for predictive modeling rather than a full RNA-seq normalization workflow.

Configured steps:
- optional log transform: \(x \mapsto \log_2(x+1)\)
- remove genes with very low mean expression
- retain top variable genes

Current default parameters:
- `min_expression: 1.0`
- `top_variable_genes: 2000`

Rationale:
- Log transform stabilizes variance and reduces domination by very large values.
- Low-expression filtering removes weak, noisy genes.
- Variance filtering reduces dimensionality without using survival labels, limiting overfitting risk while keeping the procedure unsupervised.

Important caveat:
- For raw count RNA-seq data, one would normally discuss library-size normalization and count-aware models more explicitly.
- Here, because the default source is Xena-hosted processed expression tables, the workflow assumes the downloaded matrix is already in a usable expression scale and applies only lightweight downstream transformations.

Tradeoff:
- Top-variance filtering may drop genes with prognostic but low-variance behavior.
- However, for a general-purpose benchmark it is a simple and robust first-pass reduction.

### 5.3 Somatic mutation

Mutation data are converted into a patient-by-gene feature matrix, typically binary.

Configured steps:
- optionally require nonsynonymous mutations only
- keep top \(N\) most frequently mutated genes

Current defaults:
- `top_n_genes: 300`
- `require_nonsynonymous: true`

Rationale:
- Mutation matrices are sparse; most genes are never mutated in most patients.
- Restricting to the most frequently altered genes creates a denser and more learnable matrix.
- Nonsynonymous filtering prioritizes likely functional alterations over silent noise.

Tradeoff:
- Rare driver mutations can be discarded.
- A binary encoding ignores allele fraction, clonality, and mutation type severity.

### 5.4 Copy-number variation (CNV)

CNV features come from GISTIC-style values and can be treated either continuously or binarized.

Current default:
- `cnv_mode: continuous`
- `top_variable_features: 2000`

Optional binary setting:
- mark a feature altered if \(|x| > \tau\) for threshold \(\tau\)

Rationale for continuous default:
- Continuous GISTIC-like scores preserve more information about amplification/deletion magnitude.
- Variance filtering again gives an unsupervised dimensionality reduction step.

Rationale for binary option:
- more interpretable altered/not-altered encoding
- potentially more robust if continuous scores are noisy or platform-specific

### 5.5 DNA methylation

Methylation arrays are very high-dimensional and frequently contain missing values.

Configured steps:
- drop probes with too much missingness
- keep top variable probes

Current defaults:
- `max_missing_frac: 0.2`
- `top_variable_probes: 5000`

Rationale:
- Excessively missing probes are unreliable and complicate downstream modeling.
- Variance filtering is a scalable unsupervised compression strategy.
- Keeping more probes than RNA by default can be reasonable because methylation often contains distributed prognostic signal across many loci.

Tradeoff:
- Probe-level features are difficult to interpret biologically without later aggregation to genes/pathways/regions.

### 5.6 miRNA

Configured steps:
- optional log transform
- retain top variable miRNAs

Current defaults:
- `log_transform: true`
- `top_variable: 200`

Rationale:
- miRNA panels are lower-dimensional than gene expression but still benefit from dynamic-range compression and variance-based pruning.

---

## 6. Data integration strategy

There are two integration regimes in this repo.

### 6.1 Early fusion for baseline and DeepSurv models

For Cox, RSF, and DeepSurv, selected modalities are concatenated into one matrix:
\[
X_{\text{concat}} = [X_{\text{clinical}} \; | \; X_{\text{rna}} \; | \; X_{\text{mutation}} \; | \; \cdots]
\]

Advantages:
- simple
- easy to benchmark
- works with classical models directly

Disadvantages:
- ignores modality-specific structure
- a very large modality can dominate optimization
- does not explicitly learn per-modality embeddings

Code-level note: feature sets are loaded from processed parquet files and resolved by name via `select_feature_set()` in `src/modeling/_data.py`.

### 6.2 Intermediate fusion for MultiBranch

The multi-branch model creates one encoder per modality:
\[
z^{(m)}_i = g^{(m)}_{\theta_m}(x^{(m)}_i)
\]
for modality \(m\). The embeddings are concatenated:
\[
z_i = [z^{(1)}_i | z^{(2)}_i | \cdots | z^{(M)}_i]
\]
and passed through a fusion network:
\[
\eta_i = h_\phi(z_i).
\]

Advantages:
- each modality gets its own representation subnetwork
- allows different input dimensionalities naturally
- more appropriate for heterogeneous data types

Disadvantages:
- more parameters
- higher overfitting risk on modest TCGA sample sizes
- still assumes all selected modalities can be aligned at patient level

---

## 7. Model details

## 7.1 Cox elastic-net

### Formulation

The model is
\[
\eta_i = x_i^\top \beta
\]
with objective
\[
\min_{\beta} -\ell(\beta) + \lambda \left( \alpha \lVert \beta \rVert_1 + \frac{1-\alpha}{2} \lVert \beta \rVert_2^2 \right)
\]
where:
- \(\ell(\beta)\) is the Cox partial log-likelihood
- \(\lambda\) is regularization strength
- \(\alpha\) is the elastic-net mixing parameter (`l1_ratio`)

This is the standard regularized Cox setup used in high-dimensional survival modeling, especially when the number of molecular features is large relative to the number of patients.

### Why this model?

Cox elastic-net is a strong tabular survival baseline because:
- it handles censoring correctly
- it is interpretable via coefficients
- regularization is important in high dimensions
- it often performs surprisingly well relative to deep models on TCGA-scale cohorts

### Hyperparameter rationale

Current defaults:
- `l1_ratio: 0.5`
- `n_alphas: 50`
- `cv_folds: 5`

Reasoning:
- `l1_ratio=0.5` balances sparsity and stability.
- Searching a path of 50 alphas gives a practical compromise between coverage and speed.
- Five-fold CV is conventional and stable enough for cohort sizes typical of TCGA.

Code-level note: alpha selection is done by a custom K-fold procedure `_cv_pick_alpha()` in `src/modeling/cox_elastic_net.py`, then the final model is fit with a single chosen alpha.

### Interpretability

Nonzero coefficients indicate selected features.
- positive coefficient -> higher predicted risk
- negative coefficient -> lower predicted risk

The repository saves top coefficients for downstream review.

---

## 7.2 Random Survival Forest (RSF)

### Concept

RSF, introduced by Ishwaran et al. (2008), extends Breiman-style random forests to right-censored survival data. Each tree is grown from a randomized sample, candidate splits are evaluated with a survival-aware criterion such as log-rank separation, and terminal nodes store survival information derived from the training subjects that land there.

Predictions can be expressed in several related ways depending on implementation, including estimated cumulative hazard functions, survival functions, or scalar risk summaries. In this repository, the scikit-survival implementation is used and the downstream benchmarking reduces predictions to a risk score for comparison by C-index and Kaplan-Meier separation.

### Why this model?

RSF is a useful nonlinear baseline because it:
- captures interactions without manual feature engineering
- makes fewer linearity assumptions than Cox
- is often robust on mixed clinical/molecular feature sets

### Hyperparameter rationale

Current defaults:
- `n_estimators: 300`
- `min_samples_split: 10`
- `min_samples_leaf: 15`
- `max_features: sqrt`

Reasoning:
- 300 trees gives a reasonably stable ensemble without being excessive.
- A larger minimum leaf size regularizes the forest, important for moderate sample sizes and high-dimensional inputs.
- `sqrt` feature subsampling is a standard variance-reduction choice.

### Feature importance

The repository computes permutation importance on the test set, optionally capped to a subset of features for runtime reasons.

Caution:
- importance in correlated omics features can be unstable
- RSF importance should be treated as heuristic, not mechanistic proof

Code-level note: permutation importance is computed against the model’s `score()` output in `src/modeling/random_survival_forest.py`, and for speed may be restricted to the top-variance test features.

---

## 7.3 DeepSurv

### Context in the literature

DeepSurv, popularized by Katzman et al. (2018), is a neural-network generalization of the Cox proportional hazards model. Conceptually it extends earlier neural survival work such as Faraggi and Simon (1995), replacing the linear Cox predictor with a deeper nonlinear function approximator [4,6,7].

### Architecture

DeepSurv replaces the linear Cox predictor with a multilayer perceptron:
\[
\eta_i = f_\theta(x_i)
\]
where each hidden layer applies
\[
h^{(l+1)} = \text{Dropout}(\text{ReLU}(\text{BN}(W^{(l)} h^{(l)} + b^{(l)})))
\]
when batch normalization is enabled.

Default architecture:
- hidden layers: `[128, 64]`
- dropout: `0.3`
- batch norm: `true`

ASCII diagram:

```text
Concatenated feature vector
          |
      Linear(d -> 128)
          |
      BatchNorm
          |
         ReLU
          |
       Dropout
          |
      Linear(128 -> 64)
          |
      BatchNorm
          |
         ReLU
          |
       Dropout
          |
       Linear(64 -> 1)
          |
      Risk score eta
```

### Loss

DeepSurv retains the Cox partial-likelihood objective but replaces the linear predictor with a nonlinear function approximator. So the model still optimizes a proportional-hazards-style ranking objective, but allows nonlinear covariate effects and interactions.

Code-level note: the implementation in `src/modeling/deepsurv.py` uses the shared `cox_ph_loss()` in `src/modeling/_torch_surv.py`.

### Hyperparameter rationale

Current defaults:
- `lr: 1e-3`
- `weight_decay: 1e-4`
- `epochs: 200`
- `patience: 25`
- `val_size: 0.2`
- `batch_size: 256`

Reasoning:
- `1e-3` is a standard stable Adam learning rate.
- `weight_decay=1e-4` gives mild regularization.
- `[128,64]` is expressive enough for nonlinear structure without being excessively large for TCGA sample sizes.
- `dropout=0.3` is a moderate anti-overfitting choice.
- early stopping avoids wasting epochs after validation loss saturates.

### Important caveat on minibatching

The exact Cox partial likelihood is naturally a full-risk-set objective. In the classical derivation, each event time compares the event subject against the entire contemporaneous risk set. Minibatching approximates that objective because each batch only sees a subset of the true risk sets.

In practice:
- if the training set is smaller than `batch_size`, the code uses full-batch training
- otherwise it uses minibatches, skipping batches with no events

This is practical, but it means optimization is an approximation to the global Cox objective on larger cohorts.

### Implementation correction

The current code selects the best checkpoint by **lowest validation loss**, not by highest validation C-index. Validation C-index is tracked and plotted, but the saved `best_state` in `src/modeling/deepsurv.py` is updated when `val_loss` improves.

---

## 7.4 MultiBranch survival model

### Motivation

Different modalities have different dimensionalities, scales, sparsity patterns, and biological meanings. This is a standard motivation in modern multimodal learning and is especially relevant in cancer studies that combine clinical and multiple omics modalities. A single shared encoder on concatenated features may be suboptimal.

The multi-branch model therefore learns one encoder per modality and then fuses modality embeddings.

### Mathematical form

For each modality \(m\):
\[
z_i^{(m)} = g^{(m)}_{\theta_m}(x_i^{(m)}) \in \mathbb{R}^{k}
\]
where \(k\) is the branch embedding size.

The fused representation is
\[
z_i = \text{concat}(z_i^{(1)}, z_i^{(2)}, \dots, z_i^{(M)}).
\]

Then
\[
\eta_i = h_\phi(z_i)
\]
and training minimizes
\[
\mathcal{L}_{\text{Cox}}(\eta).
\]

### Default architecture

Branch settings:
- `branch_hidden: [64, 32]`
- `branch_embedding: 16`

Fusion settings:
- `fusion_hidden: [64, 32]`

ASCII diagram:

```text
Clinical features  -> Branch MLP -> z_clinical --\
RNA features       -> Branch MLP -> z_rna       -- \
Mutation features  -> Branch MLP -> z_mut       ----> Concatenate -> Fusion MLP -> Risk score eta
CNV features       -> Branch MLP -> z_cnv       -- /
Methylation        -> Branch MLP -> z_meth      --/
miRNA features     -> Branch MLP -> z_mirna    -/
```

More explicitly for one branch:

```text
Modality input x^(m)
      |
 Linear(d_m -> 64)
      |
 BatchNorm + ReLU + Dropout
      |
 Linear(64 -> 32)
      |
 BatchNorm + ReLU + Dropout
      |
 Linear(32 -> 16)
      |
     ReLU
      |
 embedding z^(m)
```

### Hyperparameter rationale

Current defaults:
- `branch_hidden: [64, 32]`
- `branch_embedding: 16`
- `fusion_hidden: [64, 32]`
- `dropout: 0.3`
- `lr: 1e-3`
- `weight_decay: 1e-4`
- `epochs: 300`
- `patience: 30`

Reasoning:
- Small branch networks reduce overfitting risk while still allowing modality-specific nonlinear compression.
- Embedding size 16 is a deliberate bottleneck to force compact modality summaries.
- Fusion network size mirrors branch size for simplicity and symmetry.
- Slightly longer max training (`300`) and patience (`30`) are reasonable because the architecture is deeper and optimization may take longer than DeepSurv.

### Training design choice

Per-branch standardization is fit on inner-train data only, then applied to validation/test. This prevents leakage.

The best model checkpoint is selected by **validation C-index**, while early stopping monitors **validation loss** as a smoother optimization signal.

This is a sensible compromise:
- C-index is the actual downstream target metric
- validation loss is often less noisy epoch-to-epoch

Code-level note: this description matches `src/modeling/multibranch_survival.py`. In contrast to DeepSurv, `best_state` here is updated when `val_cindex` improves, while early stopping still uses validation loss.

### Missing modalities

The implementation can skip requested modalities if their processed files are absent. However, it requires at least the clinical branch.

This is useful because real-world multi-omics cohorts are often incomplete.

### Important scope note

This multi-branch model is a **project-specific architecture**, not a direct implementation of one single canonical named method from the literature. It is best described as a multimodal neural Cox model inspired by standard multimodal deep-learning design patterns.

---

## 8. Why these models as a progression?

The repository uses a deliberate progression:

1. **Cox elastic-net**: interpretable linear baseline
2. **RSF**: nonlinear tree baseline
3. **DeepSurv**: nonlinear neural baseline on early-fused features
4. **MultiBranch**: modality-aware neural architecture

This progression is methodologically useful because it tests whether extra modeling complexity actually helps.

A key principle of the repo is that the deep models are not assumed to win. If Cox or RSF performs best on held-out data, that is the correct conclusion.

---

## 9. Train/validation/test handling and leakage control

### Global split

The main train/test split is fixed once and reused across all models.

Why this matters:
- model comparisons are fair
- performance differences are not confounded by different sample splits

### Inner validation

DeepSurv and MultiBranch carve a validation set out of training data only.

### Scaling and preprocessing

The code follows this general rule:
- unsupervised modality filtering may occur before model fitting
- any learned scaling used by the model is fit on train only

This is an important leakage control measure.

Code-level note: for early-fusion models, train-only scaling is handled by `fit_scaler_train()` in `src/modeling/_data.py`; for the multi-branch model, a separate `StandardScaler` is fit per modality on inner-train patients only.

One subtlety noted in the README: top-variance prefiltering is computed without survival labels, so it is unsupervised, but if computed on the full dataset it can still slightly couple train and test distributions. This is less severe than supervised leakage, but worth documenting explicitly.

---

## 10. Known assumptions and limitations

### 10.1 Proportional hazards assumption

Cox-based methods assume the effect of covariates is multiplicative on hazard and stable over time. This applies directly to Cox elastic-net and, in modeling spirit, also to DeepSurv and the MultiBranch model because both optimize a Cox-style proportional hazards objective. If strong time-varying effects exist, performance or interpretability may suffer.

### 10.2 Small-n / large-p setting

TCGA often has far fewer patients than molecular features. Even with filtering, overfitting is a major risk. This is why regularized classical baselines are essential.

### 10.3 Early fusion limitations

Concatenation can overweight large modalities and ignore structured modality relationships.

### 10.4 Missingness and cohort heterogeneity

Different cohorts and Xena versions may differ in completeness and annotation conventions. The clinical preprocessing tries to be robust, but portability is never perfect.

A concrete example from the GDC-backed cohort sweep is `TCGA-LAML`, where `TCGAbiolinks` currently fails in cohort-specific clinical metadata handling because expected fields such as `disease_response` are absent in some internal preparation paths. That failure occurs before modeling and should be treated as a data-access / downloader compatibility issue rather than evidence about model performance on LAML.

### 10.5 Tie handling and minibatch approximation

The Cox-style losses use practical approximations for ties, and neural training may approximate full-risk-set optimization when minibatching is used.

### 10.6 Discrimination vs calibration

A model can have a good C-index while still producing poorly calibrated survival probabilities. Since this repository emphasizes ranking metrics, users should avoid overinterpreting the outputs as fully calibrated survival distributions.

---

## 11. Practical interpretation guide

If you run this pipeline, interpret outcomes in this order:

1. **Test C-index**: primary ranking metric
2. **Train vs test gap**: overfitting check
3. **KM separation/log-rank p**: clinical interpretability of risk groups
4. **Feature importance / coefficients**: exploratory biological interpretation

A stronger deep model should show improvement on **held-out test C-index**, not just on training curves.

---

## 12. Suggested future extensions

Potential future improvements include:
- nested cross-validation rather than one train/test split
- time-dependent AUC and integrated Brier score in addition to C-index
- Efron tie handling in custom neural losses
- explicit missing-modality masking
- pathway-level or gene-set aggregation for omics features
- autoencoder or variational pretraining per modality
- transformer or graph-based multimodal fusion
- calibration of survival curves rather than only ranking risk

---

## 13. Relation to published methods

This repository combines:
- a **penalized Cox model** in the tradition of regularized Cox regression for high-dimensional data
- a **Random Survival Forest** baseline following the RSF literature
- a **DeepSurv-style** neural Cox model
- a **custom multimodal multi-branch neural Cox model** inspired by multimodal deep learning practice

That last point matters: the multi-branch model in this repo is best described as a **project-specific architecture**, not as a direct implementation of a single standard published method.

---

## 14. References

1. Cox DR. Regression Models and Life-Tables. *Journal of the Royal Statistical Society: Series B (Methodological)*. 1972;34(2):187–220.
2. Breslow NE. Covariance Analysis of Censored Survival Data. *Biometrics*. 1974;30(1):89–99.
3. Efron B. The Efficiency of Cox’s Likelihood Function for Censored Data. *Journal of the American Statistical Association*. 1977;72(359):557–565.
4. Faraggi D, Simon R. A Neural Network Model for Survival Data. *Statistics in Medicine*. 1995;14(1):73–82.
5. Ishwaran H, Kogalur UB, Blackstone EH, Lauer MS. Random Survival Forests. *The Annals of Applied Statistics*. 2008;2(3):841–860.
6. Katzman JL, Shaham U, Cloninger A, Bates J, Jiang T, Kluger Y. DeepSurv: Personalized Treatment Recommender System Using a Cox Proportional Hazards Deep Neural Network. *BMC Medical Research Methodology*. 2018;18:24.
7. Wiegrebe S, Kopper P, Sonabend R, Bischl B, Rahnenführer J. Deep Learning for Survival Analysis: A Review. *Artificial Intelligence Review*. 2024;57:65.
8. Pölsterl S. scikit-survival documentation. Cox and survival-model implementation notes, accessed 2026.

---

## 15. File map for methods-relevant code

Primary implementation files:
- `src/modeling/cox_elastic_net.py`
- `src/modeling/random_survival_forest.py`
- `src/modeling/deepsurv.py`
- `src/modeling/multibranch_survival.py`
- `src/modeling/_torch_surv.py`
- `src/modeling/_data.py`
- `src/evaluation/metrics.py`
- `src/evaluation/survival_curves.py`
- `src/preprocessing/`

Configuration source:
- `config.yaml`

Top-level overview:
- `README.md`

---

## 16. Executive summary

This repository is a benchmark-oriented survival modeling pipeline for TCGA multi-modal data. It uses a common evaluation framework across:
- penalized Cox regression
- random survival forest
- DeepSurv
- a custom multi-branch deep Cox model

The central modeling objective is the Cox partial likelihood, which learns patient risk rankings under censoring. Feature engineering is modality-specific but intentionally simple, mostly based on log transforms, missingness filtering, and unsupervised variance/frequency selection. The multi-branch model is the most biologically structured architecture in the repository, but the project is explicitly designed so that simpler baselines remain the standards to beat.
