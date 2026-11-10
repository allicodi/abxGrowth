# abxGrowth

> Doubly robust methods for estimating the effect of enteric infection on linear growth, mediated by antibiotic treatment

**Author:** Allison Codi

## Overview

`abxGrowth` is an R package for estimating how an enteric infection (e.g. *Shigella*-attributable diarrhea) affects linear growth (e.g. height-for-age *z*-score). The methods implemented in this package allow us to estimate the causal effect appropriately accounting for the role of antibiotic treatment as a mediator of Shigella and growth, as well as the role of disease severity as a confounder of antibiotics and growth and mediator of Shigella and growth. 

It accompanies the analysis repository [shigella_abx_growth_analysis](https://github.com/allicodi/shigella_abx_growth_analysis), which holds the code for the manuscript:

> *Effects of* Shigella *diarrhea on linear growth with and without antibiotic treatment: a meta-analysis of five multisite studies among children in low-resource settings.*

The package supports two comparison designs:

| Design | Comparison group | Estimands reported for each antibiotic level *a* |
|---|---|---|
| **No-etiology** (`case_control = FALSE`) | Diarrhea episodes with no identified etiology (or all other diarrhea) | E[Growth \| Infection = 1, Abx = *a*], E[Growth \| Infection = 0, Abx = *a*], and their difference |
| **Case-control** (`case_control = TRUE`) | Healthy (non-diarrheal) controls | E[Growth \| Case = 1, Abx = *a*], E[Growth \| Healthy control], and their difference |

Estimation uses augmented inverse probability weighting (AIPW) with [SuperLearner](https://cran.r-project.org/package=SuperLearner) for the nuisance models:

- outcome regression, with an optional second-stage regression that handles post-infection severity covariates
- antibiotic-treatment propensity
- infection propensity
- outcome-missingness propensity

Standard errors come from influence functions. Optional features include adjustment for co-infecting pathogen quantities, propensity-score truncation, and a marginal structural model (MSM) for effect heterogeneity (for example, by age).

## Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("allicodi/abxGrowth")
```

To reproduce the manuscript, install the version recorded in the analysis repository's `renv.lock`. Running `renv::restore()` from that repository does this for you; to install it by hand:

```r
devtools::install_github("allicodi/abxGrowth")
```

The package uses `SuperLearner`, `dplyr`, `fastDummies`, `table1`, `gtsummary`, `labelled`, `ggplot2`, and `splines`. Any SuperLearner learners you pass in (for example `SL.glmnet`, `SL.ranger`, `SL.earth`) need their own packages installed.

## Main functions

| Function | Purpose |
|---|---|
| `agaipw()` | Main estimator. Fits the AIPW estimator for either design and returns an object of class `agaipw_res` with a `print()` method. |
| `one_hot_encode()` | Prepares a dataset for `agaipw()`: one-hot encodes factor and character covariates, and returns the encoded data plus the matching covariate, severity, and site column names. |
| `impute_covariates()` | Mean or mode imputation of selected covariates, optionally within site, and adds an `<var>_imputed` indicator for each one. |
| `average_agaipw_across_seeds()` | Combines `agaipw()` results fitted with different seeds to reduce cross-validation randomness. |
| `aggcomp()`, `abx_growth_gcomp()`, `abx_growth_gcomp_case_control()` | Earlier g-computation estimators with bootstrap confidence intervals. These were not used for the primary manuscript results. |

Full argument documentation is available through `?agaipw` and the other help pages.

## Reproducing the manuscript

Every analysis in the manuscript calls `agaipw()` from the scripts in `code/run_analyses/` of [shigella_abx_growth_analysis](https://github.com/allicodi/shigella_abx_growth_analysis). That repository's README covers system requirements, the `renv` setup, the simulated demo data, and expected run times.

## License

MIT © Allison Codi. See [LICENSE](LICENSE).