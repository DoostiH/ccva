#' ccva: Climate Change Vulnerability Assessment with Constrained Latent Models
#'
#' A general framework for climate change vulnerability assessment (CCVA).
#' The workflow has four stages:
#'
#' 1. **Data** -- [as_ccva()] tags each variable of any dataset as an
#'    indicator of exposure, sensitivity or adaptive capacity, gives its
#'    direction, and names the observable outcomes (vulnerability proxies).
#'    [ccva_simulate()] generates data with a known latent structure.
#' 2. **Components** -- [ccva_exposure()], [ccva_sensitivity()] and
#'    [ccva_adaptive()] reduce each component's own predictors to a score
#'    (Ofori-style ordinal scoring, continuous mean, PCA or penalised
#'    regression for high-dimensional blocks such as genomic data).
#' 3. **Vulnerability** -- [ccva_combine()] combines the components with the
#'    weighted index of Ofori et al. (2017), a shape-constrained additive
#'    model, a latent-variable structural equation model, or a penalised
#'    two-stage model. Ecological knowledge enters through
#'    [ccva_constraints()]. [ccva()] runs the whole chain.
#' 4. **Assessment** -- [ccva_check()], [ccva_validate()], [ccva_compare()],
#'    [ccva_uncertainty()], [ccva_classify()], [ccva_rank()] and
#'    [ccva_plot()] diagnose, validate, quantify uncertainty and visualise
#'    every step.
#'
#' @references
#' Ofori, B. Y., Stow, A. J., Baumgartner, J. B. and Beaumont, L. J. (2017).
#' Influence of adaptive capacity on the outcome of climate change
#' vulnerability assessment. *Scientific Reports*, 7, 12979.
#'
#' Bollen, K. A. (1989). *Structural Equations with Latent Variables*. Wiley.
#'
#' Pya, N. and Wood, S. N. (2015). Shape constrained additive models.
#' *Statistics and Computing*, 25, 543--559.
#'
#' @keywords internal
"_PACKAGE"

## Shared constants ---------------------------------------------------------

.ccva_components <- c("exposure", "sensitivity", "adaptive")

.ccva_labels <- c(exposure = "Exposure", sensitivity = "Sensitivity",
                  adaptive = "Adaptive capacity",
                  vulnerability = "Vulnerability")

# Latent factor names used inside lavaan models
.ccva_latent <- c(exposure = "Expos", sensitivity = "Sens",
                  adaptive = "Adapt")

utils::globalVariables(c(
  "component", "indicator", "value", "score", "unit", "lower", "upper",
  "rank_low", "rank_high", "estimate", "method", "truth", "x", "y",
  "xend", "yend", "label", "from", "to", "rho", "outcome", "coef",
  "var1", "var2", "r", "xlab", "ylab", "class", "fit", "group",
  "fill_col", "lab_col", "rank_a", "rank_b", "id", "pos", "alpha", "lx",
  "ly", "p", "std", "probability", "predicted", "observed", "n",
  "share_high", "mean_vulnerability", "n_species", "lab", "site"
))
