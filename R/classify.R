#' Classify vulnerability scores into Low / Moderate / High
#'
#' Conventional CCVAs cut scores at 1/3 and 2/3, which silently assumes the
#' scores are uniformly distributed and has no link to any observed outcome.
#' This function offers that convention and three alternatives.
#'
#' @param score Numeric vector of vulnerability scores.
#' @param cutoffs One of
#'   * `"uniform"` -- fixed cuts at 1/3 and 2/3 (Ofori et al. 2017);
#'   * `"quantile"` -- sample tertiles (equal-sized classes);
#'   * `"mixture"` -- natural clusters from a three-component Gaussian
#'     mixture fitted by EM;
#'   * `"roc"` -- outcome-calibrated cuts maximising Youden's J against an
#'     ordinal outcome (see `roc_outcome`, `roc_levels`);
#'   * a numeric vector of two cut points.
#' @param data A `ccva_data` object; needed for `"roc"`.
#' @param roc_outcome Outcome column used for `"roc"` (default: the first
#'   outcome). It is oriented so higher = more vulnerable.
#' @param roc_levels Two outcome levels on the oriented scale: units at or
#'   above the first are "at least Moderate", at or above the second are
#'   "High". Defaults to the second and third distinct values.
#' @param labels Class labels.
#' @return An ordered factor with attribute `cutoffs`.
#' @export
#' @examples
#' fit <- ccva_combine(ccva_sim_data(), "weighted")
#' table(ccva_classify(fit$score, "uniform"))
#' table(ccva_classify(fit$score, "roc", data = fit$data))
ccva_classify <- function(score, cutoffs = "uniform", data = NULL,
                          roc_outcome = NULL, roc_levels = NULL,
                          labels = c("Low", "Moderate", "High")) {
  br <- if (is.numeric(cutoffs)) {
    if (length(cutoffs) != 2) stop("Give two numeric cut points.",
                                   call. = FALSE)
    sort(cutoffs)
  } else switch(cutoffs,
    uniform = c(1 / 3, 2 / 3),
    quantile = stats::quantile(score, c(1 / 3, 2 / 3), names = FALSE),
    mixture = .cut_mixture(score),
    roc = .cut_roc(score, data, roc_outcome, roc_levels),
    stop("Unknown `cutoffs`: ", cutoffs, call. = FALSE))
  cls <- cut(score, c(-Inf, br, Inf), labels = labels, ordered_result = TRUE,
             right = FALSE)
  attr(cls, "cutoffs") <- stats::setNames(br, c("low_moderate",
                                                "moderate_high"))
  attr(cls, "method") <- if (is.numeric(cutoffs)) "manual" else cutoffs
  cls
}

.cut_mixture <- function(score, k = 3, iter = 500, tol = 1e-8) {
  # Three-component univariate Gaussian mixture with equal variances,
  # fitted by EM; cut points where the posterior class changes.
  mu <- stats::quantile(score, c(1, 3, 5) / 6, names = FALSE)
  s2 <- stats::var(score) / k
  pi_k <- rep(1 / k, k)
  ll_old <- -Inf
  for (it in seq_len(iter)) {
    dens <- vapply(seq_len(k), function(j)
      pi_k[j] * stats::dnorm(score, mu[j], sqrt(s2)), numeric(length(score)))
    tot <- rowSums(dens)
    ll <- sum(log(tot))
    z <- dens / tot
    nk <- colSums(z)
    pi_k <- nk / length(score)
    mu <- colSums(z * score) / nk
    s2 <- sum(z * outer(score, mu, "-")^2) / length(score)
    if (abs(ll - ll_old) < tol) break
    ll_old <- ll
  }
  o <- order(mu)
  cls <- match(max.col(z, ties.method = "first"), o)
  if (length(unique(cls)) < 3)
    stop("Mixture model found fewer than three groups; use another ",
         "`cutoffs` method.", call. = FALSE)
  c(min(score[cls >= 2]), min(score[cls == 3]))
}

.cut_roc <- function(score, data, outcome, levels) {
  if (is.null(data) || !length(data$outcomes))
    stop("cutoffs = \"roc\" needs `data` with outcomes.", call. = FALSE)
  outcome <- outcome %||% data$outcomes[1]
  y <- data$data[[outcome]] * data$direction[[outcome]]
  u <- sort(unique(y))
  if (length(u) < 3 && is.null(levels))
    stop("ROC cut-offs need an ordinal outcome with at least 3 levels.",
         call. = FALSE)
  levels <- levels %||% u[2:3]
  youden <- function(pos) {
    cand <- sort(unique(score))
    j <- vapply(cand, function(t) {
      pred <- score >= t
      mean(pred[pos]) + mean(!pred[!pos]) - 1
    }, numeric(1))
    cand[which.max(j)]
  }
  lo <- youden(y >= levels[1]); hi <- youden(y >= levels[2])
  if (hi <= lo) hi <- max(score[score > lo], lo)
  c(lo, hi)
}
