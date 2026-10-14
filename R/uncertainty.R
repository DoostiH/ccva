#' Quantify uncertainty in vulnerability scores and rankings
#'
#' Published CCVAs usually report a sorted list of point estimates, which
#' implies a precision the data rarely support. Two complementary sources of
#' uncertainty are available:
#'
#' * `type = "bootstrap"` -- **sampling uncertainty**: units are resampled
#'   with replacement, the *whole* pipeline (component scoring, constraints,
#'   combination model) is refitted on each resample, and the original units
#'   are scored with the refitted model. For the latent model each unit's
#'   score is a prediction of an unobserved quantity, so (with
#'   `score_error = TRUE`) its factor-score standard error is added to every
#'   replicate; without it, intervals reflect parameter uncertainty only and
#'   are far too narrow.
#' * `type = "weights"` -- **weighting uncertainty** for the weighted index
#'   (a partial-identification view): component weights are drawn uniformly
#'   from the simplex, so the intervals show how much a ranking depends on
#'   the arbitrary choice of weights.
#'
#' For each unit the result gives an interval for its score and for its
#' rank; for each pair of units, the probability that one is more vulnerable
#' than the other. A ranking of two units is called *robust* when that
#' probability is above `1 - robust` or below `robust`.
#'
#' @param fit A `ccva_fit`.
#' @param type `"bootstrap"` or `"weights"`.
#' @param R Number of replicates.
#' @param level Interval coverage.
#' @param robust Threshold for calling a pairwise ranking robust.
#' @param score_error For the latent model, add each unit's factor-score
#'   prediction error to the bootstrap draws (recommended).
#' @param seed Optional random seed.
#' @param progress Show a progress bar.
#' @return An object of class `ccva_uncertainty` with `summary` (one row
#'   per unit: estimate, interval, rank interval), `pairwise` (matrix of
#'   P(row unit more vulnerable than column unit)), `robust_share`,
#'   `overlap_share` and the raw `draws`.
#' @export
#' @examples
#' fit <- ccva_combine(ccva_sim_data(), "weighted")
#' u <- ccva_uncertainty(fit, "weights", R = 500, seed = 1)
#' u
#' \donttest{
#' ub <- ccva_uncertainty(fit, "bootstrap", R = 50, seed = 1)
#' head(ub$summary)
#' }
ccva_uncertainty <- function(fit, type = c("bootstrap", "weights"),
                             R = 200, level = 0.95, robust = 0.05,
                             score_error = TRUE, seed = NULL,
                             progress = interactive()) {
  type <- match.arg(type)
  if (!inherits(fit, "ccva_fit")) stop("`fit` must be a ccva_fit.",
                                       call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  x <- fit$data
  n <- length(x$id)
  failures <- 0L

  if (type == "weights") {
    if (fit$method != "weighted")
      stop("type = \"weights\" applies to the weighted index only.",
           call. = FALSE)
    C <- as.matrix(fit$component_scores)
    s <- fit$details$signs[colnames(C)]
    k <- ncol(C)
    W <- matrix(stats::rexp(R * k), R, k)
    W <- W / rowSums(W)
    draws <- t(apply(W, 1, function(w)
      (C %*% (w * s) + sum(w[s < 0])) / sum(w)))
  } else {
    if (progress) pb <- utils::txtProgressBar(0, R, style = 3)
    draws <- matrix(NA_real_, R, n)
    for (b in seq_len(R)) {
      rows <- sample.int(n, n, replace = TRUE)
      xb <- .subset_ccva(x, rows)
      p <- tryCatch(suppressMessages(suppressWarnings({
        fb <- do.call(ccva_combine, c(list(x = xb), fit$args))
        raw <- fb$predict_raw(x$data)
        se <- fb$details$score_se
        if (score_error && !is.null(se)) raw <- raw + stats::rnorm(n, 0, se)
        .rescale01(raw, fb$raw_range[1], fb$raw_range[2])
      })), error = function(e) NULL)
      if (is.null(p) || any(!is.finite(p))) {
        failures <- failures + 1L
      } else {
        draws[b, ] <- p
      }
      if (progress) utils::setTxtProgressBar(pb, b)
    }
    if (progress) close(pb)
    draws <- draws[stats::complete.cases(draws), , drop = FALSE]
    if (nrow(draws) < 10)
      stop("Too few successful bootstrap replicates (", nrow(draws), ").",
           call. = FALSE)
    if (failures)
      message(failures, " of ", R, " bootstrap replicates failed and were ",
              "dropped.")
  }
  colnames(draws) <- x$id

  a <- (1 - level) / 2
  ranks <- t(apply(-draws, 1, rank, ties.method = "min"))
  summ <- data.frame(
    unit = x$id, estimate = as.numeric(fit$score),
    lower = apply(draws, 2, stats::quantile, a, names = FALSE),
    upper = apply(draws, 2, stats::quantile, 1 - a, names = FALSE),
    rank = rank(-fit$score, ties.method = "min"),
    rank_low = apply(ranks, 2, stats::quantile, a, names = FALSE),
    rank_high = apply(ranks, 2, stats::quantile, 1 - a, names = FALSE),
    row.names = NULL)
  summ$width <- summ$upper - summ$lower

  # pairwise probabilities P(i more vulnerable than j)
  P <- matrix(0, n, n, dimnames = list(x$id, x$id))
  for (r in seq_len(nrow(draws))) P <- P + outer(draws[r, ], draws[r, ], ">")
  P <- P / nrow(draws)
  ut <- upper.tri(P)
  robust_share <- mean(P[ut] > 1 - robust | P[ut] < robust)
  ov <- outer(summ$lower, summ$upper, "<") & outer(summ$upper, summ$lower, ">")
  overlap_share <- mean(ov[ut])

  structure(list(summary = summ[order(summ$rank), ], pairwise = P,
                 robust_share = robust_share, overlap_share = overlap_share,
                 draws = draws, type = type, level = level, robust = robust,
                 R = nrow(draws), failures = failures, method = fit$method),
            class = "ccva_uncertainty")
}

#' @export
print.ccva_uncertainty <- function(x, ...) {
  lab <- if (x$type == "bootstrap") "bootstrap (sampling)" else
    "random weights (weighting)"
  cat(sprintf("<ccva_uncertainty> %s, %d replicates, %d%% intervals\n",
              lab, x$R, round(100 * x$level)))
  s <- x$summary
  cat(sprintf("  mean interval width: %.3f (score range %.3f)\n",
              mean(s$width), diff(range(s$estimate))))
  cat(sprintf("  pairs with overlapping intervals: %.1f%%\n",
              100 * x$overlap_share))
  cat(sprintf("  robust pairwise rankings (P > %.2f): %.1f%%\n",
              1 - x$robust, 100 * x$robust_share))
  top <- utils::head(s, 5)
  cat("  most vulnerable:\n")
  for (i in seq_len(nrow(top)))
    cat(sprintf("    %-16s %.3f [%.3f, %.3f]  rank %d (%g-%g)\n",
                top$unit[i], top$estimate[i], top$lower[i], top$upper[i],
                top$rank[i], round(top$rank_low[i]), round(top$rank_high[i])))
  invisible(x)
}
