#' Diagnose the indicators before combining them
#'
#' Summing indicators into a component assumes they measure one coherent
#' construct; combining components assumes they are not redundant. Neither
#' assumption is usually checked. This function checks both:
#'
#' * **internal consistency** of each component -- Cronbach's alpha (raw and
#'   standardised), the mean inter-item correlation (which, unlike alpha,
#'   does not rise with the number of indicators) and corrected item-total
#'   correlations of its oriented indicators. Alpha is only expected to be
#'   high for *reflective* components, whose indicators are interchangeable
#'   measures; for components declared `formative` in [as_ccva()] (distinct,
#'   complementary facets) a low alpha is expected and is not flagged;
#' * **cross-component correlation** of the component scores (Spearman);
#' * **redundancy** -- indicator pairs, within or across components, with
#'   absolute Spearman correlation above `redundancy`;
#' * **dimensionality** -- predictors per unit, flagging blocks where
#'   predictors approach or exceed the number of units (use a penalised
#'   method there).
#'
#' @param x A `ccva_data` object.
#' @param redundancy Absolute correlation above which a pair is flagged.
#' @param component_method Method used to score components for the
#'   cross-component correlations.
#' @return An object of class `ccva_check`.
#' @export
#' @examples
#' chk <- ccva_check(ccva_sim_data())
#' chk
ccva_check <- function(x, redundancy = 0.9, component_method = "mean") {
  .check_ccva(x)
  comps <- .active_components(x)
  n <- length(x$id)

  alpha <- do.call(rbind, lapply(comps, function(cmp) {
    X <- .oriented(x, cmp)
    k <- ncol(X)
    if (k < 2 || k > n) {
      return(data.frame(component = cmp, n_indicators = k, alpha = NA,
                        alpha_std = NA, mean_r = NA, p_over_n = k / n))
    }
    R <- stats::cor(X)
    data.frame(component = cmp, n_indicators = k,
               alpha = .cronbach(X), alpha_std = .cronbach(scale(X)),
               mean_r = mean(R[upper.tri(R)]), p_over_n = k / n)
  }))
  alpha$verdict <- cut(alpha$alpha, c(-Inf, 0.6, 0.7, 0.8, Inf),
                       labels = c("questionable", "acceptable-",
                                  "acceptable", "good"), right = FALSE)
  alpha$verdict <- as.character(alpha$verdict)
  alpha$verdict[alpha$component %in% x$formative] <-
    "formative: alpha not required"
  alpha$verdict[alpha$p_over_n >= 0.5] <- "high-dimensional: use penalised"

  item_total <- do.call(rbind, lapply(comps, function(cmp) {
    X <- .oriented(x, cmp)
    if (ncol(X) < 2 || ncol(X) > 200) return(NULL)
    tot <- rowSums(scale(X))
    data.frame(component = cmp, indicator = colnames(X),
               r_item_total = vapply(seq_len(ncol(X)), function(j)
                 stats::cor(X[, j], tot - scale(X)[, j]), numeric(1)),
               row.names = NULL)
  }))

  cs <- as.data.frame(lapply(stats::setNames(comps, comps), function(cmp)
    ccva_component(x, cmp, method = component_method)$score))
  comp_cor <- stats::cor(cs, method = "spearman")

  vars <- unlist(x$indicators[comps])
  red <- NULL
  if (length(vars) <= 500) {
    allX <- do.call(cbind, lapply(comps, function(cmp) .oriented(x, cmp)))
    R <- stats::cor(allX, method = "spearman")
    idx <- which(abs(R) > redundancy & upper.tri(R), arr.ind = TRUE)
    if (nrow(idx))
      red <- data.frame(var1 = rownames(R)[idx[, 1]],
                        var2 = colnames(R)[idx[, 2]],
                        r = R[idx],
                        comp1 = vapply(rownames(R)[idx[, 1]], .component_of,
                                       character(1), x = x),
                        comp2 = vapply(colnames(R)[idx[, 2]], .component_of,
                                       character(1), x = x),
                        row.names = NULL)
  }
  structure(list(alpha = alpha, item_total = item_total,
                 component_cor = comp_cor, redundant = red,
                 n = n, threshold = redundancy),
            class = "ccva_check")
}

#' @export
print.ccva_check <- function(x, ...) {
  cat("<ccva_check> internal consistency (Cronbach's alpha)\n")
  a <- x$alpha
  for (i in seq_len(nrow(a)))
    cat(sprintf("  %-12s %3d indicators  alpha %s  (std %s)  mean r %s  %s\n",
                a$component[i], a$n_indicators[i],
                if (is.na(a$alpha[i])) "  -  " else sprintf("%.3f", a$alpha[i]),
                if (is.na(a$alpha_std[i])) "  -  " else
                  sprintf("%.3f", a$alpha_std[i]),
                if (is.na(a$mean_r[i])) "  -  " else sprintf("%.2f", a$mean_r[i]),
                a$verdict[i]))
  cat("  guide (reflective components): >= 0.8 good, >= 0.7 acceptable, < 0.6 questionable;\n",
      " alpha rises with the number of indicators, so compare mean r too\n")
  cat("\ncross-component Spearman correlations\n")
  print(round(x$component_cor, 3))
  if (is.null(x$redundant)) {
    cat(sprintf("\nno indicator pairs with |r| > %.2f\n", x$threshold))
  } else {
    cat(sprintf("\n%d redundant indicator pair(s) (|r| > %.2f):\n",
                nrow(x$redundant), x$threshold))
    print(utils::head(x$redundant, 10), row.names = FALSE)
  }
  low <- x$item_total[!is.na(x$item_total$r_item_total) &
                        x$item_total$r_item_total < 0.2, ]
  if (!is.null(low) && nrow(low))
    cat("\nweak indicators (item-total r < 0.2):",
        paste(low$indicator, collapse = ", "), "\n")
  invisible(x)
}

.cronbach <- function(X) {
  k <- ncol(X)
  (k / (k - 1)) * (1 - sum(apply(X, 2, stats::var)) /
                     stats::var(rowSums(X)))
}

#' Validate vulnerability scores against observed outcomes
#'
#' A CCVA score that has never been compared with any observable outcome is
#' unvalidated. This reports the rank correlation of a fit's scores with
#' every outcome proxy (oriented so a positive correlation is good) and,
#' for simulated data, with the latent truth.
#'
#' Note: when the outcomes were also used to fit the model (latent,
#' shape and penalised methods), correlation with them measures fit, not
#' out-of-sample validity; use held-out outcomes or the truth for that.
#'
#' @param fit A `ccva_fit`.
#' @param truth Optional numeric vector of true vulnerability. If `NULL`,
#'   `fit$data$truth$V_star` is used when present.
#' @return A data frame with one row per target: Spearman and Kendall
#'   correlations.
#' @export
#' @examples
#' fit <- ccva_combine(ccva_sim_data(), "weighted")
#' ccva_validate(fit)
ccva_validate <- function(fit, truth = NULL) {
  x <- fit$data
  rows <- list()
  O <- .oriented_outcomes(x)
  if (!is.null(O)) for (o in colnames(O))
    rows[[o]] <- .val_row(o, "outcome", fit$score, O[, o])
  tv <- truth %||% x$truth$V_star
  if (!is.null(tv)) rows[["truth"]] <- .val_row("V* (truth)", "truth",
                                                fit$score, tv)
  if (!length(rows)) stop("No outcomes or truth to validate against.",
                          call. = FALSE)
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out
}

.val_row <- function(name, type, s, y) {
  data.frame(target = name, type = type,
             spearman = stats::cor(s, y, method = "spearman"),
             kendall = stats::cor(s, y, method = "kendall"))
}

#' Compare several CCVA fits
#'
#' @param ... Named `ccva_fit` objects (or a single named list of them),
#'   all fitted to the same units.
#' @param truth Optional true vulnerability; defaults to the simulated
#'   truth when present.
#' @return An object of class `ccva_comparison` with a `table` (validation
#'   correlations per method, class agreement with the first fit),
#'   `rank_cor` (Spearman correlations between methods' scores) and the
#'   fits.
#' @export
#' @examples
#' d <- ccva_sim_data()
#' cmp <- ccva_compare(
#'   ofori  = ccva_combine(d, "weighted", component_method = "ordinal"),
#'   latent = ccva_combine(d, "latent",
#'                shared_error = list(c("E1_risk", "E2_risk", "E3_risk"))))
#' cmp
ccva_compare <- function(..., truth = NULL) {
  fits <- list(...)
  if (length(fits) == 1 && is.list(fits[[1]]) &&
      !inherits(fits[[1]], "ccva_fit")) fits <- fits[[1]]
  if (is.null(names(fits)) || any(names(fits) == ""))
    stop("Give every fit a name, e.g. ccva_compare(a = fit1, b = fit2).",
         call. = FALSE)
  if (!all(vapply(fits, inherits, logical(1), "ccva_fit")))
    stop("All inputs must be ccva_fit objects.", call. = FALSE)
  ids <- names(fits[[1]]$score)
  if (!all(vapply(fits, function(f) identical(names(f$score), ids),
                  logical(1))))
    stop("All fits must cover the same units in the same order.",
         call. = FALSE)
  tbl <- do.call(rbind, lapply(names(fits), function(nm) {
    v <- ccva_validate(fits[[nm]], truth = truth)
    w <- stats::setNames(v$spearman, v$target)
    data.frame(method = nm, t(w), check.names = FALSE,
               class_agreement = mean(as.character(fits[[nm]]$class) ==
                                        as.character(fits[[1]]$class)),
               row.names = NULL)
  }))
  S <- sapply(fits, `[[`, "score")
  structure(list(table = tbl, rank_cor = stats::cor(S, method = "spearman"),
                 fits = fits, truth = truth %||% fits[[1]]$data$truth$V_star),
            class = "ccva_comparison")
}

#' @export
print.ccva_comparison <- function(x, ...) {
  cat("<ccva_comparison> Spearman correlation of scores with each target\n")
  tb <- x$table
  num <- vapply(tb, is.numeric, logical(1))
  tb[num] <- lapply(tb[num], round, 3)
  print(tb, row.names = FALSE)
  cat("\nrank correlation between methods\n")
  print(round(x$rank_cor, 3))
  invisible(x)
}
