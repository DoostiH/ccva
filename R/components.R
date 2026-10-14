#' Score the exposure, sensitivity and adaptive-capacity components
#'
#' Each component is computed from its own indicators (its predictors) by a
#' chosen method. Indicators are first oriented with the `direction` given in
#' [as_ccva()], so that a higher value always means more of the component.
#'
#' Methods:
#' \describe{
#'   \item{`"ordinal"`}{The Ofori et al. (2017) scheme: each indicator is
#'     scored 1--3 (by `breaks`, or by sample tertiles when `breaks` is
#'     `NULL`), and the scores are summed with optional `weights` (the
#'     maximum score of each indicator) and divided by the maximum. Scores
#'     therefore lie in \[1/3, 1\], exactly as in the original study. Use it
#'     to reproduce conventional assessments and to measure what
#'     discretisation loses.}
#'   \item{`"mean"`}{Continuous: each indicator rescaled to \[0, 1\], then a
#'     weighted mean. No information is discarded.}
#'   \item{`"pca"`}{First principal component of the standardised
#'     indicators, signed so loadings are positive on average, rescaled to
#'     \[0, 1\]. Data-driven weights.}
#'   \item{`"penalised"`}{Penalised regression (LASSO by default; set
#'     `alpha`) of an `anchor` on the indicators, with coefficient signs
#'     constrained to the ecological direction. Designed for
#'     high-dimensional blocks where predictors outnumber units, such as
#'     genomic markers. Needs the \pkg{glmnet} package.}
#' }
#' Methods added with [ccva_register()] are also accepted.
#'
#' @param x A `ccva_data` object.
#' @param component One of `"exposure"`, `"sensitivity"`, `"adaptive"`.
#' @param method Scoring method; see Details.
#' @param weights Optional named numeric vector of indicator weights (for
#'   `"ordinal"`: maximum scores, e.g. Ofori's 11, 10, 9, ... for
#'   sensitivity). Missing indicators get weight 1.
#' @param breaks For `"ordinal"`, an optional named list giving, per
#'   indicator, two cut points on the raw scale. Values above the upper cut
#'   score 3 (more of the component) when the direction is `+1`, and 1 when
#'   it is `-1`.
#' @param anchor For `"penalised"`: an outcome column name, a numeric vector
#'   (higher = more vulnerable), or `NULL` to use the first principal
#'   component of all outcomes.
#' @param alpha Elastic-net mixing for `"penalised"`: 1 = LASSO, 0 = ridge.
#' @param constraints A [ccva_constraints()] object. Predictor-level rules
#'   are honoured: fixed-at-zero indicators are dropped; signs and bounds
#'   restrict the penalised coefficients.
#' @param ... Passed to registered methods.
#'
#' @return An object of class `ccva_component` with elements `score`
#'   (named numeric vector), `importance` (per-indicator weight, loading or
#'   coefficient), `method`, `component` and `predict` (a function scoring
#'   new data on the same scale).
#' @seealso [ccva_combine()], which calls these for you.
#' @export
#' @examples
#' d <- ccva_sim_data()
#' e_ord  <- ccva_exposure(d, method = "ordinal")
#' e_cont <- ccva_exposure(d, method = "mean")
#' cor(e_ord$score, e_cont$score, method = "spearman")
#'
#' s <- ccva_sensitivity(d, method = "pca")
#' head(sort(s$importance, decreasing = TRUE))
ccva_component <- function(x, component = c("exposure", "sensitivity",
                                             "adaptive"),
                           method = "mean", weights = NULL, breaks = NULL,
                           anchor = NULL, alpha = 1, constraints = NULL,
                           ...) {
  .check_ccva(x)
  component <- match.arg(component)
  vars <- x$indicators[[component]]
  if (!length(vars))
    stop("No ", component, " indicators in the data.", call. = FALSE)
  cons <- .resolve_constraints(constraints, x)

  # drop indicators fixed at zero
  drop <- names(Filter(function(v) v == 0, cons$ind_fixed))
  vars <- setdiff(vars, drop)
  if (!length(vars))
    stop("All ", component, " indicators were fixed at zero.", call. = FALSE)

  w <- stats::setNames(rep(1, length(vars)), vars)
  if (!is.null(weights)) {
    unk <- setdiff(names(weights), x$indicators[[component]])
    if (length(unk)) stop("Unknown weight name(s): ",
                          paste(unk, collapse = ", "), call. = FALSE)
    w[intersect(names(weights), vars)] <-
      weights[intersect(names(weights), vars)]
    if (any(w < 0)) stop("`weights` must be non-negative.", call. = FALSE)
  }

  X <- .oriented(x, component)[, vars, drop = FALSE]
  raw <- as.matrix(x$data[vars])
  dir <- x$direction[vars]

  fit <- switch(method,
    ordinal = .comp_ordinal(raw, dir, w, breaks),
    mean = .comp_mean(X, w),
    pca = .comp_pca(X),
    penalised = .comp_penalised(X, .anchor(x, anchor), component, cons,
                                vars, alpha),
    {
      reg <- .ccva_registry$component[[method]]
      if (is.null(reg))
        stop("Unknown component method '", method, "'. See ccva_methods().",
             call. = FALSE)
      out <- reg$fun(X, anchor = .anchor(x, anchor, required = FALSE),
                     signs = cons$ind_sign[vars], ...)
      if (is.null(out$score) || is.null(out$predict))
        stop("Registered method must return `score` and `predict`.",
             call. = FALSE)
      out$importance <- out$importance %||% stats::setNames(
        rep(NA_real_, length(vars)), vars)
      out
    })

  orient_fun <- if (method == "ordinal") {
    function(newdata) as.matrix(newdata[vars])
  } else {
    function(newdata) sweep(as.matrix(newdata[vars]), 2, dir, `*`)
  }
  pred <- fit$predict
  structure(list(
    component = component, method = method, indicators = vars,
    dropped = drop,
    score = stats::setNames(as.numeric(fit$score), x$id),
    importance = fit$importance,
    details = fit$details,
    predict = function(newdata) as.numeric(pred(orient_fun(newdata)))
  ), class = "ccva_component")
}

#' @rdname ccva_component
#' @export
ccva_exposure <- function(x, method = "mean", ...)
  ccva_component(x, "exposure", method = method, ...)

#' @rdname ccva_component
#' @export
ccva_sensitivity <- function(x, method = "mean", ...)
  ccva_component(x, "sensitivity", method = method, ...)

#' @rdname ccva_component
#' @export
ccva_adaptive <- function(x, method = "mean", ...)
  ccva_component(x, "adaptive", method = method, ...)

#' @export
print.ccva_component <- function(x, ...) {
  cat(sprintf("<ccva_component> %s via '%s' from %d indicator(s)\n",
              .ccva_labels[[x$component]], x$method, length(x$indicators)))
  if (length(x$dropped))
    cat("  dropped (fixed at 0):", paste(x$dropped, collapse = ", "), "\n")
  cat(sprintf("  score: mean %.3f, range [%.3f, %.3f]\n", mean(x$score),
              min(x$score), max(x$score)))
  imp <- x$importance[order(-abs(x$importance))]
  imp <- imp[!is.na(imp)]
  if (length(imp))
    cat("  top indicators:",
        paste(sprintf("%s (%.2f)", names(imp)[seq_len(min(5, length(imp)))],
                      imp[seq_len(min(5, length(imp)))]), collapse = ", "),
        "\n")
  invisible(x)
}

## Method implementations -----------------------------------------------------

.comp_ordinal <- function(raw, dir, w, breaks) {
  vars <- colnames(raw)
  cuts <- lapply(vars, function(v) {
    if (!is.null(breaks[[v]])) {
      b <- sort(breaks[[v]])
      if (length(b) != 2) stop("`breaks$", v, "` needs two cut points.",
                               call. = FALSE)
      b
    } else {
      stats::quantile(raw[, v], c(1 / 3, 2 / 3), names = FALSE)
    }
  })
  names(cuts) <- vars
  score_items <- function(R) {
    vapply(vars, function(v) {
      it <- findInterval(R[, v], cuts[[v]], left.open = TRUE) + 1
      if (dir[[v]] < 0) 4 - it else it
    }, numeric(nrow(R)))
  }
  sc <- function(R) {
    items <- score_items(R)
    if (is.null(dim(items))) items <- matrix(items, nrow = 1)
    as.numeric((items / 3) %*% w / sum(w))
  }
  list(score = sc(raw), predict = sc,
       importance = w / sum(w),
       details = list(breaks = cuts, items = score_items(raw)))
}

.comp_mean <- function(X, w) {
  lo <- apply(X, 2, min); hi <- apply(X, 2, max)
  rng <- ifelse(hi - lo > 0, hi - lo, 1)
  sc <- function(M) as.numeric(sweep(sweep(M, 2, lo), 2, rng, `/`) %*% w /
                                 sum(w))
  list(score = sc(X), predict = sc, importance = w / sum(w),
       details = list(min = lo, max = hi))
}

.comp_pca <- function(X) {
  if (ncol(X) == 1) return(.comp_mean(X, c(1)))
  ctr <- colMeans(X); scl <- apply(X, 2, stats::sd)
  scl[scl == 0] <- 1
  pc <- stats::prcomp(X, center = ctr, scale. = scl)
  load <- pc$rotation[, 1]
  if (mean(load) < 0) load <- -load
  raw_sc <- function(M) as.numeric(sweep(sweep(M, 2, ctr), 2, scl, `/`) %*%
                                     load)
  s0 <- raw_sc(X); lo <- min(s0); hi <- max(s0)
  ve <- summary(pc)$importance[2, 1]
  list(score = .rescale01(s0), predict = function(M) .rescale01(raw_sc(M),
                                                                lo, hi),
       importance = load, details = list(var_explained = ve))
}

.comp_penalised <- function(X, y, component, cons, vars, alpha) {
  .need("glmnet", "method = \"penalised\"")
  if (is.null(y))
    stop("method = \"penalised\" needs an anchor: outcomes in the data or ",
         "an explicit `anchor`.", call. = FALSE)
  csign <- cons$comp_sign[[component]]
  if (is.na(csign)) csign <- if (component == "adaptive") -1 else 1
  y <- y * csign  # anchor aligned with the component
  sgn <- cons$ind_sign[vars]
  lower <- ifelse(!is.na(sgn) & sgn > 0, 0, -Inf)
  upper <- ifelse(!is.na(sgn) & sgn < 0, 0, Inf)
  for (v in intersect(names(cons$ind_bounds), vars)) {
    lower[v] <- max(lower[v], cons$ind_bounds[[v]][1])
    upper[v] <- min(upper[v], cons$ind_bounds[[v]][2])
  }
  ctr <- colMeans(X); scl <- apply(X, 2, stats::sd); scl[scl == 0] <- 1
  Z <- sweep(sweep(X, 2, ctr), 2, scl, `/`)
  nf <- min(10, nrow(Z))
  cv <- glmnet::cv.glmnet(Z, y, alpha = alpha, nfolds = nf,
                          lower.limits = lower, upper.limits = upper,
                          standardize = FALSE)
  beta <- as.numeric(stats::coef(cv, s = "lambda.min"))[-1]
  lam <- cv$lambda.min
  if (all(beta == 0)) {
    nz <- which(cv$glmnet.fit$df > 0)
    if (!length(nz))
      stop("Penalised fit selected no indicators for ", component, ".",
           call. = FALSE)
    lam <- cv$lambda[nz[1]]
    beta <- as.numeric(stats::coef(cv$glmnet.fit, s = lam))[-1]
    warning("lambda.min selected no ", component, " indicators; using the ",
            "largest lambda with a non-zero coefficient.", call. = FALSE)
  }
  names(beta) <- vars
  raw_sc <- function(M) as.numeric(sweep(sweep(M, 2, ctr), 2, scl, `/`) %*%
                                     beta)
  s0 <- raw_sc(X); lo <- min(s0); hi <- max(s0)
  list(score = .rescale01(s0),
       predict = function(M) .rescale01(raw_sc(M), lo, hi),
       importance = beta,
       details = list(lambda = lam, alpha = alpha,
                      n_selected = sum(beta != 0), cv = cv))
}

# Anchor: oriented vulnerability proxy (higher = more vulnerable)
.anchor <- function(x, anchor, required = TRUE) {
  if (is.numeric(anchor)) {
    if (length(anchor) != length(x$id))
      stop("`anchor` must have one value per unit.", call. = FALSE)
    return(anchor)
  }
  if (is.character(anchor)) {
    if (!anchor %in% names(x$data))
      stop("Anchor column '", anchor, "' not found.", call. = FALSE)
    d <- if (anchor %in% names(x$direction)) x$direction[[anchor]] else 1
    return(x$data[[anchor]] * d)
  }
  O <- .oriented_outcomes(x)
  if (is.null(O)) {
    if (required) return(NULL)
    return(NULL)
  }
  if (ncol(O) == 1) return(as.numeric(O[, 1]))
  pc <- stats::prcomp(O, scale. = TRUE)
  s <- pc$x[, 1]
  if (mean(pc$rotation[, 1]) < 0) s <- -s
  as.numeric(s)
}

.need <- function(pkg, what) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is needed for ", what,
         ". Install it with install.packages(\"", pkg, "\").", call. = FALSE)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
