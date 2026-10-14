#' Calibrated risk predictions
#'
#' Rankings say which species are *more* vulnerable than others; managers
#' need to know how *likely* a species is to decline. This calibrates a
#' fit's vulnerability scores against an observed outcome (for example IUCN
#' status of Vulnerable or worse, a declining population trend, or range
#' loss) by logistic regression, and returns for every unit the probability
#' of that outcome with a confidence interval.
#'
#' Intervals come from a bootstrap over units that refits the calibration
#' model. If a [ccva_uncertainty()] result is supplied, each replicate also
#' uses a different draw of the vulnerability scores, so uncertainty in the
#' scores themselves is propagated into the risk.
#'
#' When the outcome was also used to fit the vulnerability model (latent,
#' shape and penalised methods use the outcomes), this calibration is
#' in-sample and optimistic. For an honest assessment, calibrate against an
#' outcome the model did not see.
#'
#' @param fit A `ccva_fit`.
#' @param outcome Outcome column (default: the first outcome). Oriented so
#'   higher = more vulnerable.
#' @param threshold Units with an oriented outcome at or above this value
#'   count as an event. Default: for an ordinal outcome (10 or fewer
#'   distinct values) the second-highest value, e.g. IUCN status
#'   Vulnerable or worse on a 0--3 scale; otherwise the median.
#' @param event Optional logical vector of events, overriding `outcome`.
#' @param uncertainty Optional [ccva_uncertainty()] result.
#' @param R Number of bootstrap replicates.
#' @param level Interval coverage.
#' @param seed Optional random seed.
#' @return An object of class `ccva_calibration` with `units` (probability
#'   and interval per unit), `curve` (probability across the score range),
#'   `bins` (observed vs predicted), `metrics` (AUC, Brier score) and the
#'   `model`.
#' @export
#' @examples
#' fit <- ccva_combine(ccva_sim_data(), "weighted",
#'                     component_method = "ordinal")
#' cal <- ccva_calibrate(fit, outcome = "iucn", R = 200, seed = 1)
#' cal
ccva_calibrate <- function(fit, outcome = NULL, threshold = NULL,
                           event = NULL, uncertainty = NULL, R = 500,
                           level = 0.95, seed = NULL) {
  if (!inherits(fit, "ccva_fit")) stop("`fit` must be a ccva_fit.",
                                       call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  x <- fit$data
  n <- length(fit$score)
  if (is.null(event)) {
    if (!length(x$outcomes))
      stop("No outcomes in the data; supply `event`.", call. = FALSE)
    outcome <- outcome %||% x$outcomes[1]
    if (!outcome %in% names(x$data))
      stop("Outcome '", outcome, "' not found.", call. = FALSE)
    d <- if (outcome %in% names(x$direction)) x$direction[[outcome]] else 1
    y <- x$data[[outcome]] * d
    u <- sort(unique(y))
    if (is.null(threshold))
      threshold <- if (length(u) <= 10) u[length(u) - 1] else stats::median(y)
    event <- y >= threshold
    in_sample <- fit$method %in% c("latent", "shape", "penalised")
  } else {
    if (length(event) != n) stop("`event` needs one value per unit.",
                                 call. = FALSE)
    event <- as.logical(event); outcome <- "event"; in_sample <- FALSE
  }
  if (length(unique(event)) < 2)
    stop("The outcome defines only one class; change `threshold`.",
         call. = FALSE)

  s <- as.numeric(fit$score)
  glm_fit <- function(sc, ev) suppressWarnings(
    stats::glm(ev ~ sc, family = stats::binomial()))
  m <- glm_fit(s, event)
  prob <- as.numeric(stats::predict(m, type = "response"))
  grid <- seq(min(s), max(s), length.out = 60)

  draws <- if (!is.null(uncertainty)) uncertainty$draws else NULL
  P <- matrix(NA_real_, R, n); G <- matrix(NA_real_, R, length(grid))
  for (b in seq_len(R)) {
    idx <- sample.int(n, n, replace = TRUE)
    sb <- if (is.null(draws)) s else as.numeric(draws[(b - 1) %% nrow(draws) + 1, ])
    mb <- glm_fit(sb[idx], event[idx])
    cf <- stats::coef(mb)
    if (any(!is.finite(cf))) next
    P[b, ] <- stats::plogis(cf[1] + cf[2] * sb)
    G[b, ] <- stats::plogis(cf[1] + cf[2] * grid)
  }
  ok <- stats::complete.cases(P)
  P <- P[ok, , drop = FALSE]; G <- G[ok, , drop = FALSE]
  a <- (1 - level) / 2
  q <- function(M) apply(M, 2, stats::quantile, c(a, 1 - a), names = FALSE)
  qp <- q(P); qg <- q(G)

  units <- data.frame(unit = names(fit$score), score = s, probability = prob,
                      lower = qp[1, ], upper = qp[2, ], event = event,
                      row.names = NULL)
  units <- units[order(-units$probability), ]
  curve <- data.frame(score = grid,
                      probability = stats::plogis(stats::coef(m)[1] +
                                                    stats::coef(m)[2] * grid),
                      lower = qg[1, ], upper = qg[2, ])
  bins <- .calib_bins(prob, event)
  auc <- .auc(s, event)
  structure(list(units = units, curve = curve, bins = bins,
                 metrics = c(auc = auc, brier = mean((prob - event)^2),
                             event_rate = mean(event),
                             slope = unname(stats::coef(m)[2])),
                 model = m, outcome = outcome, threshold = threshold,
                 level = level, R = sum(ok), in_sample = in_sample,
                 propagated = !is.null(draws), method = fit$method),
            class = "ccva_calibration")
}

.calib_bins <- function(prob, event, k = 5) {
  br <- unique(stats::quantile(prob, seq(0, 1, length.out = k + 1)))
  if (length(br) < 3) br <- range(prob) + c(-1e-9, 1e-9)
  g <- cut(prob, br, include.lowest = TRUE)
  data.frame(predicted = tapply(prob, g, mean),
             observed = tapply(as.numeric(event), g, mean),
             n = as.numeric(table(g)), row.names = NULL)
}

.auc <- function(s, event) {
  r <- rank(s)
  n1 <- sum(event); n0 <- sum(!event)
  (sum(r[event]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' @export
print.ccva_calibration <- function(x, ...) {
  cat(sprintf("<ccva_calibration> P(%s %s %s) from %s vulnerability\n",
              x$outcome, ">=", format(x$threshold), x$method))
  m <- x$metrics
  cat(sprintf("  AUC %.3f | Brier %.3f | event rate %.2f | %d bootstrap fits%s\n",
              m[["auc"]], m[["brier"]], m[["event_rate"]], x$R,
              if (x$propagated) " (score uncertainty propagated)" else ""))
  if (x$in_sample)
    cat("  note: this outcome was used to fit the model; calibration is in-sample\n")
  top <- utils::head(x$units, 5)
  cat(sprintf("  highest risk (%d%% intervals):\n", round(100 * x$level)))
  for (i in seq_len(nrow(top)))
    cat(sprintf("    %-16s %.2f [%.2f, %.2f]\n", top$unit[i],
                top$probability[i], top$lower[i], top$upper[i]))
  invisible(x)
}

#' Match methods to the data available
#'
#' Turns general guidance -- shape-constrained parametric models
#' for small samples, regularised regression for regional compilations,
#' flexible non-parametric methods for global databases, with internal
#' consistency checks, structural or two-stage handling of dependencies
#' among components, and honest uncertainty in every case -- into concrete
#' advice for a given dataset, based on the number of units and the number
#' of predictors per component.
#'
#' The sample-size thresholds (50 and 150 units) reflect rough feasibility:
#' correlated additive latent models need about 50 units, interaction
#' surfaces about 150. They are guides, not rules.
#'
#' @param x A `ccva_data` object.
#' @return An object of class `ccva_recommendation` with the `tier`, a
#'   `profile` of the data and a table of `steps`.
#' @export
#' @examples
#' ccva_recommend(ccva_sim_data())
#' ccva_recommend(ccva_simulate(40, seed = 1))
ccva_recommend <- function(x) {
  .check_ccva(x)
  n <- length(x$id)
  comps <- .active_components(x)
  p <- vapply(x$indicators[comps], length, integer(1))
  hd <- names(p)[p / n >= 0.5]
  has_out <- length(x$outcomes) > 0
  tier <- if (n < 50) "small" else if (n < 150) "regional" else "global"

  step <- function(stage, recommendation, call)
    data.frame(stage = stage, recommendation = recommendation, call = call)
  steps <- list(step("diagnose",
    "Check internal consistency of each component (Cronbach's alpha, mean inter-item r); declare complementary-facet components as formative",
    "ccva_check(x)"))
  if (length(hd))
    steps <- c(steps, list(step("components",
      sprintf("Predictors outnumber or approach units for %s: use regularised (penalised) scoring", paste(hd, collapse = ", ")),
      sprintf("component_method = list(%s)",
              paste(sprintf("%s = \"penalised\"", hd), collapse = ", ")))))
  main <- switch(tier,
    small = step("vulnerability",
      "Small sample: shape-constrained parametric model with ecological sign constraints; weighted index as baseline",
      if (has_out) "ccva_combine(x, \"shape\")" else
        "ccva_combine(x, \"weighted\")"),
    regional = step("vulnerability",
      if (length(hd)) "Regional compilation with many predictors: regularised two-stage model"
      else "Regional compilation: constrained latent SEM (vulnerability inferred from outcomes)",
      if (!has_out) "ccva_combine(x, \"weighted\")" else if (length(hd))
        "ccva_combine(x, \"penalised\")" else "ccva_combine(x, \"latent\")"),
    global = step("vulnerability",
      "Large database: constrained latent SEM, plus flexible models -- interaction surfaces (adaptive capacity as moderator) or unconstrained smooths",
      if (has_out) "ccva_combine(x, \"latent\"); ccva_combine(x, \"shape\", constraints = ccva_constraints(constrain_moderator()))"
      else "ccva_combine(x, \"weighted\")"))
  steps <- c(steps, list(main))
  steps <- c(steps, list(step("dependencies",
    "Model dependencies among components you hypothesise (e.g. adaptive capacity altering exposure via range shifts, or sensitivity via plasticity)",
    "constraints = ccva_constraints(constrain_link(\"adaptive\", \"exposure\"))")))
  steps <- c(steps, list(step("uncertainty",
    "Propagate uncertainty to scores and ranks; report intervals, not a bare ranking",
    "ccva_uncertainty(fit, \"bootstrap\")")))
  if (has_out)
    steps <- c(steps, list(step("calibration",
      "Turn scores into calibrated risk (probability of an observed outcome) with intervals",
      "ccva_calibrate(fit, uncertainty = u)")))
  else
    steps <- c(steps, list(step("outcomes",
      "No outcome proxies: latent, shape and calibration steps are unavailable; add IUCN status, population trend or range change if possible",
      "as_ccva(..., outcomes = )")))
  if (n < 50)
    steps <- c(steps, list(step("caution",
      "With fewer than 50 units expect wide intervals and few robust rankings; constraints are what keep latent fits stable",
      "")))
  steps <- do.call(rbind, steps)
  structure(list(tier = tier,
                 profile = data.frame(component = comps, predictors = p,
                                      p_over_n = round(p / n, 2),
                                      row.names = NULL),
                 n = n, outcomes = x$outcomes, steps = steps),
            class = "ccva_recommendation")
}

#' @export
print.ccva_recommendation <- function(x, ...) {
  lab <- c(small = "small sample (< 50 units)",
           regional = "regional compilation (50-149 units)",
           global = "global database (>= 150 units)")[[x$tier]]
  cat("<ccva_recommendation>", lab, "\n")
  cat("  predictors:", paste(sprintf("%s %d (p/n %.2f)", x$profile$component,
                                     x$profile$predictors, x$profile$p_over_n),
                             collapse = ", "), "\n")
  cat("  outcomes:  ", if (length(x$outcomes))
    paste(x$outcomes, collapse = ", ") else "none", "\n\n")
  for (i in seq_len(nrow(x$steps))) {
    cat(sprintf("  %d. [%s] %s\n", i, x$steps$stage[i],
                x$steps$recommendation[i]))
    if (nzchar(x$steps$call[i])) cat("     ", x$steps$call[i], "\n")
  }
  invisible(x)
}

#' Adjust trait measurements for assay methodology
#'
#' Physiological sensitivity traits such as heat tolerance (CTmax) are
#' measured with different protocols -- static or dynamic assays, different
#' ramping rates, starting temperatures and acclimation regimes -- and the
#' protocol shifts the value obtained. Before species are compared, the
#' methodological effect should be removed. This fits a linear model of the
#' measurements on the method variables (with a unit effect when units are
#' measured more than once) and returns every measurement adjusted to a
#' common reference protocol, plus a per-unit summary ready to use as an
#' indicator in [as_ccva()].
#'
#' Method effects can only be separated from true differences between
#' units when protocols vary within units, or are not confounded with them
#' across units; aliased effects are reported.
#'
#' @param data A data frame of measurements (one row per measurement).
#' @param value Name of the measured trait column.
#' @param method_vars Names of the method columns (numeric, e.g. ramp rate;
#'   or categorical, e.g. static / dynamic).
#' @param unit Optional name of the unit column (species or population).
#'   Including it lets method effects be estimated within units.
#' @param reference Optional named list of reference values for the method
#'   variables. Defaults: the median of numeric variables, the most common
#'   level of categorical ones.
#' @param weights Optional name of a column of weights (e.g. sample size).
#' @return An object of class `ccva_adjustment` with `data` (including the
#'   adjusted column `<value>_adj`), `by_unit` (mean adjusted value per
#'   unit), `effects` (method coefficients), `reference` and `model`.
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(species = rep(paste0("sp", 1:20), each = 3),
#'                 ramp = rep(c(0.1, 0.25, 0.5), 20),
#'                 assay = rep(c("dynamic", "dynamic", "static"), 20))
#' d$ctmax <- 40 + rep(rnorm(20), each = 3) + 4 * d$ramp + rnorm(60, 0, 0.3)
#' adj <- ccva_adjust_method(d, "ctmax", c("ramp", "assay"), unit = "species")
#' adj
#' head(adj$by_unit)
ccva_adjust_method <- function(data, value, method_vars, unit = NULL,
                               reference = NULL, weights = NULL) {
  data <- as.data.frame(data)
  need <- c(value, method_vars, unit, weights)
  miss <- setdiff(need, names(data))
  if (length(miss)) stop("Columns not found: ", paste(miss, collapse = ", "),
                         call. = FALSE)
  if (!is.numeric(data[[value]])) stop("`value` must be numeric.",
                                       call. = FALSE)
  cc <- stats::complete.cases(data[need])
  if (!all(cc)) message("Dropping ", sum(!cc), " row(s) with missing values.")
  d <- data[cc, , drop = FALSE]
  for (v in method_vars) if (is.character(d[[v]]) || is.logical(d[[v]]))
    d[[v]] <- factor(d[[v]])
  if (!is.null(unit)) d[[unit]] <- factor(d[[unit]])

  ref <- lapply(method_vars, function(v) {
    if (!is.null(reference[[v]])) return(reference[[v]])
    if (is.factor(d[[v]])) names(which.max(table(d[[v]]))) else
      stats::median(d[[v]])
  })
  names(ref) <- method_vars

  rhs <- c(if (!is.null(unit) && anyDuplicated(d[[unit]])) unit, method_vars)
  fml <- stats::as.formula(paste(value, "~", paste(rhs, collapse = " + ")))
  w <- if (is.null(weights)) NULL else d[[weights]]
  m <- stats::lm(fml, data = d, weights = w)
  aliased <- names(stats::coef(m))[is.na(stats::coef(m))]

  nd <- d
  for (v in method_vars)
    nd[[v]] <- if (is.factor(d[[v]])) factor(ref[[v]], levels(d[[v]])) else
      ref[[v]]
  shift <- suppressWarnings(stats::predict(m, newdata = d) -
                              stats::predict(m, newdata = nd))
  shift[is.na(shift)] <- 0
  adj_name <- paste0(value, "_adj")
  d[[adj_name]] <- d[[value]] - shift

  cf <- suppressWarnings(stats::coef(summary(m)))
  eff_rows <- grep(paste0("^(", paste(method_vars, collapse = "|"), ")"),
                   rownames(cf))
  effects <- data.frame(term = rownames(cf)[eff_rows],
                        estimate = cf[eff_rows, 1], se = cf[eff_rows, 2],
                        pvalue = cf[eff_rows, 4], row.names = NULL)
  by_unit <- NULL
  if (!is.null(unit)) {
    by_unit <- data.frame(
      unit = levels(d[[unit]]),
      adjusted = as.numeric(tapply(d[[adj_name]], d[[unit]], mean)),
      raw = as.numeric(tapply(d[[value]], d[[unit]], mean)),
      n = as.numeric(table(d[[unit]])), row.names = NULL)
    names(by_unit)[2:3] <- c(adj_name, value)
  }
  if (length(aliased))
    warning("Method effect(s) not identifiable (confounded with units): ",
            paste(aliased, collapse = ", "), call. = FALSE)
  structure(list(data = d, by_unit = by_unit, effects = effects,
                 reference = ref, model = m, value = value,
                 aliased = aliased),
            class = "ccva_adjustment")
}

#' @export
print.ccva_adjustment <- function(x, ...) {
  cat(sprintf("<ccva_adjustment> %s adjusted to reference protocol (%s)\n",
              x$value, paste(sprintf("%s = %s", names(x$reference),
                                     unlist(x$reference)), collapse = ", ")))
  e <- x$effects
  for (i in seq_len(nrow(e)))
    cat(sprintf("  %-22s %+8.3f (se %.3f, p %.3g)\n", e$term[i],
                e$estimate[i], e$se[i], e$pvalue[i]))
  if (length(x$aliased))
    cat("  not identifiable:", paste(x$aliased, collapse = ", "), "\n")
  if (!is.null(x$by_unit))
    cat(sprintf("  %d units; use $by_unit as an indicator in as_ccva()\n",
                nrow(x$by_unit)))
  invisible(x)
}

#' Aggregate species vulnerability to sites or communities
#'
#' For spatial prioritisation: given which species occur where, summarise
#' the vulnerability of each site's community -- mean and maximum
#' vulnerability, the number and share of highly vulnerable species, and
#' the mean exposure, sensitivity and adaptive capacity of its species. With
#' a [ccva_uncertainty()] result, the site means get intervals too.
#'
#' @param fit A `ccva_fit`.
#' @param occurrence Either a data frame with one row per species presence
#'   at a site (columns named by `site` and `unit`), or a sites x species
#'   0/1 matrix with dimnames.
#' @param site,unit Column names in `occurrence`.
#' @param uncertainty Optional [ccva_uncertainty()] result.
#' @param coords Optional data frame with columns `site`, `x`, `y` for maps.
#' @param level Interval coverage.
#' @return An object of class `ccva_sites`: a data frame per site, ordered
#'   from most to least vulnerable community.
#' @export
#' @examples
#' fit <- ccva_combine(ccva_sim_data(), "weighted")
#' set.seed(1)
#' occ <- data.frame(site = rep(paste0("site", 1:8), each = 25),
#'                   species = sample(names(fit$score), 200, replace = TRUE))
#' occ <- unique(occ)
#' ccva_aggregate(fit, occ)
ccva_aggregate <- function(fit, occurrence, site = "site", unit = "species",
                           uncertainty = NULL, coords = NULL, level = 0.95) {
  if (!inherits(fit, "ccva_fit")) stop("`fit` must be a ccva_fit.",
                                       call. = FALSE)
  if (is.matrix(occurrence)) {
    if (is.null(dimnames(occurrence)[[1]]) || is.null(dimnames(occurrence)[[2]]))
      stop("The occurrence matrix needs site and species dimnames.",
           call. = FALSE)
    idx <- which(occurrence > 0, arr.ind = TRUE)
    occurrence <- data.frame(rownames(occurrence)[idx[, 1]],
                             colnames(occurrence)[idx[, 2]])
    names(occurrence) <- c(site, unit)
  }
  occurrence <- as.data.frame(occurrence)
  if (!all(c(site, unit) %in% names(occurrence)))
    stop("`occurrence` needs columns '", site, "' and '", unit, "'.",
         call. = FALSE)
  occ <- unique(occurrence[c(site, unit)])
  known <- occ[[unit]] %in% names(fit$score)
  if (!all(known))
    warning(sum(!known), " occurrence(s) of species not in the fit were ",
            "ignored.", call. = FALSE)
  occ <- occ[known, , drop = FALSE]
  if (!nrow(occ)) stop("No occurrences match the fitted species.",
                       call. = FALSE)
  sp <- as.character(occ[[unit]]); st <- as.character(occ[[site]])
  v <- fit$score[sp]
  cls <- stats::setNames(as.character(fit$class), names(fit$score))
  high <- cls[sp] == "High"
  cs <- fit$component_scores[sp, , drop = FALSE]
  sites <- unique(st)
  out <- data.frame(
    site = sites,
    n_species = as.numeric(tapply(v, st, length)[sites]),
    mean_vulnerability = as.numeric(tapply(v, st, mean)[sites]),
    max_vulnerability = as.numeric(tapply(v, st, max)[sites]),
    n_high = as.numeric(tapply(high, st, sum)[sites]),
    row.names = NULL)
  out$share_high <- out$n_high / out$n_species
  for (cmp in names(cs))
    out[[paste0("mean_", cmp)]] <- as.numeric(tapply(cs[[cmp]], st, mean)[sites])
  if (!is.null(uncertainty)) {
    D <- uncertainty$draws[, sp, drop = FALSE]
    a <- (1 - level) / 2
    M <- t(apply(D, 1, function(r) tapply(r, st, mean)[sites]))
    out$lower <- apply(M, 2, stats::quantile, a, names = FALSE)
    out$upper <- apply(M, 2, stats::quantile, 1 - a, names = FALSE)
    rk <- t(apply(-M, 1, rank, ties.method = "min"))
    out$rank_low <- apply(rk, 2, stats::quantile, a, names = FALSE)
    out$rank_high <- apply(rk, 2, stats::quantile, 1 - a, names = FALSE)
  }
  if (!is.null(coords)) {
    if (!all(c("site", "x", "y") %in% names(coords)))
      stop("`coords` needs columns site, x, y.", call. = FALSE)
    out <- merge(out, coords[c("site", "x", "y")], by = "site", all.x = TRUE)
  }
  out <- out[order(-out$mean_vulnerability), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "level") <- level
  class(out) <- c("ccva_sites", "data.frame")
  out
}

#' @export
print.ccva_sites <- function(x, ...) {
  cat(sprintf("<ccva_sites> %d sites, %d species occurrences\n", nrow(x),
              sum(x$n_species)))
  y <- as.data.frame(unclass(x))
  attr(y, "level") <- NULL
  keep <- intersect(c("site", "n_species", "mean_vulnerability", "lower",
                      "upper", "share_high", "rank_low", "rank_high"),
                    names(y))
  y <- y[keep]
  num <- vapply(y, is.numeric, logical(1))
  y[num] <- lapply(y[num], round, 3)
  print(utils::head(y, 10), row.names = FALSE)
  if (nrow(y) > 10) cat("  ...", nrow(y) - 10, "more sites\n")
  invisible(x)
}
