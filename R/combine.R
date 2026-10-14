#' Combine the components into vulnerability
#'
#' Scores the three components and combines them into a vulnerability score
#' for every unit, under a chosen model. Constraints from
#' [ccva_constraints()] are imposed on predictors and on each component's
#' role in every model that can express them.
#'
#' @section Methods:
#' \describe{
#'   \item{`"weighted"`}{The conventional index (Ofori et al. 2017) with the
#'     corrected normalisation \eqn{V = (w_E E + w_S S - w_A A + w_A) /
#'     (w_E + w_S + w_A)}, which equals \eqn{(E + S - A + 1)/3} for equal
#'     weights and maps the full range onto \[0, 1\] without clipping.
#'     Without adaptive-capacity indicators it is \eqn{(E + S)/2}. With
#'     [constrain_moderator()] adaptive capacity has no main effect:
#'     \eqn{V \propto w_E E (1 - \gamma_E A) + w_S S (1 - \gamma_S A)}.
#'     Use `component_method = "ordinal"` to reproduce Ofori's 1--3 scoring.}
#'   \item{`"shape"`}{Shape-constrained additive model (Pya and Wood 2015,
#'     package \pkg{scam}) of an observed `response` on the component
#'     scores. Component signs become monotone increasing / decreasing
#'     splines, so the *shape* is learned from data while the *direction* is
#'     fixed by ecology. With [constrain_moderator()] it fits the
#'     interaction surfaces \eqn{f_{EA}(E, A) + f_{SA}(S, A)}, monotone in
#'     each argument (adaptive capacity as moderator only).}
#'   \item{`"latent"`}{Structural equation model (\pkg{lavaan}) in which
#'     vulnerability \eqn{V^*} is a latent variable measured by the outcome
#'     proxies and caused by latent exposure, sensitivity and adaptive
#'     capacity, each measured by its own indicators (a full MIMIC model).
#'     Every constraint type is imposed as a linear (in)equality on
#'     loadings and paths. Indicators sharing an error source (e.g. SDM
#'     metrics with a common denominator) are declared in `shared_error`
#'     and modelled with an orthogonal method factor.}
#'   \item{`"penalised"`}{Two-stage model for high-dimensional predictors:
#'     (1) latent vulnerability is estimated from the outcome proxies alone;
#'     (2) it is regressed on *all* indicators with a LASSO / elastic net
#'     (\pkg{glmnet}) whose coefficient signs follow the component and
#'     indicator constraints. Selects which predictors drive vulnerability.}
#' }
#' Methods added with [ccva_register()] are also accepted.
#'
#' @param x A `ccva_data` object.
#' @param method Combination method; see Methods.
#' @param component_method Scoring method for the components: one string
#'   for all, or a named list / vector per component (see
#'   [ccva_component()]).
#' @param component_args Optional named list (per component) of extra
#'   arguments for [ccva_component()], e.g. `weights` or `breaks`.
#' @param weights Component weights for `"weighted"`, named `exposure`,
#'   `sensitivity`, `adaptive`. Default equal.
#' @param constraints A [ccva_constraints()] object. Defaults to
#'   [ccva_default_constraints()]; use [ccva_unconstrained()] for free
#'   estimation.
#' @param response For `"shape"`: an outcome column, a numeric vector
#'   (higher = more vulnerable), or `NULL` to use the first principal
#'   component of all outcomes.
#' @param shared_error For `"latent"`: a list of character vectors, each a
#'   group of indicators with a common error source.
#' @param max_indicators For `"latent"`: a component with more indicators
#'   than this (or with only one) enters as its observed component score
#'   instead of a latent factor. Default 30.
#' @param k Basis dimension for `"shape"` splines; chosen from `n` if `NULL`.
#' @param alpha Elastic-net mixing for `"penalised"` (1 = LASSO).
#' @param cutoffs Classification of scores into Low / Moderate / High; see
#'   [ccva_classify()].
#' @param ... Passed to registered combine methods.
#'
#' @return An object of class `ccva_fit`: `score` (vulnerability in
#'   \[0, 1\], higher = more vulnerable), `class`, `components` (fitted
#'   [ccva_component()] objects), `component_scores`, `model`, `details`,
#'   `notes` (any constraint a method could not express) and `predict`.
#' @references
#' Ofori, B. Y. et al. (2017) *Scientific Reports* 7, 12979.
#' Pya, N. and Wood, S. N. (2015) *Statistics and Computing* 25, 543--559.
#' Rosseel, Y. (2012) *Journal of Statistical Software* 48(2).
#' @seealso [ccva()], [ccva_compare()], [ccva_uncertainty()], [ccva_plot()]
#' @export
#' @examples
#' d <- ccva_sim_data()
#'
#' # Conventional Ofori-style weighted index
#' w <- ccva_combine(d, "weighted", component_method = "ordinal")
#' w
#'
#' # Constrained latent model (the default constraints)
#' l <- ccva_combine(d, "latent",
#'                   shared_error = list(c("E1_risk", "E2_risk", "E3_risk")))
#' l
#' cor(l$score, d$truth$V_star, method = "spearman")
ccva_combine <- function(x, method = c("weighted", "shape", "latent",
                                       "penalised"),
                         component_method = "mean", component_args = list(),
                         weights = NULL, constraints = NULL,
                         response = NULL, shared_error = NULL,
                         max_indicators = 30, k = NULL, alpha = 1,
                         cutoffs = "uniform", ...) {
  .check_ccva(x)
  args <- list(method = if (is.character(method)) method[1] else method,
               component_method = component_method,
               component_args = component_args, weights = weights,
               constraints = constraints, response = response,
               shared_error = shared_error, max_indicators = max_indicators,
               k = k, alpha = alpha, cutoffs = cutoffs, ...)
  method <- args$method
  if (is.null(constraints)) constraints <- ccva_default_constraints()
  cons <- .resolve_constraints(constraints, x)
  comps <- .active_components(x)

  ## 1. component scores ------------------------------------------------------
  cm <- .per_component(component_method, comps)
  components <- lapply(comps, function(cmp) {
    a <- c(list(x = x, component = cmp, method = cm[[cmp]],
                constraints = constraints, alpha = alpha),
           component_args[[cmp]])
    do.call(ccva_component, a)
  })
  names(components) <- comps
  cs <- as.data.frame(lapply(components, `[[`, "score"))
  rownames(cs) <- x$id

  ## 2. links between components (two-stage outside the SEM) -----------------
  notes <- character(0)
  if (length(cons$comp_shape) && method != "shape")
    notes <- c(notes, "constrain_shape() applies only to method = \"shape\" and was not used")
  link_info <- NULL
  if (length(cons$links) && method != "latent") {
    if (method == "penalised") {
      notes <- c(notes, "component links are not used by the penalised model (it works on indicators)")
    } else {
      lk <- .two_stage_links(components, cons$links)
      components <- lk$components
      link_info <- lk$info
      cs <- as.data.frame(lapply(components, `[[`, "score"))
      rownames(cs) <- x$id
      notes <- c(notes, sprintf("two-stage link: %s adjusted for %s (R2 = %.2f)",
                                lk$info$to, lk$info$from, lk$info$r2))
    }
  }

  ## 3. combination -----------------------------------------------------------
  res <- switch(method,
    weighted = .comb_weighted(cs, comps, weights, cons, components),
    shape = .comb_shape(x, cs, comps, cons, response, k, components),
    latent = .comb_latent(x, components, comps, cons, shared_error,
                          max_indicators),
    penalised = .comb_penalised(x, comps, cons, alpha),
    {
      reg <- .ccva_registry$combine[[method]]
      if (is.null(reg))
        stop("Unknown combine method '", method, "'. See ccva_methods().",
             call. = FALSE)
      out <- reg$fun(x, components = cs, constraints = cons, ...)
      if (is.null(out$score))
        stop("Registered combine method must return `score`.", call. = FALSE)
      out
    })
  notes <- c(notes, res$notes)

  score <- stats::setNames(.rescale01(as.numeric(res$score)), x$id)
  lo <- min(res$score); hi <- max(res$score)
  if (method == "weighted") {
    # the weighted index is already on its natural [0, 1] scale
    score <- stats::setNames(as.numeric(res$score), x$id)
    lo <- 0; hi <- 1
  }
  cls <- ccva_classify(score, cutoffs, data = x)

  fit <- structure(list(
    data = x, method = method, args = args,
    score = score, class = cls,
    components = components, component_scores = cs,
    model = res$model, details = res$details,
    constraints = cons, notes = notes,
    raw_range = c(lo, hi)
  ), class = "ccva_fit")
  if (!is.null(link_info)) fit$details$links <- link_info
  fit$predict_raw <- res$predict
  fit
}

# Two-stage resolution of directed links: the target component is replaced
# by its part not explained by the source (rescaled to [0, 1]).
.two_stage_links <- function(components, links) {
  info <- list()
  for (l in links) {
    src <- components[[l$from]]; tgt <- components[[l$to]]
    fit <- stats::lm(tgt$score ~ src$score)
    b <- stats::coef(fit)[2]
    ms <- mean(src$score)
    adj <- tgt$score - b * (src$score - ms)
    lo <- min(adj); hi <- max(adj)
    p_t <- tgt$predict; p_s <- src$predict
    tgt$score <- stats::setNames(.rescale01(adj), names(tgt$score))
    tgt$predict <- function(newdata)
      .rescale01(p_t(newdata) - b * (p_s(newdata) - ms), lo, hi)
    tgt$adjusted_for <- l$from
    components[[l$to]] <- tgt
    info[[length(info) + 1]] <- data.frame(
      from = l$from, to = l$to, slope = unname(b),
      r2 = summary(fit)$r.squared, method = "two-stage")
  }
  list(components = components, info = do.call(rbind, info))
}

#' Run a complete CCVA
#'
#' Convenience wrapper: builds the `ccva_data` (if given a data frame) and
#' calls [ccva_combine()].
#'
#' @param data A `ccva_data` object, or a data frame together with the
#'   component arguments of [as_ccva()].
#' @param method Combination method, see [ccva_combine()].
#' @param ... Arguments for [as_ccva()] (when `data` is a data frame) and
#'   [ccva_combine()].
#' @return A `ccva_fit` object.
#' @export
#' @examples
#' fit <- ccva(ccva_sim_data(), method = "weighted",
#'             component_method = "ordinal")
#' head(ccva_rank(fit))
ccva <- function(data, method = "latent", ...) {
  dots <- list(...)
  if (!inherits(data, "ccva_data")) {
    spec <- intersect(names(dots), names(formals(as_ccva)))
    data <- do.call(as_ccva, c(list(data = data), dots[spec]))
    dots <- dots[setdiff(names(dots), spec)]
  }
  do.call(ccva_combine, c(list(x = data, method = method), dots))
}

#' @export
print.ccva_fit <- function(x, ...) {
  lab <- c(weighted = "weighted index", shape = "shape-constrained",
           latent = "latent SEM", penalised = "penalised two-stage")
  m <- if (x$method %in% names(lab)) lab[[x$method]] else x$method
  cat(sprintf("<ccva_fit> %s | %d units | constraints: %d rule(s)\n", m,
              length(x$score), x$constraints$n_rules))
  cm <- vapply(x$components, `[[`, character(1), "method")
  cat("  components:", paste(sprintf("%s (%s)", names(cm), cm),
                             collapse = ", "), "\n")
  tb <- table(x$class)
  cat("  classes:   ", paste(sprintf("%s %d", names(tb), tb),
                             collapse = " | "), "\n")
  if (!is.null(x$details$fit_measures)) {
    fm <- x$details$fit_measures
    cat(sprintf("  fit: CFI %.3f, TLI %.3f, RMSEA %.3f, SRMR %.3f%s\n",
                fm["cfi"], fm["tli"], fm["rmsea"], fm["srmr"],
                if (isTRUE(x$details$converged)) "" else "  (NOT converged)"))
  }
  if (!is.null(x$details$effects)) {
    e <- x$details$effects
    cat("  effects on V*:", paste(sprintf("%s %+.2f", e$component,
                                          e$std), collapse = ", "),
        "(standardised)\n")
  }
  if (!is.null(x$details$links)) {
    lk <- x$details$links
    for (i in seq_len(nrow(lk))) {
      if (lk$method[i] == "structural path")
        cat(sprintf("  link %s -> %s: %+.2f (standardised)\n", lk$from[i],
                    lk$to[i], lk$std[i]))
      else
        cat(sprintf("  link %s -> %s: two-stage, %.0f%% of %s explained by %s removed\n",
                    lk$from[i], lk$to[i], 100 * lk$r2[i], lk$to[i], lk$from[i]))
    }
  }
  if (!is.null(x$details$total_effects)) {
    te <- x$details$total_effects
    for (i in seq_len(nrow(te)))
      cat(sprintf("  total effect of %s on V*: %+.2f = direct %+.2f + indirect %+.2f (%s)\n",
                  te$component[i], te$total[i], te$direct[i], te$indirect[i],
                  te$routes[i]))
  }
  if (!is.null(x$details$n_selected))
    cat(sprintf("  selected predictors: %d of %d\n", x$details$n_selected,
                nrow(x$details$coefficients)))
  if (isTRUE(x$details$sign_flipped))
    cat("  note: latent scores were sign-aligned to the outcomes\n")
  if (length(x$notes)) cat(paste0("  note: ", x$notes, "\n"), sep = "")
  invisible(x)
}

#' Predict vulnerability for new units
#'
#' @param object A `ccva_fit`.
#' @param newdata A `ccva_data` object (or data frame with the same columns)
#'   for the new units.
#' @param ... Unused.
#' @return Numeric vector of vulnerability on the scale of the training fit
#'   (values outside \[0, 1\] indicate units beyond the training range).
#' @export
predict.ccva_fit <- function(object, newdata, ...) {
  nd <- if (inherits(newdata, "ccva_data")) newdata$data else newdata
  if (is.null(object$predict_raw))
    stop("This method does not support prediction.", call. = FALSE)
  raw <- object$predict_raw(nd)
  .rescale01(raw, object$raw_range[1], object$raw_range[2])
}

#' Rank units by vulnerability
#'
#' @param fit A `ccva_fit`.
#' @param uncertainty Optional [ccva_uncertainty()] result to add intervals.
#' @return Data frame ordered from most to least vulnerable.
#' @export
#' @examples
#' fit <- ccva_combine(ccva_sim_data(), "weighted")
#' head(ccva_rank(fit))
ccva_rank <- function(fit, uncertainty = NULL) {
  out <- data.frame(unit = names(fit$score), score = fit$score,
                    rank = rank(-fit$score, ties.method = "min"),
                    class = fit$class, fit$component_scores,
                    row.names = NULL)
  if (!is.null(uncertainty)) {
    u <- uncertainty$summary
    out <- cbind(out, u[match(out$unit, u$unit),
                        c("lower", "upper", "rank_low", "rank_high")])
  }
  out[order(-out$score), , drop = FALSE]
}

## Combination implementations ---------------------------------------------------

.per_component <- function(spec, comps) {
  if (is.character(spec) && length(spec) == 1)
    return(stats::setNames(as.list(rep(spec, length(comps))), comps))
  spec <- as.list(spec)
  out <- stats::setNames(as.list(rep("mean", length(comps))), comps)
  out[intersect(names(spec), comps)] <- spec[intersect(names(spec), comps)]
  out
}

.comb_weighted <- function(cs, comps, weights, cons, components) {
  w <- stats::setNames(rep(1, length(comps)), comps)
  if (!is.null(weights)) {
    unk <- setdiff(names(weights), .ccva_components)
    if (length(unk)) stop("Unknown component weight(s): ",
                          paste(unk, collapse = ", "), call. = FALSE)
    w[intersect(names(weights), comps)] <-
      weights[intersect(names(weights), comps)]
  }
  for (cmp in names(cons$comp_fixed))
    if (cmp %in% comps) w[cmp] <- cons$comp_fixed[[cmp]]
  if (any(w < 0)) stop("Component weights must be non-negative.",
                       call. = FALSE)
  s <- cons$comp_sign
  s[is.na(s)] <- ifelse(names(s)[is.na(s)] == "adaptive", -1, 1)
  notes <- character(0)
  if (length(cons$comp_bounds) || length(cons$equal) ||
      length(cons$ind_bounds))
    notes <- "bounds / equality rules are not expressible in the weighted index and were not applied"

  index <- function(C) {
    C <- as.matrix(C[comps])
    if (!is.null(cons$moderator) && "adaptive" %in% comps) {
      g <- cons$moderator
      A <- C[, "adaptive"]
      num <- 0; den <- 0
      if ("exposure" %in% comps) {
        num <- num + w[["exposure"]] * C[, "exposure"] * (1 - g[1] * A)
        den <- den + w[["exposure"]]
      }
      num <- num + w[["sensitivity"]] * C[, "sensitivity"] * (1 - g[2] * A)
      den <- den + w[["sensitivity"]]
      return(num / den)
    }
    (C %*% (w * s) + sum(w[s < 0])) / sum(w)
  }
  list(score = as.numeric(index(cs)), model = list(weights = w, signs = s),
       details = list(weights = w, signs = s,
                      moderator = cons$moderator),
       notes = notes,
       predict = function(newdata)
         as.numeric(index(.component_newdata(components, newdata))))
}

.comb_shape <- function(x, cs, comps, cons, response, k, components) {
  .need("scam", "method = \"shape\"")
  y <- .anchor(x, response)
  if (is.null(y)) stop("method = \"shape\" needs outcomes or a `response`.",
                       call. = FALSE)
  n <- nrow(cs)
  df <- cs; df$.y <- y
  if (is.null(k)) k <- max(4, min(10, floor(n / (3 * 4))))
  notes <- character(0)
  if (length(cons$comp_bounds) || length(cons$equal))
    notes <- "bounds / equality rules are not expressible in the shape model and were not applied"

  terms <- character(0); used <- comps
  if (!is.null(cons$moderator)) {
    if (!all(c("sensitivity", "adaptive") %in% comps))
      stop("Moderator structure needs sensitivity and adaptive capacity.",
           call. = FALSE)
    kk <- if (n >= 150) 5 else 4
    df$.Ainv <- -df$adaptive  # decreasing in A == increasing in -A
    if ("exposure" %in% comps)
      terms <- c(terms, sprintf("s(exposure, .Ainv, bs = \"tedmi\", k = c(%d, %d))", kk, kk))
    terms <- c(terms, sprintf("s(sensitivity, .Ainv, bs = \"tedmi\", k = c(%d, %d))", kk, kk))
    if (n < 150)
      notes <- c(notes, "interaction surfaces need many units (the report suggests n >= 150); interpret with care")
  } else {
    for (cmp in comps) {
      fx <- cons$comp_fixed[[cmp]]
      if (!is.null(fx) && fx == 0) { used <- setdiff(used, cmp); next }
      sg <- cons$comp_sign[[cmp]]
      cv <- cons$comp_shape[[cmp]]
      bs <- if (is.null(cv)) {
        if (is.na(sg)) "cr" else if (sg > 0) "mpi" else "mpd"
      } else {
        paste0(if (is.na(sg)) "" else if (sg > 0) "mi" else "md",
               if (cv == "concave") "cv" else "cx")
      }
      terms <- c(terms, sprintf("s(%s, bs = \"%s\", k = %d)", cmp, bs, k))
    }
  }
  fml <- stats::as.formula(paste(".y ~", paste(terms, collapse = " + ")))
  fit <- scam::scam(fml, data = df)
  sc <- as.numeric(stats::fitted(fit))

  # partial effect curves for plotting (others at their mean)
  curves <- do.call(rbind, lapply(used, function(cmp) {
    g <- as.data.frame(lapply(cs[comps], function(v) rep(mean(v), 50)))
    g[[cmp]] <- seq(min(cs[[cmp]]), max(cs[[cmp]]), length.out = 50)
    if (!is.null(cons$moderator)) g$.Ainv <- -g$adaptive
    data.frame(component = cmp, x = g[[cmp]],
               y = as.numeric(stats::predict(fit, newdata = g)))
  }))

  comp_list <- components
  list(score = sc, model = fit,
       details = list(formula = fml, curves = curves,
                      r2 = summary(fit)$r.sq,
                      moderator = !is.null(cons$moderator)),
       notes = notes,
       predict = function(newdata) {
         nd <- .component_newdata(comp_list, newdata)
         if (!is.null(cons$moderator)) nd$.Ainv <- -nd$adaptive
         as.numeric(stats::predict(fit, newdata = nd))
       })
}

# Component scores of new units from the fitted component scorers
.component_newdata <- function(components, newdata) {
  if (is.null(components)) stop("Prediction for this fit needs components.",
                                call. = FALSE)
  as.data.frame(lapply(components, function(cp) cp$predict(newdata)))
}

.lv_name <- function(v) paste0("l_", gsub("[^A-Za-z0-9]", "_", v))

.comb_latent <- function(x, components, comps, cons, shared_error,
                         max_indicators) {
  if (!is.null(cons$moderator))
    stop("The latent SEM cannot express adaptive capacity as a moderator. ",
         "Use method = \"shape\" (interaction surfaces) or \"weighted\".",
         call. = FALSE)
  outs <- x$outcomes
  if (length(outs) < 1)
    stop("method = \"latent\" needs at least one outcome proxy.",
         call. = FALSE)
  n <- length(x$id)

  # oriented, standardised data
  cols <- c(unlist(x$indicators[comps]), outs)
  D <- sweep(as.matrix(x$data[cols]), 2, x$direction[cols], `*`)
  mu <- colMeans(D); sdv <- apply(D, 2, stats::sd); sdv[sdv == 0] <- 1
  Z <- as.data.frame(sweep(sweep(D, 2, mu), 2, sdv, `/`))

  lines <- character(0); ineq <- character(0); notes <- character(0)
  observed_comp <- character(0)
  loading_term <- function(v) {
    fx <- cons$ind_fixed[[v]]
    if (!is.null(fx)) return(if (fx == 0) NULL else sprintf("%s*%s", fx, v))
    lab <- .lv_name(v)
    for (g in cons$equal) if (v %in% g) lab <- .lv_name(paste(g, collapse = "_eq_"))
    sg <- cons$ind_sign[[v]]
    if (!is.na(sg)) ineq <<- c(ineq, sprintf("%s %s 0", lab,
                                             if (sg > 0) ">" else "<"))
    b <- cons$ind_bounds[[v]]
    if (!is.null(b)) {
      if (is.finite(b[1])) ineq <<- c(ineq, sprintf("%s > %s", lab, b[1]))
      if (is.finite(b[2])) ineq <<- c(ineq, sprintf("%s < %s", lab, b[2]))
    }
    sprintf("%s*%s", lab, v)
  }

  # measurement models for components
  struct_rhs <- character(0)
  for (cmp in comps) {
    vars <- components[[cmp]]$indicators
    lv <- .ccva_latent[[cmp]]
    if (cmp %in% x$formative || length(vars) < 2 ||
        length(vars) > max_indicators) {
      Z[[lv]] <- as.numeric(scale(components[[cmp]]$score))
      observed_comp <- c(observed_comp, cmp)
      notes <- c(notes, if (cmp %in% x$formative)
        sprintf("%s is formative: entered as a composite of its %d indicators",
                cmp, length(vars))
        else sprintf("%s entered as its observed component score (%d indicator%s)",
                     cmp, length(vars), if (length(vars) == 1) "" else "s"))
    } else {
      tm <- unlist(lapply(vars, loading_term))
      lines <- c(lines, sprintf("%s =~ %s", lv, paste(tm, collapse = " + ")))
    }
    struct_rhs <- c(struct_rhs, lv)
  }
  # vulnerability
  if (length(outs) >= 2) {
    tm <- unlist(lapply(outs, loading_term))
    lines <- c(lines, sprintf("Vuln =~ %s", paste(tm, collapse = " + ")))
    vname <- "Vuln"
  } else {
    vname <- outs
    notes <- c(notes, "only one outcome: vulnerability is that observed outcome")
  }
  # shared error (method factors)
  latent_comps <- .ccva_latent[setdiff(comps, observed_comp)]
  if (!is.null(shared_error)) {
    for (i in seq_along(shared_error)) {
      g <- intersect(shared_error[[i]], names(Z))
      if (length(g) < 2) next
      mf <- paste0("Shared", i)
      lines <- c(lines,
                 sprintf("%s =~ %s", mf, paste0("1*", g, collapse = " + ")),
                 sprintf("%s ~~ start(0.1)*%s", mf, mf))
      if (length(latent_comps))
        lines <- c(lines, sprintf("%s ~~ %s", mf,
                                  paste0("0*", latent_comps,
                                         collapse = " + ")))
    }
  }
  # structural paths
  sterms <- character(0)
  for (j in seq_along(comps)) {
    cmp <- comps[j]; lv <- struct_rhs[j]
    fx <- cons$comp_fixed[[cmp]]
    if (!is.null(fx)) {
      if (fx != 0) sterms <- c(sterms, sprintf("%s*%s", fx, lv))
      next
    }
    lab <- paste0("b_", cmp)
    sg <- cons$comp_sign[[cmp]]
    if (!is.na(sg)) ineq <- c(ineq, sprintf("%s %s 0", lab,
                                            if (sg > 0) ">" else "<"))
    b <- cons$comp_bounds[[cmp]]
    if (!is.null(b)) {
      if (is.finite(b[1])) ineq <- c(ineq, sprintf("%s > %s", lab, b[1]))
      if (is.finite(b[2])) ineq <- c(ineq, sprintf("%s < %s", lab, b[2]))
    }
    sterms <- c(sterms, sprintf("%s*%s", lab, lv))
  }
  lines <- c(lines, sprintf("%s ~ %s", vname, paste(sterms, collapse = " + ")))
  # directed links between components, with indirect and total effects
  lv_of <- stats::setNames(struct_rhs, comps)
  linked_pairs <- character(0)
  defs <- character(0)
  for (l in cons$links) {
    klab <- sprintf("k_%s_%s", l$from, l$to)
    lines <- c(lines, sprintf("%s ~ %s*%s", lv_of[[l$to]], klab,
                              lv_of[[l$from]]))
    if (l$sign != "free")
      ineq <- c(ineq, sprintf("%s %s 0", klab, if (l$sign == "+") ">" else "<"))
    linked_pairs <- c(linked_pairs, paste(sort(c(lv_of[[l$from]],
                                                 lv_of[[l$to]])),
                                          collapse = "|"))
    if (is.null(cons$comp_fixed[[l$to]]))
      defs <- c(defs, sprintf("ind_%s_%s := %s*b_%s", l$from, l$to, klab,
                              l$to))
  }
  # correlated latent components (not for pairs joined by a directed link)
  if (length(latent_comps) >= 2) {
    pr <- utils::combn(latent_comps, 2)
    keep <- !apply(pr, 2, function(v) paste(sort(v), collapse = "|")) %in%
      linked_pairs
    if (any(keep))
      lines <- c(lines, sprintf("%s ~~ %s", pr[1, keep], pr[2, keep]))
  }
  ineq <- c(ineq, defs)
  model <- paste(c(lines, unique(ineq)), collapse = "\n")

  fit <- .lavaan_fit(model, Z)
  conv <- lavaan::lavInspect(fit, "converged")
  if (!conv) notes <- c(notes, "lavaan did not converge; treat results with caution")

  anchor <- .anchor(x, NULL)
  lat <- tryCatch(lavaan::lavPredict(fit, se = "standard"),
                  error = function(e) lavaan::lavPredict(fit))
  score_se <- NULL
  if (vname == "Vuln" && !is.null(attr(lat, "se"))) {
    se <- attr(lat, "se")[[1]]
    score_se <- as.numeric(se[1, "Vuln"])
  }
  vs <- if (vname == "Vuln") as.numeric(lat[, "Vuln"]) else Z[[vname]]
  flip <- FALSE
  if (stats::cor(vs, anchor) < 0) { vs <- -vs; flip <- TRUE }

  pe <- lavaan::parameterEstimates(fit, standardized = TRUE)
  st <- pe[pe$op == "~" & pe$lhs == vname, ]
  effects <- data.frame(component = comps[match(st$rhs, struct_rhs)],
                        estimate = st$est, se = st$se,
                        std = st$std.all, pvalue = st$pvalue,
                        row.names = NULL)
  link_tab <- NULL; total_tab <- NULL
  if (length(cons$links)) {
    link_tab <- do.call(rbind, lapply(cons$links, function(l) {
      k <- pe[pe$op == "~" & pe$lhs == lv_of[[l$to]] &
                pe$rhs == lv_of[[l$from]], ]
      ind <- pe[pe$op == ":=" & pe$lhs == sprintf("ind_%s_%s", l$from, l$to), ]
      data.frame(from = l$from, to = l$to, estimate = k$est, std = k$std.all,
                 pvalue = k$pvalue,
                 indirect = if (nrow(ind)) ind$std.all else NA_real_,
                 method = "structural path", row.names = NULL)
    }))
    # total effect of each linked source component on V*:
    # direct path + one indirect path through every component it modifies
    srcs <- unique(link_tab$from)
    total_tab <- do.call(rbind, lapply(srcs, function(cmp) {
      d <- effects$std[effects$component == cmp]
      d <- if (length(d)) d else 0
      lk <- link_tab[link_tab$from == cmp, ]
      via <- stats::setNames(lk$indirect, paste0("via_", lk$to))
      data.frame(component = cmp, direct = d,
                 indirect = sum(lk$indirect, na.rm = TRUE),
                 total = d + sum(lk$indirect, na.rm = TRUE),
                 routes = paste(sprintf("via %s %+.2f", lk$to, lk$indirect),
                                collapse = "; "),
                 row.names = NULL)
    }))
  }
  ld <- pe[pe$op == "=~" & !grepl("^Shared", pe$lhs), ]
  loadings <- data.frame(factor = ld$lhs, indicator = ld$rhs,
                         estimate = ld$est, std = ld$std.all,
                         row.names = NULL)
  fm <- tryCatch(lavaan::fitMeasures(fit, c("cfi", "tli", "rmsea", "srmr")),
                 error = function(e) c(cfi = NA, tli = NA, rmsea = NA,
                                       srmr = NA))
  latent_scores <- NULL
  if (length(latent_comps)) {
    latent_scores <- as.data.frame(lat[, latent_comps, drop = FALSE])
    names(latent_scores) <- names(latent_comps)
  }

  list(score = vs, model = fit,
       details = list(syntax = model, converged = conv,
                      fit_measures = unlist(fm), effects = effects,
                      loadings = loadings, sign_flipped = flip,
                      latent_scores = latent_scores,
                      score_se = score_se, links = link_tab,
                      total_effects = total_tab,
                      n_inequality = length(unique(ineq)),
                      observed_components = observed_comp),
       notes = notes,
       predict = function(newdata) {
         Dn <- sweep(as.matrix(newdata[cols]), 2, x$direction[cols], `*`)
         Zn <- as.data.frame(sweep(sweep(Dn, 2, mu), 2, sdv, `/`))
         for (cmp in observed_comp) {
           s_new <- components[[cmp]]$predict(newdata)
           s_old <- components[[cmp]]$score
           Zn[[.ccva_latent[[cmp]]]] <- (s_new - mean(s_old)) / stats::sd(s_old)
         }
         p <- lavaan::lavPredict(fit, newdata = Zn)
         v <- if (vname == "Vuln") as.numeric(p[, "Vuln"]) else Zn[[vname]]
         if (flip) -v else v
       })
}

.lavaan_fit <- function(model, data) {
  warn <- character(0)
  fit <- withCallingHandlers(
    lavaan::sem(model, data = data, std.lv = TRUE, estimator = "ML",
                control = list(iter.max = 5000)),
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  attr(fit, "ccva_warnings") <- warn
  fit
}

.comb_penalised <- function(x, comps, cons, alpha) {
  .need("glmnet", "method = \"penalised\"")
  outs <- x$outcomes
  if (!length(outs)) stop("method = \"penalised\" needs outcome proxies.",
                          call. = FALSE)
  # Stage 1: latent vulnerability from outcomes alone
  O <- .oriented_outcomes(x)
  Oz <- as.data.frame(scale(O))
  stage1 <- NULL
  if (ncol(O) >= 3) {
    tm <- paste0("l_", seq_along(outs), "*", outs)
    model <- paste0("Vuln =~ ", paste(tm, collapse = " + "), "\n",
                    paste0("l_", seq_along(outs), " > 0", collapse = "\n"))
    stage1 <- tryCatch(.lavaan_fit(model, Oz), error = function(e) NULL)
  }
  if (!is.null(stage1) && lavaan::lavInspect(stage1, "converged")) {
    v1 <- as.numeric(lavaan::lavPredict(stage1)[, "Vuln"])
    s1 <- "latent factor of outcomes"
  } else {
    v1 <- .anchor(x, NULL); s1 <- "first principal component of outcomes"
  }
  if (stats::cor(v1, .anchor(x, NULL)) < 0) v1 <- -v1

  # Stage 2: penalised regression on all indicators
  vars <- unlist(x$indicators[comps])
  drop <- names(Filter(function(v) v == 0, cons$ind_fixed))
  vars <- setdiff(vars, drop)
  X <- do.call(cbind, lapply(comps, function(cmp) .oriented(x, cmp)))
  X <- X[, vars, drop = FALSE]
  comp_of <- vapply(vars, .component_of, character(1), x = x)
  csign <- cons$comp_sign[comp_of]
  csign[is.na(csign)] <- ifelse(comp_of[is.na(csign)] == "adaptive", -1, 1)
  isign <- cons$ind_sign[vars]
  eff <- ifelse(is.na(isign), NA, csign * isign)
  lower <- ifelse(!is.na(eff) & eff > 0, 0, -Inf)
  upper <- ifelse(!is.na(eff) & eff < 0, 0, Inf)
  ctr <- colMeans(X); scl <- apply(X, 2, stats::sd); scl[scl == 0] <- 1
  Zx <- sweep(sweep(X, 2, ctr), 2, scl, `/`)
  cv <- glmnet::cv.glmnet(Zx, v1, alpha = alpha, nfolds = min(10, nrow(Zx)),
                          lower.limits = lower, upper.limits = upper,
                          standardize = FALSE)
  beta <- as.numeric(stats::coef(cv, s = "lambda.min"))
  names(beta) <- c("(Intercept)", vars)
  b <- beta[-1]
  sc <- as.numeric(Zx %*% b)
  coefs <- data.frame(indicator = vars, component = comp_of, coef = b,
                      row.names = NULL)
  list(score = sc, model = cv,
       details = list(stage1 = s1, coefficients = coefs,
                      n_selected = sum(b != 0), lambda = cv$lambda.min,
                      alpha = alpha,
                      cv_r2 = 1 - min(cv$cvm) / stats::var(v1)),
       notes = character(0),
       predict = function(newdata) {
         Xn <- do.call(cbind, lapply(comps, function(cmp)
           .oriented(x, cmp, data = newdata)))[, vars, drop = FALSE]
         as.numeric(sweep(sweep(Xn, 2, ctr), 2, scl, `/`) %*% b)
       })
}
