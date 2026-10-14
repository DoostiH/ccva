#' Constraints on predictors and components
#'
#' Ecological knowledge enters a CCVA as explicit constraints, at two levels:
#'
#' * **Predictor level** -- the role of an indicator (or outcome proxy)
#'   within its component: the sign of its loading or weight, bounds,
#'   equality with other indicators, or a fixed value (0 drops it).
#' * **Component level** -- how a component enters vulnerability: exposure
#'   and sensitivity raise it, adaptive capacity lowers it; bounds on the
#'   effects; or adaptive capacity acting only as a *moderator* of exposure
#'   and sensitivity, with no main effect.
#' * **Between components** -- the three components are not independent.
#'   [constrain_link()] declares a directed dependency, e.g. dispersal-driven
#'   adaptive capacity shaping SDM-derived exposure. In the latent SEM it is
#'   a structural path, and the indirect and total effects of the source
#'   component on vulnerability are reported; in the other models it is
#'   resolved in two stages, by removing from the target component the part
#'   explained by the source, so it is not counted twice.
#'
#' Constraints are built from the `constrain_*()` helpers and collected with
#' `ccva_constraints()`. The same object is honoured by every combination
#' method of [ccva_combine()] that can express it: sign, bound, equality and
#' fixed constraints become inequality / equality restrictions on loadings
#' and paths in the latent SEM; signs become monotone-increasing or
#' -decreasing splines in the shape-constrained model and coefficient limits
#' in the penalised model; the moderator role becomes an interaction surface
#' in the shape and weighted models. A constraint a method cannot express is
#' reported, never silently ignored.
#'
#' Sign constraints on loadings also solve the sign-indeterminacy of latent
#' factors: with every loading pinned to its ecological direction, the
#' latent scores cannot come out reversed.
#'
#' @param ... `ccva_rule` objects from the `constrain_*()` helpers.
#' @param defaults Logical. If `TRUE` (default), start from
#'   [ccva_default_constraints()] (component signs plus indicator signs
#'   taken from the data's `direction`) and add `...` on top; later rules
#'   for the same target override earlier ones.
#'
#' @return An object of class `ccva_constraints`.
#' @seealso [ccva_combine()], [ccva_register()]
#' @export
#' @examples
#' # Default ecological constraints
#' ccva_constraints()
#'
#' # Adaptive capacity as moderator only, thermal tolerance loading fixed
#' # equal to range size, and generation length dropped
#' ccva_constraints(
#'   constrain_moderator(),
#'   constrain_equal(c("s1_thermal", "s2_range")),
#'   constrain_fixed("s10_genlen", 0)
#' )
#'
#' # No constraints at all (for comparison)
#' ccva_unconstrained()
ccva_constraints <- function(..., defaults = TRUE) {
  rules <- list(...)
  ok <- vapply(rules, inherits, logical(1), "ccva_rule")
  if (!all(ok))
    stop("All arguments to ccva_constraints() must come from constrain_*().",
         call. = FALSE)
  base <- if (defaults) ccva_default_constraints()$rules else list()
  structure(list(rules = c(base, rules)), class = "ccva_constraints")
}

#' @rdname ccva_constraints
#' @export
ccva_default_constraints <- function() {
  structure(list(rules = list(
    constrain_sign("exposure", "+"),
    constrain_sign("sensitivity", "+"),
    constrain_sign("adaptive", "-"),
    constrain_sign(".indicators", "direction")
  )), class = "ccva_constraints")
}

#' @rdname ccva_constraints
#' @export
ccva_unconstrained <- function() {
  structure(list(rules = list()), class = "ccva_constraints")
}

#' @param target For component-level rules, one of `"exposure"`,
#'   `"sensitivity"`, `"adaptive"`. For predictor-level rules, an indicator
#'   or outcome column name. `".indicators"` applies to every indicator and
#'   outcome.
#' @param sign `"+"`, `"-"`, or for `".indicators"` the special value
#'   `"direction"` (use each column's direction from [as_ccva()]).
#' @rdname ccva_constraints
#' @export
constrain_sign <- function(target, sign = c("+", "-", "direction")) {
  sign <- match.arg(sign)
  .rule("sign", target, sign = sign)
}

#' @param lower,upper Bounds on the effect (component level) or loading /
#'   weight (predictor level). Use `-Inf` / `Inf` for one-sided bounds.
#' @rdname ccva_constraints
#' @export
constrain_bounds <- function(target, lower = -Inf, upper = Inf) {
  if (lower >= upper) stop("`lower` must be below `upper`.", call. = FALSE)
  .rule("bounds", target, lower = lower, upper = upper)
}

#' @param targets Character vector of two or more indicators (within the
#'   same component) whose loadings / weights are constrained equal.
#' @rdname ccva_constraints
#' @export
constrain_equal <- function(targets) {
  if (length(targets) < 2) stop("Need at least two targets.", call. = FALSE)
  .rule("equal", targets)
}

#' @param value Fixed value. `0` removes the indicator or component effect.
#' @rdname ccva_constraints
#' @export
constrain_fixed <- function(target, value) {
  .rule("fixed", target, value = value)
}

#' @param from,to Components for a directed structural link: `from`
#'   affects `to` (e.g. adaptive capacity through dispersal shapes the
#'   SDM-derived exposure a species will face).
#' @param link_sign Sign of the link: `"free"`, `"+"` or `"-"`.
#' @rdname ccva_constraints
#' @export
constrain_link <- function(from, to, link_sign = c("free", "+", "-")) {
  link_sign <- match.arg(link_sign)
  ok <- c(from, to) %in% .ccva_components
  if (!all(ok) || length(from) != 1 || length(to) != 1)
    stop("`from` and `to` must each be one of: ",
         paste(.ccva_components, collapse = ", "), call. = FALSE)
  if (from == to) stop("A component cannot be linked to itself.",
                       call. = FALSE)
  .rule("link", c(from, to), sign = link_sign)
}

#' @param shape For the shape-constrained model: `"concave"` (e.g. an
#'   exposure effect that rises steeply then plateaus) or `"convex"`
#'   (accelerating, e.g. highly sensitive species disproportionately
#'   affected). Combined with the component's sign this gives monotone
#'   concave / convex splines. Unimodal shapes are not available.
#' @rdname ccva_constraints
#' @export
constrain_shape <- function(target, shape = c("concave", "convex")) {
  shape <- match.arg(shape)
  if (!target %in% .ccva_components)
    stop("constrain_shape() applies to a component: ",
         paste(.ccva_components, collapse = ", "), call. = FALSE)
  .rule("shape", target, shape = shape)
}

#' @param gamma Moderation strengths `c(exposure, sensitivity)` used by the
#'   weighted method: adaptive capacity reduces the exposure and sensitivity
#'   effects by up to these proportions.
#' @rdname ccva_constraints
#' @export
constrain_moderator <- function(gamma = c(0.5, 0.4)) {
  .rule("moderator", "adaptive", gamma = gamma)
}

.rule <- function(type, target, ...) {
  structure(list(type = type, target = target, ...), class = "ccva_rule")
}

#' @export
print.ccva_constraints <- function(x, ...) {
  if (!length(x$rules)) {
    cat("<ccva_constraints> none (free estimation)\n")
    return(invisible(x))
  }
  cat("<ccva_constraints> ", length(x$rules), " rule(s)\n", sep = "")
  for (r in x$rules) cat("  ", format(r), "\n", sep = "")
  invisible(x)
}

#' @export
format.ccva_rule <- function(x, ...) {
  tgt <- if (identical(x$target, ".indicators")) "all indicators & outcomes"
         else paste(x$target, collapse = " = ")
  if (x$type == "link") tgt <- paste(x$target, collapse = " -> ")
  lvl <- if (all(x$target %in% .ccva_components)) "component" else "predictor"
  switch(x$type,
    sign = if (x$sign == "direction")
      sprintf("[%s] %s: sign follows declared direction", lvl, tgt)
      else sprintf("[%s] %s: %s effect", lvl, tgt,
                   if (x$sign == "+") "positive" else "negative"),
    bounds = sprintf("[%s] %s in [%s, %s]", lvl, tgt, x$lower, x$upper),
    equal = sprintf("[predictor] equal: %s", tgt),
    fixed = sprintf("[%s] %s fixed at %s%s", lvl, tgt, x$value,
                    if (x$value == 0) " (dropped)" else ""),
    moderator = "[component] adaptive capacity acts only as a moderator",
    shape = sprintf("[component] %s: %s shape", tgt, x$shape),
    link = sprintf("[link] %s -> %s%s", x$target[1], x$target[2],
                   switch(x$sign, free = "", "+" = " (positive)",
                          "-" = " (negative)")),
    sprintf("[%s] %s", x$type, tgt))
}

## Resolve a constraints object against a dataset ---------------------------
# Returns a list the backends can consume directly.
.resolve_constraints <- function(constraints, x) {
  if (is.null(constraints)) constraints <- ccva_default_constraints()
  if (!inherits(constraints, "ccva_constraints"))
    stop("`constraints` must come from ccva_constraints().", call. = FALSE)
  comps <- .active_components(x)
  cols <- c(unlist(x$indicators), x$outcomes)

  comp_sign <- stats::setNames(rep(NA_real_, 3), .ccva_components)
  comp_bounds <- list()
  comp_fixed <- list()
  ind_sign <- stats::setNames(rep(NA_real_, length(cols)), cols)
  ind_bounds <- list()
  ind_fixed <- list()
  equal <- list()
  moderator <- NULL
  links <- list()
  comp_shape <- list()

  for (r in constraints$rules) {
    tg <- r$target
    is_comp <- all(tg %in% .ccva_components) || r$type == "link"
    if (!is_comp && !identical(tg, ".indicators")) {
      bad <- setdiff(tg, cols)
      if (length(bad))
        stop("Constraint target(s) not in the data: ",
             paste(bad, collapse = ", "), call. = FALSE)
    }
    switch(r$type,
      sign = {
        if (identical(tg, ".indicators")) {
          # internally every column is oriented, so "follow direction"
          # is +1 on the oriented scale
          ind_sign[] <- if (r$sign == "direction") rep(1, length(cols)) else
            (if (r$sign == "+") 1 else -1) * unname(x$direction[cols])
        } else if (is_comp) {
          comp_sign[tg] <- if (r$sign == "+") 1 else -1
        } else {
          # predictor sign is stated on the raw scale; convert to the
          # oriented scale used internally
          ind_sign[tg] <- (if (r$sign == "+") 1 else -1) * x$direction[tg]
        }
      },
      bounds = {
        if (is_comp) comp_bounds[[tg]] <- c(r$lower, r$upper)
        else ind_bounds[[tg]] <- c(r$lower, r$upper)
      },
      fixed = {
        if (is_comp) comp_fixed[[tg]] <- r$value
        else ind_fixed[[tg]] <- r$value
      },
      equal = {
        grp <- vapply(tg, .component_of, character(1), x = x)
        if (length(unique(grp)) != 1)
          stop("constrain_equal() targets must belong to one component.",
               call. = FALSE)
        equal[[length(equal) + 1]] <- tg
      },
      shape = { comp_shape[[tg]] <- r$shape },
      link = {
        miss <- setdiff(tg, comps)
        if (length(miss))
          stop("constrain_link(): no indicators for ",
               paste(miss, collapse = ", "), call. = FALSE)
        for (l in links)
          if (identical(l$from, tg[2]) && identical(l$to, tg[1]))
            stop("Links ", tg[1], " -> ", tg[2], " and the reverse form a ",
                 "cycle.", call. = FALSE)
        links[[paste(tg, collapse = "->")]] <-
          list(from = tg[1], to = tg[2], sign = r$sign)
      },
      moderator = {
        if (!"adaptive" %in% comps)
          stop("constrain_moderator() needs adaptive-capacity indicators.",
               call. = FALSE)
        moderator <- r$gamma
      })
  }
  list(comp_sign = comp_sign[comps], comp_bounds = comp_bounds,
       comp_fixed = comp_fixed, ind_sign = ind_sign,
       ind_bounds = ind_bounds, ind_fixed = ind_fixed, equal = equal,
       moderator = moderator, links = unname(links),
       comp_shape = comp_shape,
       n_rules = length(constraints$rules),
       object = constraints)
}

.component_of <- function(col, x) {
  for (cmp in .ccva_components)
    if (col %in% x$indicators[[cmp]]) return(cmp)
  if (col %in% x$outcomes) return("outcome")
  NA_character_
}

.active_components <- function(x) {
  .ccva_components[vapply(x$indicators[.ccva_components], length,
                          integer(1)) > 0]
}

## Method registry ------------------------------------------------------------

.ccva_registry <- new.env(parent = emptyenv())
.ccva_registry$component <- list()
.ccva_registry$combine <- list()

#' Register a new component or combination method
#'
#' The framework grows by contribution. A registered method is used exactly
#' like a built-in one, through the `method` argument of the component
#' scorers or of [ccva_combine()].
#'
#' A **component** method is a function `function(X, anchor, signs, ...)`
#' where `X` is the oriented indicator matrix (higher = more of the
#' component), `anchor` an optional numeric vector (oriented vulnerability
#' proxy) and `signs` the resolved indicator sign constraints. It must
#' return a list with `score` (numeric, one per row) and `predict`, a
#' function of a new `X` returning scores on the same scale.
#'
#' A **combine** method is a function `function(x, components, constraints,
#' ...)` receiving the `ccva_data`, the fitted component scores (data frame)
#' and the resolved constraints. It must return a list with `score`
#' (vulnerability, higher = more vulnerable) and optionally `model`,
#' `details` and `predict` (a function of new `ccva_data`).
#'
#' @param name Method name.
#' @param fun The method function.
#' @param type `"component"` or `"combine"`.
#' @param description One-line description shown by [ccva_methods()].
#' @return Invisibly, the updated registry listing.
#' @export
#' @examples
#' # A median-based component score
#' ccva_register("median", function(X, anchor, signs, ...) {
#'   z <- scale(X)
#'   list(score = apply(z, 1, stats::median),
#'        predict = function(newX) apply(scale(newX,
#'          attr(z, "scaled:center"), attr(z, "scaled:scale")), 1,
#'          stats::median))
#' }, type = "component", description = "Median of standardised indicators")
#' ccva_methods()
#' ccva_unregister("median", "component")
ccva_register <- function(name, fun, type = c("component", "combine"),
                          description = "") {
  type <- match.arg(type)
  if (!is.function(fun)) stop("`fun` must be a function.", call. = FALSE)
  builtin <- if (type == "component") .builtin_component else .builtin_combine
  if (name %in% builtin)
    stop("'", name, "' is a built-in method name.", call. = FALSE)
  reg <- .ccva_registry[[type]]
  reg[[name]] <- list(fun = fun, description = description)
  .ccva_registry[[type]] <- reg
  invisible(ccva_methods())
}

#' @rdname ccva_register
#' @export
ccva_unregister <- function(name, type = c("component", "combine")) {
  type <- match.arg(type)
  reg <- .ccva_registry[[type]]
  reg[[name]] <- NULL
  .ccva_registry[[type]] <- reg
  invisible(ccva_methods())
}

#' List available methods
#'
#' @return A data frame of built-in and registered methods.
#' @export
#' @examples
#' ccva_methods()
ccva_methods <- function() {
  b <- data.frame(
    type = c(rep("component", length(.builtin_component)),
             rep("combine", length(.builtin_combine))),
    method = c(.builtin_component, .builtin_combine),
    source = "built-in",
    description = c(
      "Ofori-style 1-3 scores summed (optionally weighted)",
      "Weighted mean of indicators rescaled to [0, 1]",
      "First principal component",
      "LASSO / ridge / elastic net against an anchor (p >= n)",
      "Weighted index (Ofori et al. 2017), corrected normalisation",
      "Shape-constrained additive model (scam)",
      "Latent-variable SEM with constraints (lavaan)",
      "Two-stage: latent V*, then penalised regression on indicators"),
    stringsAsFactors = FALSE)
  reg <- lapply(c("component", "combine"), function(tp) {
    r <- .ccva_registry[[tp]]
    if (!length(r)) return(NULL)
    data.frame(type = tp, method = names(r), source = "registered",
               description = vapply(r, `[[`, character(1), "description"),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, c(list(b), reg))
  rownames(out) <- NULL
  out
}

.builtin_component <- c("ordinal", "mean", "pca", "penalised")
.builtin_combine <- c("weighted", "shape", "latent", "penalised")
