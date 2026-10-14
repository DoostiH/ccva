#' Plot any step of a CCVA
#'
#' One entry point for the figures of every stage. `ccva_plot()` dispatches
#' on the class of `object`; `type` picks the figure. Every plot is a
#' `ggplot` object, so it can be themed, combined or saved with
#' `ggplot2::ggsave()`.
#'
#' | object | `type` | shows |
#' |---|---|---|
#' | `ccva_data` | `"indicators"` | distribution of every indicator, grouped by component |
#' | | `"correlation"` | indicator correlation matrix, ordered by component |
#' | `ccva_check` | `"alpha"` | Cronbach's alpha per component against thresholds |
#' | | `"correlation"` | cross-component correlations |
#' | `ccva_component` | `"importance"` | weight / loading / coefficient of each indicator |
#' | | `"distribution"` | distribution of the component score |
#' | `ccva_fit` | `"structure"` | path diagram: components to vulnerability, with effects |
#' | | `"components"` | the three component scores side by side |
#' | | `"scores"` | vulnerability distribution with class cut-offs |
#' | | `"ranking"` | most and least vulnerable units with their components |
#' | | `"contributions"` | how each component pushes each unit's score (weighted) |
#' | | `"effects"` | shape curves, latent loadings or penalised coefficients |
#' | | `"map"` | vulnerability over coordinates |
#' | `ccva_comparison` | `"truth"` | each method against the true vulnerability |
#' | | `"ranks"` | rank of each unit under two methods |
#' | | `"validation"` | correlation of each method with each target |
#' | `ccva_uncertainty` | `"intervals"` | score intervals, ranked |
#' | | `"pairwise"` | probability one unit is more vulnerable than another |
#'
#' @param object A `ccva_data`, `ccva_check`, `ccva_component`, `ccva_fit`,
#'   `ccva_comparison` or `ccva_uncertainty` object.
#' @param type Figure type (see table). `NULL` gives the first listed.
#' @param n For rankings and interval plots, the number of units shown at
#'   each end (`"ranking"`) or in total (`"intervals"`, `"pairwise"`).
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
#' @examples
#' d <- ccva_sim_data()
#' ccva_plot(d, "indicators")
#' fit <- ccva_combine(d, "weighted", component_method = "ordinal")
#' ccva_plot(fit, "structure")
#' ccva_plot(fit, "ranking", n = 8)
ccva_plot <- function(object, type = NULL, n = 15, ...) {
  UseMethod("ccva_plot")
}

#' @export
ccva_plot.default <- function(object, type = NULL, n = 15, ...) {
  stop("No ccva_plot() method for class '", class(object)[1], "'.",
       call. = FALSE)
}

#' @export
plot.ccva_data <- function(x, ...) ccva_plot(x, ...)
#' @export
plot.ccva_check <- function(x, ...) ccva_plot(x, ...)
#' @export
plot.ccva_component <- function(x, ...) ccva_plot(x, ...)
#' @export
plot.ccva_fit <- function(x, ...) ccva_plot(x, ...)
#' @export
plot.ccva_comparison <- function(x, ...) ccva_plot(x, ...)
#' @export
plot.ccva_uncertainty <- function(x, ...) ccva_plot(x, ...)

## Palette and theme -----------------------------------------------------------

#' Colours and theme used by ccva plots
#'
#' `ccva_colours` maps each component (and vulnerability) to a fixed colour,
#' so a component looks the same in every figure. The palette was checked
#' for colour-vision-deficiency separation; facet and direct labels carry
#' identity too, so colour is never the only cue.
#'
#' @param base_size Base font size.
#' @return `theme_ccva()` returns a ggplot2 theme.
#' @export
#' @examples
#' ccva_colours
theme_ccva <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#e6e6e3",
                                               linewidth = 0.3),
      axis.ticks = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0,
                                         colour = "#33332f"),
      plot.title = ggplot2::element_text(face = "bold", colour = "#1f1f1c"),
      plot.subtitle = ggplot2::element_text(colour = "#5c5c57"),
      plot.caption = ggplot2::element_text(colour = "#7a7a74", hjust = 0),
      axis.title = ggplot2::element_text(colour = "#45453f"),
      axis.text = ggplot2::element_text(colour = "#5c5c57"),
      legend.position = "top", legend.justification = "left",
      legend.title = ggplot2::element_text(colour = "#45453f"),
      plot.title.position = "plot", plot.caption.position = "plot")
}

#' @rdname theme_ccva
#' @format `ccva_colours`: a named character vector.
#' @export
ccva_colours <- c(exposure = "#2a6fb0", sensitivity = "#c4573a",
                  adaptive = "#2f8f6b", vulnerability = "#8a5aa8",
                  outcome = "#6b6b66")

.class_colours <- c(Low = "#b9d3ea", Moderate = "#5f93c6",
                    High = "#1f4f80")

.comp_scale <- function(aes = "colour") {
  ggplot2::scale_discrete_manual(aes, values = ccva_colours,
                                 labels = function(b) unname(.lab(b)),
                                 name = NULL)
}

.comp_factor <- function(v) factor(v, levels = c(.ccva_components,
                                                 "outcome"),
                                   labels = c(.ccva_labels[1:3], "Outcome"))

## ccva_data ---------------------------------------------------------------------

#' @export
ccva_plot.ccva_data <- function(object, type = NULL, n = 15, ...) {
  type <- match.arg(type %||% "indicators", c("indicators", "correlation"))
  x <- object
  comps <- .active_components(x)
  vars <- unlist(x$indicators[comps])
  if (length(vars) > 60) {
    vars <- unlist(lapply(comps, function(cmp)
      utils::head(x$indicators[[cmp]], 20)))
    cap <- "Showing up to 20 indicators per component."
  } else cap <- NULL
  grp <- vapply(vars, .component_of, character(1), x = x)
  if (type == "indicators") {
    Z <- scale(sweep(as.matrix(x$data[vars]), 2, x$direction[vars], `*`))
    long <- data.frame(indicator = rep(vars, each = nrow(Z)),
                       component = rep(grp, each = nrow(Z)),
                       value = as.numeric(Z))
    long$indicator <- factor(long$indicator, levels = rev(vars))
    long$component <- factor(long$component, levels = comps)
    return(
      ggplot2::ggplot(long, ggplot2::aes(value, indicator, fill = component)) +
        ggplot2::geom_violin(colour = NA, alpha = 0.35, scale = "width") +
        ggplot2::geom_boxplot(width = 0.25, outlier.size = 0.4,
                              colour = "#33332f", alpha = 0.9) +
        ggplot2::facet_grid(component ~ ., scales = "free_y", space = "free_y",
                            labeller = ggplot2::as_labeller(.ccva_labels)) +
        .comp_scale("fill") +
        ggplot2::labs(title = "Indicators by component",
                      subtitle = "Standardised and oriented: right = more of the component",
                      x = "Standardised value", y = NULL, caption = cap) +
        theme_ccva() + ggplot2::theme(legend.position = "none"))
  }
  allv <- c(vars, x$outcomes)
  M <- sweep(as.matrix(x$data[allv]), 2, x$direction[allv], `*`)
  R <- stats::cor(M, method = "spearman")
  lv <- allv
  df <- data.frame(var1 = factor(rep(lv, length(lv)), levels = lv),
                   var2 = factor(rep(lv, each = length(lv)), levels = rev(lv)),
                   r = as.numeric(R))
  ggplot2::ggplot(df, ggplot2::aes(var1, var2, fill = r)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::scale_fill_gradient2(low = "#c4573a", mid = "#f4f4f1",
                                  high = "#2a6fb0", limits = c(-1, 1),
                                  name = "Spearman r") +
    ggplot2::labs(title = "Correlations among indicators and outcomes",
                  subtitle = "Oriented: positive = both point towards more of their component / more vulnerability",
                  x = NULL, y = NULL, caption = cap) +
    ggplot2::coord_equal() + theme_ccva() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                       hjust = 1, vjust = 0.5),
                   panel.grid = ggplot2::element_blank())
}

## ccva_check ----------------------------------------------------------------------

#' @export
ccva_plot.ccva_check <- function(object, type = NULL, n = 15, ...) {
  type <- match.arg(type %||% "alpha", c("alpha", "correlation"))
  if (type == "alpha") {
    a <- object$alpha
    a$component <- factor(a$component, levels = rev(a$component))
    return(
      ggplot2::ggplot(a, ggplot2::aes(alpha, component)) +
        ggplot2::annotate("rect", xmin = -Inf, xmax = 0.6, ymin = -Inf,
                          ymax = Inf, fill = "#f6e3dc", alpha = 0.6) +
        ggplot2::geom_vline(xintercept = c(0.6, 0.7, 0.8), linetype = 2,
                            colour = "#9a9a94", linewidth = 0.3) +
        ggplot2::geom_segment(ggplot2::aes(x = 0, xend = alpha,
                                           yend = component,
                                           colour = as.character(component)),
                              linewidth = 1.2, na.rm = TRUE) +
        ggplot2::geom_point(ggplot2::aes(colour = as.character(component)),
                            size = 4, na.rm = TRUE) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", alpha)),
                           nudge_y = 0.3, size = 3.4, colour = "#33332f",
                           na.rm = TRUE) +
        .comp_scale() +
        ggplot2::scale_y_discrete(labels = .ccva_labels) +
        ggplot2::scale_x_continuous(limits = c(min(0, a$alpha, na.rm = TRUE), 1)) +
        ggplot2::labs(title = "Do the indicators of each component agree?",
                      subtitle = "Cronbach's alpha; shaded = questionable (< 0.6), dashed lines at 0.6 / 0.7 / 0.8",
                      x = "Cronbach's alpha", y = NULL) +
        theme_ccva() + ggplot2::theme(legend.position = "none"))
  }
  .cor_heat(object$component_cor, "Correlation between components",
            "Spearman; strong correlation means the additive index double-counts")
}

.cor_heat <- function(R, title, subtitle) {
  lv <- colnames(R)
  df <- data.frame(var1 = factor(rep(lv, length(lv)), levels = lv),
                   var2 = factor(rep(lv, each = length(lv)), levels = rev(lv)),
                   r = as.numeric(R))
  ggplot2::ggplot(df, ggplot2::aes(var1, var2, fill = r)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", r)), size = 4,
                       colour = "#1f1f1c") +
    ggplot2::scale_fill_gradient2(low = "#c4573a", mid = "#f4f4f1",
                                  high = "#2a6fb0", limits = c(-1, 1),
                                  name = "r") +
    ggplot2::scale_x_discrete(labels = function(v) .lab(v)) +
    ggplot2::scale_y_discrete(labels = function(v) .lab(v)) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::coord_equal() + theme_ccva() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

.lab <- function(v) {
  v <- as.character(v)
  ifelse(v %in% names(.ccva_labels), .ccva_labels[v], v)
}

## ccva_component -------------------------------------------------------------------

#' @export
ccva_plot.ccva_component <- function(object, type = NULL, n = 15, ...) {
  type <- match.arg(type %||% "importance", c("importance", "distribution"))
  col <- ccva_colours[[object$component]]
  if (type == "distribution") {
    df <- data.frame(score = object$score)
    return(
      ggplot2::ggplot(df, ggplot2::aes(score)) +
        ggplot2::geom_histogram(bins = 30, fill = col, colour = "white",
                                linewidth = 0.3) +
        ggplot2::labs(title = sprintf("%s score", .ccva_labels[[object$component]]),
                      subtitle = sprintf("Method: %s", object$method),
                      x = "Component score", y = "Units") +
        theme_ccva())
  }
  imp <- object$importance
  df <- data.frame(indicator = names(imp), value = as.numeric(imp))
  df <- df[order(-abs(df$value)), ]
  sub <- "Weight, loading or coefficient of each indicator"
  if (nrow(df) > 40) {
    nz <- sum(df$value != 0)
    df <- utils::head(df, 40)
    sub <- sprintf("Top 40 of %d indicators (%d non-zero)",
                   length(imp), nz)
  }
  df$indicator <- factor(df$indicator, levels = rev(df$indicator))
  ggplot2::ggplot(df, ggplot2::aes(value, indicator)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#9a9a94") +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = value,
                                       yend = indicator),
                          colour = col, linewidth = 1) +
    ggplot2::geom_point(colour = col, size = 2.8) +
    ggplot2::labs(title = sprintf("What drives %s?",
                                  tolower(.ccva_labels[[object$component]])),
                  subtitle = sprintf("%s (method: %s)", sub, object$method),
                  x = NULL, y = NULL) +
    theme_ccva()
}

## ccva_fit -----------------------------------------------------------------------

#' @export
ccva_plot.ccva_fit <- function(object, type = NULL, n = 15, ...) {
  type <- match.arg(type %||% "structure",
                    c("structure", "components", "scores", "ranking",
                      "contributions", "effects", "map"))
  switch(type,
         structure = .plot_structure(object),
         components = .plot_components(object),
         scores = .plot_scores(object),
         ranking = .plot_ranking(object, n),
         contributions = .plot_contributions(object, n),
         effects = .plot_effects(object),
         map = .plot_map(object))
}

.method_label <- function(fit) {
  c(weighted = "Weighted index", shape = "Shape-constrained model",
    latent = "Latent SEM", penalised = "Penalised two-stage")[[fit$method]] %||%
    fit$method
}

.plot_structure <- function(fit) {
  comps <- names(fit$components)
  ypos <- c(exposure = 3, sensitivity = 2, adaptive = 1)[comps]
  nodes <- data.frame(id = c(comps, "vulnerability"),
                      x = c(rep(1, length(comps)), 3.2),
                      y = c(ypos, mean(ypos)))
  nodes$label <- c(sprintf("%s\n%d indicator%s", .ccva_labels[comps],
                           vapply(fit$components, function(cp)
                             length(cp$indicators), integer(1)),
                           ifelse(vapply(fit$components, function(cp)
                             length(cp$indicators), integer(1)) == 1, "", "s")),
                   "Vulnerability")
  # edge labels
  elab <- switch(fit$method,
    weighted = {
      w <- fit$details$weights; s <- fit$details$signs
      sprintf("%s weight %.2f", ifelse(s[comps] > 0, "+", "-"),
              w[comps] / sum(w))
    },
    latent = {
      e <- fit$details$effects
      v <- e$std[match(comps, e$component)]
      ifelse(is.na(v), "fixed", sprintf("effect %+.2f", v))
    },
    shape = {
      s <- fit$constraints$comp_sign[comps]
      if (isTRUE(fit$details$moderator))
        ifelse(comps == "adaptive", "moderates", "increasing")
      else ifelse(s > 0, "increasing", "decreasing")
    },
    penalised = {
      cf <- fit$details$coefficients
      vapply(comps, function(cmp) {
        b <- cf$coef[cf$component == cmp]
        sprintf("%d of %d selected", sum(b != 0), length(b))
      }, character(1))
    },
    rep("", length(comps)))
  edges <- data.frame(component = comps, x = 1.55, y = ypos,
                      xend = 2.6, yend = mean(ypos), label = elab)
  # labels sit on each edge near its source, where edges do not overlap
  edges$lx <- edges$x + 0.32 * (edges$xend - edges$x)
  edges$ly <- edges$y + 0.32 * (edges$yend - edges$y)
  cons <- fit$constraints
  cap <- sprintf("%s | %d constraint rule(s)%s", .method_label(fit),
                 cons$n_rules,
                 if (!is.null(fit$details$fit_measures))
                   sprintf(" | CFI %.2f, RMSEA %.3f",
                           fit$details$fit_measures[["cfi"]],
                           fit$details$fit_measures[["rmsea"]]) else "")
  te <- fit$details$total_effects
  if (!is.null(te))
    cap <- paste0(cap, " | total effect: ",
                  paste(sprintf("%s %+.2f (direct %+.2f, indirect %+.2f)",
                                te$component, te$total, te$direct,
                                te$indirect), collapse = "; "))
  mod <- !is.null(cons$moderator) && "adaptive" %in% comps
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = if (mod) edges[edges$component != "adaptive", ]
                          else edges,
                          ggplot2::aes(x, y, xend = xend, yend = yend,
                                       colour = component),
                          linewidth = 1.1,
                          arrow = ggplot2::arrow(length = ggplot2::unit(3, "mm"),
                                                 type = "closed"))
  if (mod) {
    ea <- edges[edges$component == "adaptive", ]
    tgt <- edges[edges$component != "adaptive", ]
    p <- p + ggplot2::geom_curve(
      data = data.frame(x = ea$x, y = ea$y,
                        xend = tgt$x + 0.62 * (tgt$xend - tgt$x),
                        yend = tgt$y + 0.62 * (tgt$yend - tgt$y) - 0.04),
      ggplot2::aes(x, y, xend = xend, yend = yend), colour = ccva_colours[["adaptive"]],
      linetype = 2, curvature = 0.25, linewidth = 0.9,
      arrow = ggplot2::arrow(length = ggplot2::unit(2.5, "mm")))
  }
  lk <- fit$details$links
  linked <- if (is.null(lk)) character(0) else
    apply(lk[c("from", "to")], 1, function(v) paste(sort(v), collapse = "|"))
  if (length(comps) >= 2 && fit$method == "latent") {
    cy <- sort(ypos)
    pr <- data.frame(a = names(sort(ypos))[-length(cy)],
                     b = names(sort(ypos))[-1])
    pr <- pr[!apply(pr, 1, function(v) paste(sort(v), collapse = "|")) %in%
               linked, , drop = FALSE]
    if (nrow(pr))
      p <- p + ggplot2::geom_curve(
        data = data.frame(x = 0.42, y = ypos[pr$a] + 0.08,
                          xend = 0.42, yend = ypos[pr$b] - 0.08),
        ggplot2::aes(x, y, xend = xend, yend = yend), curvature = 0.5,
        colour = "#9a9a94", linewidth = 0.5,
        arrow = ggplot2::arrow(length = ggplot2::unit(2, "mm"),
                               ends = "both"))
  }
  # directed links between components (causal diagram)
  if (!is.null(lk) && nrow(lk)) {
    ld <- data.frame(from = lk$from, to = lk$to,
                     y = ypos[lk$from], yend = ypos[lk$to])
    ld$dist <- abs(ld$y - ld$yend)
    ld$x <- 0.46 - 0.05 * (ld$dist - 1)
    ld$lab <- if (lk$method[1] == "structural path")
      sprintf("%+.2f", lk$std) else "two-stage"
    ld$lx <- ld$x - 0.2 - 0.12 * (ld$dist - 1)
    ld$ly <- (ld$y + ld$yend) / 2
    for (i in seq_len(nrow(ld))) {
      di <- ld[i, ]
      dir <- if (di$yend > di$y) -1 else 1
      p <- p + ggplot2::geom_curve(
        data = data.frame(x = di$x, y = di$y + 0.15 * sign(di$yend - di$y),
                          xend = di$x, yend = di$yend - 0.15 * sign(di$yend - di$y)),
        ggplot2::aes(x, y, xend = xend, yend = yend),
        curvature = 0.45 * dir, colour = ccva_colours[[di$from]],
        linewidth = 1, linetype = if (di$lab == "two-stage") 2 else 1,
        arrow = ggplot2::arrow(length = ggplot2::unit(2.5, "mm"),
                               type = "closed"))
    }
    p <- p + ggplot2::geom_label(data = ld, ggplot2::aes(lx, ly, label = lab),
                                 size = 3, label.size = 0, fill = "white",
                                 colour = "#33332f")
  }
  p + ggplot2::geom_label(data = edges[!(mod & edges$component == "adaptive"), ],
                          ggplot2::aes(lx, ly, label = label),
                          size = 3.3, label.size = 0, fill = "white",
                          colour = "#33332f") +
    ggplot2::geom_label(data = nodes, ggplot2::aes(x, y, label = label,
                                                   fill = id),
                        colour = "white", size = 3.8, fontface = "bold",
                        label.padding = ggplot2::unit(0.55, "lines"),
                        label.size = 0) +
    ggplot2::scale_colour_manual(values = ccva_colours, guide = "none") +
    ggplot2::scale_fill_manual(values = ccva_colours, guide = "none") +
    ggplot2::scale_x_continuous(limits = c(if (is.null(lk)) 0.3 else -0.05,
                                           3.75)) +
    ggplot2::scale_y_continuous(limits = c(0.5, 3.5)) +
    ggplot2::labs(title = "How the components combine into vulnerability",
                  subtitle = if (!is.null(lk)) "Causal diagram: directed links between components (left) and effects on vulnerability"
                  else if (mod) "Adaptive capacity acts only as a moderator (dashed)"
                  else if (fit$method == "latent")
                    "Latent components (correlated, grey) with standardised effects on latent vulnerability"
                  else NULL,
                  caption = cap) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   plot.subtitle = ggplot2::element_text(colour = "#5c5c57"),
                   plot.caption = ggplot2::element_text(colour = "#7a7a74",
                                                        hjust = 0),
                   plot.margin = ggplot2::margin(10, 10, 10, 10))
}

.plot_components <- function(fit) {
  cs <- fit$component_scores
  long <- data.frame(component = rep(names(cs), each = nrow(cs)),
                     score = unlist(cs, use.names = FALSE))
  long$component <- factor(long$component, levels = names(cs))
  meth <- vapply(fit$components, `[[`, character(1), "method")
  ggplot2::ggplot(long, ggplot2::aes(score, fill = component)) +
    ggplot2::geom_histogram(bins = 30, colour = "white", linewidth = 0.3) +
    ggplot2::facet_wrap(~component, nrow = 1,
                        labeller = ggplot2::as_labeller(.ccva_labels)) +
    .comp_scale("fill") +
    ggplot2::labs(title = "Component scores",
                  subtitle = paste("Scoring:", paste(sprintf("%s = %s",
                    .ccva_labels[names(meth)], meth), collapse = "; ")),
                  x = "Score (higher = more of the component)", y = "Units") +
    theme_ccva() + ggplot2::theme(legend.position = "none")
}

.plot_scores <- function(fit) {
  df <- data.frame(score = fit$score, class = fit$class)
  br <- attr(fit$class, "cutoffs")
  ggplot2::ggplot(df, ggplot2::aes(score, fill = class)) +
    ggplot2::geom_histogram(bins = 35, colour = "white", linewidth = 0.3,
                            boundary = 0) +
    ggplot2::geom_vline(xintercept = br, linetype = 2, colour = "#45453f") +
    ggplot2::scale_fill_manual(values = .class_colours, name = NULL,
                               drop = FALSE) +
    ggplot2::labs(title = "Distribution of vulnerability",
                  subtitle = sprintf("%s; class cut-offs (%s) at %.2f and %.2f",
                                     .method_label(fit),
                                     attr(fit$class, "method"), br[1], br[2]),
                  x = "Vulnerability (0 = least, 1 = most)", y = "Units") +
    theme_ccva()
}

.plot_ranking <- function(fit, n) {
  rk <- ccva_rank(fit)
  n <- min(n, floor(nrow(rk) / 2))
  sel <- rbind(utils::head(rk, n), utils::tail(rk, n))
  sel$group <- rep(c("Most vulnerable", "Least vulnerable"), each = n)
  sel$group <- factor(sel$group, levels = c("Most vulnerable",
                                            "Least vulnerable"))
  comps <- names(fit$components)
  long <- do.call(rbind, lapply(comps, function(cmp)
    data.frame(unit = sel$unit, group = sel$group, component = cmp,
               value = sel[[cmp]])))
  ord <- sel$unit[order(sel$score)]
  sel$unit <- factor(sel$unit, levels = ord)
  long$unit <- factor(long$unit, levels = ord)
  long$component <- factor(long$component, levels = comps)
  ggplot2::ggplot() +
    ggplot2::geom_point(data = long, ggplot2::aes(value, unit,
                                                  colour = component),
                        size = 2, alpha = 0.75, shape = 16) +
    ggplot2::geom_point(data = sel, ggplot2::aes(score, unit),
                        shape = 21, size = 3.6, fill = ccva_colours[["vulnerability"]],
                        colour = "white", stroke = 0.8) +
    ggplot2::facet_wrap(~group, scales = "free_y", ncol = 1) +
    ggplot2::scale_colour_manual(values = ccva_colours,
                                 labels = function(b) unname(.lab(b)),
                                 name = NULL) +
    ggplot2::labs(title = sprintf("The %d most and least vulnerable units", n),
                  subtitle = "Large purple dot = vulnerability; small dots = component scores",
                  x = "Score", y = NULL, caption = .method_label(fit)) +
    theme_ccva()
}

.plot_contributions <- function(fit, n) {
  if (fit$method != "weighted")
    stop("Contributions are defined for the weighted index; use ",
         "type = \"effects\" for other methods.", call. = FALSE)
  w <- fit$details$weights; s <- fit$details$signs
  cs <- fit$component_scores
  comps <- names(cs)
  ctr <- sweep(as.matrix(cs), 2, (w * s)[comps] / sum(w), `*`)
  rk <- order(-fit$score)
  n <- min(n, floor(nrow(cs) / 2))
  idx <- c(utils::head(rk, n), utils::tail(rk, n))
  long <- do.call(rbind, lapply(comps, function(cmp)
    data.frame(unit = rownames(cs)[idx], component = cmp,
               value = ctr[idx, cmp])))
  long$unit <- factor(long$unit, levels = rev(rownames(cs)[idx]))
  long$component <- factor(long$component, levels = comps)
  ggplot2::ggplot(long, ggplot2::aes(value, unit, fill = component)) +
    ggplot2::geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = 0, colour = "#45453f") +
    ggplot2::geom_hline(yintercept = n + 0.5, colour = "#9a9a94",
                        linetype = 3) +
    .comp_scale("fill") +
    ggplot2::labs(title = "What pushes each unit's vulnerability up or down",
                  subtitle = sprintf("Weighted contributions; top %d and bottom %d units", n, n),
                  x = "Contribution to the index", y = NULL) +
    theme_ccva()
}

.plot_effects <- function(fit) {
  switch(fit$method,
    shape = {
      cv <- fit$details$curves
      cv$component <- factor(cv$component, levels = names(fit$components))
      rug <- do.call(rbind, lapply(levels(cv$component), function(cmp)
        data.frame(component = cmp, x = fit$component_scores[[cmp]])))
      rug$component <- factor(rug$component, levels = levels(cv$component))
      ggplot2::ggplot(cv, ggplot2::aes(x, y, colour = component)) +
        ggplot2::geom_rug(data = rug, ggplot2::aes(x = x), inherit.aes = FALSE,
                          alpha = 0.25, colour = "#6b6b66") +
        ggplot2::geom_line(linewidth = 1.3) +
        ggplot2::facet_wrap(~component, nrow = 1, scales = "free_x",
                            labeller = ggplot2::as_labeller(.ccva_labels)) +
        .comp_scale() +
        ggplot2::labs(title = "Learned shape of each component's effect",
                      subtitle = if (isTRUE(fit$details$moderator))
                        "Interaction surfaces, sliced with other components at their mean"
                      else "Direction fixed by constraints, shape estimated from data (others at mean)",
                      x = "Component score", y = "Predicted response") +
        theme_ccva() + ggplot2::theme(legend.position = "none")
    },
    latent = {
      ld <- fit$details$loadings
      fmap <- c(stats::setNames(names(.ccva_latent), .ccva_latent),
                Vuln = "vulnerability")
      ld$component <- fmap[ld$factor]
      ld <- ld[!is.na(ld$component), ]
      ld$component <- factor(ld$component, levels = c(names(fit$components),
                                                      "vulnerability"))
      ld <- ld[order(ld$component, ld$std), ]
      ld$indicator <- factor(ld$indicator, levels = unique(ld$indicator))
      ggplot2::ggplot(ld, ggplot2::aes(std, indicator, colour = component)) +
        ggplot2::geom_vline(xintercept = 0, colour = "#9a9a94") +
        ggplot2::geom_segment(ggplot2::aes(x = 0, xend = std,
                                           yend = indicator), linewidth = 1) +
        ggplot2::geom_point(size = 2.6) +
        ggplot2::facet_grid(component ~ ., scales = "free_y", space = "free_y",
                            labeller = ggplot2::as_labeller(.ccva_labels)) +
        .comp_scale() +
        ggplot2::labs(title = "Measurement model: how each indicator reflects its latent component",
                      subtitle = "Standardised loadings on the oriented scale (sign-constrained by default)",
                      x = "Standardised loading", y = NULL) +
        theme_ccva() + ggplot2::theme(legend.position = "none")
    },
    penalised = {
      cf <- fit$details$coefficients
      cf <- cf[cf$coef != 0, ]
      if (!nrow(cf)) stop("No predictor was selected.", call. = FALSE)
      cf <- cf[order(-abs(cf$coef)), ]
      cf <- utils::head(cf, 40)
      cf$indicator <- factor(cf$indicator, levels = rev(cf$indicator))
      cf$component <- factor(cf$component, levels = .ccva_components)
      ggplot2::ggplot(cf, ggplot2::aes(coef, indicator, colour = component)) +
        ggplot2::geom_vline(xintercept = 0, colour = "#9a9a94") +
        ggplot2::geom_segment(ggplot2::aes(x = 0, xend = coef,
                                           yend = indicator), linewidth = 1) +
        ggplot2::geom_point(size = 2.6) +
        .comp_scale() +
        ggplot2::labs(title = "Predictors selected as drivers of vulnerability",
                      subtitle = sprintf("%d of %d predictors kept; signs constrained by component (alpha = %g)",
                                         fit$details$n_selected,
                                         nrow(fit$details$coefficients),
                                         fit$details$alpha),
                      x = "Standardised coefficient on latent vulnerability",
                      y = NULL) +
        theme_ccva()
    },
    weighted = .plot_contributions(fit, 10),
    stop("No effects plot for method '", fit$method, "'.", call. = FALSE))
}

.plot_map <- function(fit) {
  co <- fit$data$coords
  if (is.null(co)) stop("No coordinates: set `coords` in as_ccva().",
                        call. = FALSE)
  df <- data.frame(x = fit$data$data[[co[1]]], y = fit$data$data[[co[2]]],
                   score = fit$score)
  ggplot2::ggplot(df, ggplot2::aes(x, y, colour = score)) +
    ggplot2::geom_point(size = 2.6, alpha = 0.9) +
    ggplot2::scale_colour_gradient(low = "#d9cbe6", high = "#4b2a6b",
                                   limits = c(0, 1), name = "Vulnerability") +
    ggplot2::coord_quickmap() +
    ggplot2::labs(title = "Where the vulnerable units are",
                  x = co[1], y = co[2]) +
    theme_ccva()
}

## ccva_comparison --------------------------------------------------------------------

#' @export
ccva_plot.ccva_comparison <- function(object, type = NULL, n = 15, ...) {
  type <- match.arg(type %||% "truth", c("truth", "ranks", "validation"))
  fits <- object$fits
  if (type == "truth") {
    tv <- object$truth
    if (is.null(tv)) stop("No truth available; use type = \"validation\".",
                          call. = FALSE)
    tv01 <- .rescale01(tv)
    long <- do.call(rbind, lapply(names(fits), function(nm)
      data.frame(method = nm, truth = tv01, score = fits[[nm]]$score)))
    long$method <- factor(long$method, levels = names(fits))
    rho <- vapply(fits, function(f) stats::cor(f$score, tv,
                                               method = "spearman"),
                  numeric(1))
    lab <- data.frame(method = factor(names(fits), levels = names(fits)),
                      label = sprintf("rho = %.3f", rho))
    return(
      ggplot2::ggplot(long, ggplot2::aes(truth, score)) +
        ggplot2::geom_abline(linetype = 2, colour = "#9a9a94") +
        ggplot2::geom_point(colour = ccva_colours[["vulnerability"]],
                            alpha = 0.45, size = 1.4) +
        ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                             colour = "#33332f", linewidth = 0.7) +
        ggplot2::geom_text(data = lab, ggplot2::aes(x = 0.02, y = 0.98,
                                                    label = label),
                           hjust = 0, vjust = 1, size = 3.6,
                           colour = "#1f1f1c") +
        ggplot2::facet_wrap(~method) +
        ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::labs(title = "Which method recovers the true vulnerability?",
                      subtitle = "Each point is a unit; dashed line = perfect recovery; rho = Spearman rank correlation",
                      x = "True vulnerability (rescaled)",
                      y = "Estimated vulnerability") +
        theme_ccva())
  }
  if (type == "ranks") {
    if (length(fits) < 2) stop("Need two fits.", call. = FALSE)
    a <- fits[[1]]$score; b <- fits[[2]]$score
    df <- data.frame(unit = names(a), rank_a = rank(-a), rank_b = rank(-b))
    df$shift <- abs(df$rank_a - df$rank_b)
    top <- utils::head(df[order(-df$shift), ], 6)
    return(
      ggplot2::ggplot(df, ggplot2::aes(rank_a, rank_b)) +
        ggplot2::geom_abline(linetype = 2, colour = "#9a9a94") +
        ggplot2::geom_point(colour = ccva_colours[["vulnerability"]],
                            alpha = 0.55, size = 1.6) +
        ggplot2::geom_point(data = top, shape = 21, size = 3,
                            colour = "#1f1f1c", fill = NA) +
        ggplot2::geom_text(data = top, ggplot2::aes(label = unit),
                           size = 3, vjust = -1, colour = "#33332f") +
        ggplot2::scale_x_reverse() + ggplot2::scale_y_reverse() +
        ggplot2::coord_equal() +
        ggplot2::labs(title = "How rankings change between methods",
                      subtitle = sprintf("Rank 1 = most vulnerable; circled = largest shifts; rank correlation %.3f",
                                         stats::cor(a, b, method = "spearman")),
                      x = paste("Rank under", names(fits)[1]),
                      y = paste("Rank under", names(fits)[2])) +
        theme_ccva())
  }
  tb <- object$table
  tcols <- setdiff(names(tb), c("method", "class_agreement"))
  long <- do.call(rbind, lapply(tcols, function(tc)
    data.frame(method = tb$method, target = tc, rho = tb[[tc]])))
  long$method <- factor(long$method, levels = rev(tb$method))
  ggplot2::ggplot(long, ggplot2::aes(rho, method)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#9a9a94") +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = rho, yend = method),
                          colour = ccva_colours[["vulnerability"]],
                          linewidth = 1) +
    ggplot2::geom_point(colour = ccva_colours[["vulnerability"]], size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rho)),
                       nudge_y = 0.32, size = 3.1, colour = "#33332f") +
    ggplot2::facet_wrap(~target, nrow = 1) +
    ggplot2::scale_x_continuous(limits = c(min(0, long$rho), 1)) +
    ggplot2::labs(title = "Validation against observed outcomes",
                  subtitle = "Spearman correlation of each method's scores with each target (oriented: higher = better)",
                  x = "Spearman rho", y = NULL) +
    theme_ccva()
}

## ccva_uncertainty -------------------------------------------------------------------

#' @export
ccva_plot.ccva_uncertainty <- function(object, type = NULL, n = 40, ...) {
  type <- match.arg(type %||% "intervals", c("intervals", "pairwise"))
  s <- object$summary
  if (type == "intervals") {
    sel <- if (nrow(s) > n) utils::head(s, n) else s
    sel$unit <- factor(sel$unit, levels = rev(sel$unit))
    return(
      ggplot2::ggplot(sel, ggplot2::aes(estimate, unit)) +
        ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper),
                                height = 0, linewidth = 1.6,
                                colour = "#cbbad9") +
        ggplot2::geom_point(colour = ccva_colours[["vulnerability"]],
                            size = 2.4) +
        ggplot2::labs(title = "How certain is each unit's vulnerability?",
                      subtitle = sprintf("%d%% %s intervals%s; %.0f%% of pairs overlap, %.0f%% of pairwise rankings robust",
                                         round(100 * object$level),
                                         if (object$type == "bootstrap") "bootstrap" else "weight",
                                         if (nrow(s) > n) sprintf(", top %d units", n) else "",
                                         100 * object$overlap_share,
                                         100 * object$robust_share),
                      x = "Vulnerability", y = NULL,
                      caption = if (object$type == "bootstrap")
                        "Each replicate is on its own refit's 0-1 scale; values beyond 1 lie outside that resample's range." else NULL) +
        theme_ccva())
  }
  ids <- utils::head(s$unit, min(n, 30))
  P <- object$pairwise[ids, ids]
  df <- data.frame(var1 = factor(rep(ids, length(ids)), levels = ids),
                   var2 = factor(rep(ids, each = length(ids)),
                                 levels = rev(ids)),
                   p = as.numeric(t(P)))
  df$p[df$var1 == df$var2] <- NA
  ggplot2::ggplot(df, ggplot2::aes(var1, var2, fill = p)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::scale_fill_gradient2(low = "#2a6fb0", mid = "#f4f4f1",
                                  high = "#8a5aa8", midpoint = 0.5,
                                  limits = c(0, 1), na.value = "white",
                                  name = "P(row > column)") +
    ggplot2::labs(title = "Is one unit really more vulnerable than another?",
                  subtitle = "Probability the row unit outranks the column unit; pale = ranking uncertain",
                  x = NULL, y = NULL) +
    ggplot2::coord_equal() + theme_ccva() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,
                                                       vjust = 0.5),
                   panel.grid = ggplot2::element_blank())
}

## ccva_calibration ---------------------------------------------------------------

#' @export
plot.ccva_calibration <- function(x, ...) ccva_plot(x, ...)
#' @export
plot.ccva_sites <- function(x, ...) ccva_plot(x, ...)

#' @export
ccva_plot.ccva_calibration <- function(object, type = NULL, n = 25, ...) {
  type <- match.arg(type %||% "curve", c("curve", "risk", "reliability"))
  lab <- sprintf("P(%s >= %s)", object$outcome, format(object$threshold))
  if (type == "curve") {
    u <- object$units
    u$y <- ifelse(u$event, 1, 0)
    return(
      ggplot2::ggplot(object$curve, ggplot2::aes(score, probability)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                             fill = ccva_colours[["vulnerability"]],
                             alpha = 0.18) +
        ggplot2::geom_line(colour = ccva_colours[["vulnerability"]],
                           linewidth = 1.2) +
        ggplot2::geom_jitter(data = u, ggplot2::aes(score, y), height = 0.025,
                             width = 0, size = 1, alpha = 0.35,
                             colour = "#45453f") +
        ggplot2::scale_y_continuous(limits = c(-0.05, 1.05),
                                    breaks = seq(0, 1, 0.25)) +
        ggplot2::labs(title = "From vulnerability score to risk",
                      subtitle = sprintf("Calibrated %s with %d%% interval; AUC %.2f",
                                         lab, round(100 * object$level),
                                         object$metrics[["auc"]]),
                      x = "Vulnerability score", y = "Probability",
                      caption = paste0("Dots: observed events (top) and non-events (bottom).",
                        if (object$in_sample)
                          " In-sample: this outcome was also used to fit the vulnerability model." else "")) +
        theme_ccva())
  }
  if (type == "reliability") {
    b <- object$bins
    return(
      ggplot2::ggplot(b, ggplot2::aes(predicted, observed)) +
        ggplot2::geom_abline(linetype = 2, colour = "#9a9a94") +
        ggplot2::geom_line(colour = ccva_colours[["vulnerability"]]) +
        ggplot2::geom_point(colour = ccva_colours[["vulnerability"]],
                            size = 3.5) +
        ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::labs(title = "Are the predicted risks reliable?",
                      subtitle = sprintf("Quintiles of predicted risk (about %d units each)",
                                         round(mean(b$n))),
                      caption = "Dashed line = perfect calibration.",
                      x = "Predicted probability", y = "Observed frequency") +
        theme_ccva())
  }
  u <- utils::head(object$units, n)
  u$unit <- factor(u$unit, levels = rev(u$unit))
  ggplot2::ggplot(u, ggplot2::aes(probability, unit)) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper),
                            height = 0, linewidth = 1.6, colour = "#cbbad9") +
    ggplot2::geom_point(colour = ccva_colours[["vulnerability"]], size = 2.4) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(title = sprintf("Risk of decline for the %d highest-risk units", nrow(u)),
                  subtitle = sprintf("%s with %d%% intervals%s", lab,
                                     round(100 * object$level),
                                     if (object$propagated) "; score uncertainty propagated" else ""),
                  x = "Probability", y = NULL) +
    theme_ccva()
}

#' @export
ccva_plot.ccva_sites <- function(object, type = NULL, n = 25, ...) {
  type <- match.arg(type %||% "ranking", c("ranking", "map"))
  d <- as.data.frame(unclass(object))
  if (type == "map") {
    if (!all(c("x", "y") %in% names(d)))
      stop("No site coordinates: pass `coords` to ccva_aggregate().",
           call. = FALSE)
    return(
      ggplot2::ggplot(d, ggplot2::aes(x, y)) +
        ggplot2::geom_point(ggplot2::aes(size = n_species,
                                         colour = mean_vulnerability),
                            alpha = 0.9) +
        ggplot2::scale_colour_gradient(low = "#d9cbe6", high = "#4b2a6b",
                                       name = "Mean\nvulnerability") +
        ggplot2::scale_size_area(max_size = 7, name = "Species") +
        ggplot2::coord_quickmap() +
        ggplot2::labs(title = "Where vulnerable communities are",
                      x = NULL, y = NULL) +
        theme_ccva())
  }
  d <- utils::head(d, n)
  d$site <- factor(d$site, levels = rev(d$site))
  p <- ggplot2::ggplot(d, ggplot2::aes(mean_vulnerability, site))
  if (!is.null(d$lower))
    p <- p + ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper),
                                     height = 0, linewidth = 1.6,
                                     colour = "#cbbad9")
  xr <- max(c(d$mean_vulnerability, d$upper), na.rm = TRUE)
  p + ggplot2::geom_point(colour = ccva_colours[["vulnerability"]], size = 2.8) +
    ggplot2::geom_text(ggplot2::aes(x = xr, label = sprintf("%d spp, %d%% high",
                                                         n_species,
                                                         round(100 * share_high))),
                       hjust = -0.1, size = 3, colour = "#5c5c57") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.35))) +
    ggplot2::labs(title = "Which communities are most vulnerable?",
                  subtitle = "Mean vulnerability of each site's species, with interval when uncertainty is supplied",
                  x = "Mean vulnerability", y = NULL) +
    theme_ccva()
}
