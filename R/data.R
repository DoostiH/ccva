#' Build a CCVA dataset
#'
#' Tags the columns of any data frame as indicators of the three CCVA
#' components, gives each indicator a direction, and names the observable
#' outcomes (vulnerability proxies such as IUCN status, population trend or
#' range change). Every other function in the package works on the resulting
#' `ccva_data` object, so this is the package's data contract: nothing about a
#' particular taxon or study is assumed.
#'
#' @param data A data frame with one row per assessed unit (species or
#'   population).
#' @param exposure,sensitivity,adaptive Character vectors naming the indicator
#'   columns of each component. Each component may use one column or many, and
#'   the sets may differ in size. `adaptive` may be `NULL` for an
#'   exposure-sensitivity (ES) assessment.
#' @param outcomes Optional character vector naming observable vulnerability
#'   proxies. Needed by the latent, shape-constrained and penalised models,
#'   and by [ccva_validate()].
#' @param id Optional name of the column holding unit identifiers. If `NULL`,
#'   rows are labelled `unit_1`, `unit_2`, ...
#' @param direction Optional named vector of `+1` / `-1` (or `"+"` / `"-"`)
#'   giving, for each indicator, whether a higher value means *more* of its
#'   component (more exposure, more sensitivity, more adaptive capacity), and
#'   for each outcome whether a higher value means *more vulnerable*. Columns
#'   not listed default to `+1`. For example a thermal tolerance ratio, where
#'   higher means less sensitive, takes `-1`; range change, where a negative
#'   value is range loss, takes `-1`.
#' @param coords Optional character vector of two column names (longitude,
#'   latitude) used by [ccva_plot()] for maps.
#' @param na_action `"omit"` (default) drops rows with missing values in any
#'   used column, with a message; `"fail"` stops instead.
#'
#' @return An object of class `ccva_data`.
#' @seealso [ccva_simulate()] for data with a known latent structure.
#' @export
#' @examples
#' d <- as_ccva(ccva_sim,
#'   exposure    = c("E1_risk", "E2_risk", "E3_risk"),
#'   sensitivity = grep("^s[0-9]+_", names(ccva_sim), value = TRUE),
#'   adaptive    = grep("^a[0-9]+_", names(ccva_sim), value = TRUE),
#'   outcomes    = c("iucn", "pop_trend", "range_change"),
#'   direction   = c(pop_trend = -1, range_change = -1),
#'   id = "species")
#' d
as_ccva <- function(data, exposure, sensitivity, adaptive = NULL,
                    outcomes = NULL, id = NULL, direction = NULL,
                    coords = NULL, na_action = c("omit", "fail")) {
  na_action <- match.arg(na_action)
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  data <- as.data.frame(data)

  ind <- list(exposure = .chr(exposure, "exposure"),
              sensitivity = .chr(sensitivity, "sensitivity"),
              adaptive = if (is.null(adaptive)) character(0) else
                .chr(adaptive, "adaptive"))
  if (length(ind$exposure) == 0L || length(ind$sensitivity) == 0L)
    stop("`exposure` and `sensitivity` need at least one column each.",
         call. = FALSE)
  outcomes <- if (is.null(outcomes)) character(0) else
    .chr(outcomes, "outcomes")

  used <- c(unlist(ind), outcomes)
  missing_cols <- setdiff(c(used, id, coords), names(data))
  if (length(missing_cols))
    stop("Columns not found in `data`: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)

  dup <- unique(unlist(ind)[duplicated(unlist(ind))])
  if (length(dup))
    warning("Indicator(s) used in more than one component: ",
            paste(dup, collapse = ", "),
            ". This double-counts information (e.g. body size in both ",
            "sensitivity and adaptive capacity).", call. = FALSE)
  overlap <- intersect(unlist(ind), outcomes)
  if (length(overlap))
    stop("Column(s) used both as indicator and outcome: ",
         paste(overlap, collapse = ", "), call. = FALSE)

  nonnum <- used[!vapply(data[used], is.numeric, logical(1))]
  if (length(nonnum))
    stop("Indicator and outcome columns must be numeric. Not numeric: ",
         paste(nonnum, collapse = ", "),
         ". Recode ordered categories as integers first.", call. = FALSE)

  if (!is.null(coords) && length(coords) != 2L)
    stop("`coords` must name two columns (longitude, latitude).",
         call. = FALSE)

  # Missing values
  cc <- stats::complete.cases(data[used])
  if (!all(cc)) {
    if (na_action == "fail")
      stop(sum(!cc), " row(s) have missing values in used columns.",
           call. = FALSE)
    message("Dropping ", sum(!cc), " row(s) with missing values.")
    data <- data[cc, , drop = FALSE]
  }
  if (nrow(data) < 3L) stop("Need at least 3 complete rows.", call. = FALSE)

  ids <- if (is.null(id)) paste0("unit_", seq_len(nrow(data))) else
    as.character(data[[id]])
  if (anyDuplicated(ids)) stop("Unit identifiers must be unique.",
                               call. = FALSE)

  dir <- .parse_direction(direction, used)

  structure(list(
    data = data,
    id = ids,
    id_col = id,
    indicators = ind,
    outcomes = outcomes,
    direction = dir,
    coords = coords
  ), class = "ccva_data")
}

#' @export
print.ccva_data <- function(x, ...) {
  cat("<ccva_data> ", length(x$id), " units\n", sep = "")
  for (cmp in .ccva_components) {
    v <- x$indicators[[cmp]]
    if (!length(v)) next
    cat(sprintf("  %-18s %2d indicator%s: %s\n", .ccva_labels[[cmp]],
                length(v), if (length(v) == 1) " " else "s",
                .abbrev(.signed(v, x$direction))))
  }
  if (length(x$outcomes))
    cat(sprintf("  %-18s %2d proxy%s:     %s\n", "Outcomes",
                length(x$outcomes),
                if (length(x$outcomes) == 1) " " else "s",
                .abbrev(.signed(x$outcomes, x$direction))))
  else cat("  Outcomes: none (latent and validation steps unavailable)\n")
  cat("  (-) marks indicators where a higher value means less of the",
      "component.\n")
  invisible(x)
}

#' @export
summary.ccva_data <- function(object, ...) {
  rows <- lapply(c(.ccva_components, "outcome"), function(cmp) {
    v <- if (cmp == "outcome") object$outcomes else object$indicators[[cmp]]
    if (!length(v)) return(NULL)
    d <- object$data[v]
    data.frame(component = cmp, indicator = v,
               direction = unname(object$direction[v]),
               mean = vapply(d, mean, numeric(1)),
               sd = vapply(d, stats::sd, numeric(1)),
               min = vapply(d, min, numeric(1)),
               max = vapply(d, max, numeric(1)),
               row.names = NULL)
  })
  do.call(rbind, rows)
}

## Internal helpers ---------------------------------------------------------

.chr <- function(x, what) {
  if (!is.character(x)) stop("`", what, "` must be a character vector of ",
                             "column names.", call. = FALSE)
  unique(x)
}

.parse_direction <- function(direction, cols) {
  dir <- stats::setNames(rep(1, length(cols)), cols)
  if (is.null(direction)) return(dir)
  if (is.null(names(direction)))
    stop("`direction` must be a named vector.", call. = FALSE)
  val <- vapply(direction, function(d) {
    if (is.character(d)) d <- switch(d, "+" = 1, "-" = -1, NA)
    as.numeric(d)
  }, numeric(1))
  if (any(!val %in% c(-1, 1)))
    stop("`direction` values must be +1 / -1 or \"+\" / \"-\".",
         call. = FALSE)
  unknown <- setdiff(names(val), cols)
  if (length(unknown))
    warning("`direction` given for unused column(s): ",
            paste(unknown, collapse = ", "), call. = FALSE)
  keep <- intersect(names(val), cols)
  dir[keep] <- val[keep]
  dir
}

.signed <- function(v, dir) ifelse(dir[v] < 0, paste0(v, " (-)"), v)

.abbrev <- function(v, n = 4) {
  if (length(v) <= n) return(paste(v, collapse = ", "))
  paste0(paste(v[seq_len(n)], collapse = ", "), ", ... +", length(v) - n)
}

# Indicator matrix for a component, oriented so that higher = more of the
# component (more exposure / sensitivity / adaptive capacity).
.oriented <- function(x, component, data = x$data) {
  v <- x$indicators[[component]]
  m <- as.matrix(data[v])
  sweep(m, 2, x$direction[v], `*`)
}

# Outcome matrix oriented so that higher = more vulnerable.
.oriented_outcomes <- function(x, data = x$data) {
  if (!length(x$outcomes)) return(NULL)
  m <- as.matrix(data[x$outcomes])
  sweep(m, 2, x$direction[x$outcomes], `*`)
}

# Row subset of a ccva_data object (used by the bootstrap).
.subset_ccva <- function(x, rows) {
  x$data <- x$data[rows, , drop = FALSE]
  x$id <- make.unique(x$id[rows])
  if (!is.null(x$truth)) x$truth <- x$truth[rows, , drop = FALSE]
  x
}

.check_ccva <- function(x) {
  if (!inherits(x, "ccva_data"))
    stop("Expected a `ccva_data` object; build one with as_ccva().",
         call. = FALSE)
  invisible(x)
}

.rescale01 <- function(v, lo = min(v), hi = max(v)) {
  if (hi - lo < .Machine$double.eps) return(rep(0.5, length(v)))
  (v - lo) / (hi - lo)
}
