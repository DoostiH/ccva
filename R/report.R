#' Write a reproducible CCVA report
#'
#' Renders an HTML report of a fit: the data and diagnostics, every
#' component, the combination model and its constraints, the ranking, and
#' optionally uncertainty and a method comparison. Needs \pkg{rmarkdown}
#' (and pandoc).
#'
#' @param fit A `ccva_fit`.
#' @param file Output HTML file.
#' @param uncertainty Optional [ccva_uncertainty()] result.
#' @param comparison Optional [ccva_compare()] result.
#' @param title Report title.
#' @param open Open the report in a browser when done.
#' @return The path to the report, invisibly.
#' @export
#' @examples
#' \dontrun{
#' fit <- ccva_combine(ccva_sim_data(), "latent")
#' ccva_report(fit, "my_ccva.html")
#' }
ccva_report <- function(fit, file = "ccva_report.html", uncertainty = NULL,
                        comparison = NULL,
                        title = "Climate change vulnerability assessment",
                        open = interactive()) {
  .need("rmarkdown", "ccva_report()")
  if (!rmarkdown::pandoc_available())
    stop("pandoc is required to render the report.", call. = FALSE)
  tpl <- system.file("report", "ccva_report.Rmd", package = "ccva")
  if (!nzchar(tpl)) stop("Report template not found.", call. = FALSE)
  file <- normalizePath(file, mustWork = FALSE)
  tmp <- tempfile(fileext = ".Rmd")
  file.copy(tpl, tmp)
  env <- new.env()
  env$fit <- fit; env$unc <- uncertainty; env$cmp <- comparison
  env$report_title <- title
  rmarkdown::render(tmp, output_file = basename(file),
                    output_dir = dirname(file), envir = env, quiet = TRUE)
  if (open) utils::browseURL(file)
  invisible(file)
}
