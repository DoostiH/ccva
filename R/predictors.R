#' Compute exposure predictors from species distribution model output
#'
#' Turns the cell counts of thresholded current and future SDM habitat maps
#' into the three exposure indices of Ofori et al. (2017), and recodes them
#' so that a higher value always means more exposure (risk direction), ready
#' to be used as exposure indicators in [as_ccva()].
#'
#' \deqn{E_1 = (N_{future} - N_{current}) / N_{current}}{E1 = (N_future - N_current) / N_current}
#' \deqn{E_2 = N_{overlap} / N_{current}}{E2 = N_overlap / N_current}
#' \deqn{E_3 = N_{protected} / N_{current}}{E3 = N_protected / N_current}
#'
#' with risk recoding \eqn{E_{1,risk} = -E_1}, \eqn{E_{2,risk} = 1 - E_2},
#' \eqn{E_{3,risk} = 1 - E_3}. The three share the denominator
#' \eqn{N_{current}}, so their residuals are correlated; declare them in
#' `shared_error` of [ccva_combine()] when using the latent model.
#'
#' The SDM itself (MaxEnt, BIOCLIM, an ensemble, ...) is fitted with the
#' user's preferred package; this function only needs the cell counts.
#'
#' @param n_current Number of suitable cells under current climate.
#' @param n_future Number of suitable cells under future climate.
#' @param n_overlap Number of cells suitable under both.
#' @param n_protected Optional number of currently suitable cells inside
#'   protected areas.
#'
#' @return A data frame with the raw indices (`E1_change`, `E2_overlap`,
#'   `E3_protected`) and their risk-coded versions (`E1_risk`, `E2_risk`,
#'   `E3_risk`).
#' @export
#' @examples
#' ccva_sdm_exposure(n_current = c(120, 80), n_future = c(60, 95),
#'                   n_overlap = c(50, 70), n_protected = c(30, 10))
ccva_sdm_exposure <- function(n_current, n_future, n_overlap,
                              n_protected = NULL) {
  if (any(n_current <= 0))
    stop("`n_current` must be positive (the species needs current habitat).",
         call. = FALSE)
  if (any(n_overlap > n_current) || any(n_overlap > n_future))
    stop("`n_overlap` cannot exceed `n_current` or `n_future`.",
         call. = FALSE)
  out <- data.frame(E1_change = (n_future - n_current) / n_current,
                    E2_overlap = n_overlap / n_current)
  out$E1_risk <- -out$E1_change
  out$E2_risk <- 1 - out$E2_overlap
  if (!is.null(n_protected)) {
    if (any(n_protected > n_current))
      stop("`n_protected` cannot exceed `n_current`.", call. = FALSE)
    out$E3_protected <- n_protected / n_current
    out$E3_risk <- 1 - out$E3_protected
  }
  out[c(grep("_risk$", names(out), invert = TRUE, value = TRUE),
        grep("_risk$", names(out), value = TRUE))]
}

#' Compute thermal sensitivity predictors
#'
#' Mechanistic sensitivity predictors from measured thermal limits:
#' warming tolerance \eqn{WT = CT_{max} - T_{hab}} and thermal safety margin
#' \eqn{TSM = T_{opt} - T_{hab}}. Smaller (or negative) values mean the
#' organism is closer to its limit, i.e. more sensitive, so give these
#' columns direction `-1` in [as_ccva()].
#'
#' @param ctmax Critical thermal maximum.
#' @param t_habitat Habitat (or operative body) temperature experienced by
#'   the unit, e.g. extracted from a climate layer at its coordinates.
#' @param topt Optional thermal optimum.
#' @return A data frame with `warming_tolerance` and, if `topt` is given,
#'   `safety_margin`.
#' @export
#' @examples
#' ccva_thermal_margin(ctmax = c(42, 45), t_habitat = c(31, 28),
#'                     topt = c(35, 36))
ccva_thermal_margin <- function(ctmax, t_habitat, topt = NULL) {
  out <- data.frame(warming_tolerance = ctmax - t_habitat)
  if (!is.null(topt)) out$safety_margin <- topt - t_habitat
  out
}
