#' Simulate CCVA data with a known latent structure
#'
#' Generates data from a three-layer data-generating process, so that every
#' method can be checked against the truth:
#'
#' * **Layer 1** -- correlated latent exposure \eqn{E^*}, sensitivity
#'   \eqn{S^*} and adaptive capacity \eqn{A^*}.
#' * **Layer 2** -- latent vulnerability
#'   \eqn{V^* = b_E E^* + b_S S^* - b_A A^* + \zeta} (linear), or a nonlinear
#'   version with an exposure plateau, convex sensitivity and adaptive
#'   capacity acting only as a moderator.
#' * **Layer 3** -- observed indicators: three SDM-derived exposure risks
#'   with a shared-denominator nuisance, eleven sensitivity traits, four
#'   adaptive-capacity traits, optional genomic markers, and three
#'   vulnerability proxies (IUCN-like status, population trend, range change).
#'
#' All indicators are continuous except the proxies `iucn` (0--3) and
#' `pop_trend` (-1, 0, 1), which are ordinal by nature. Ofori-style 1--3
#' scores are produced later by `method = "ordinal"` in the component
#' scorers, so the information lost by discretising can be measured.
#'
#' @param n Number of units (species). The full latent model has 21
#'   indicators and four latent factors; use a few hundred for stable fits.
#' @param nonlinear Logical; use the nonlinear structural model for
#'   \eqn{V^*}.
#' @param phi 3 x 3 correlation matrix of \eqn{(E^*, S^*, A^*)}.
#' @param b Structural coefficients `c(E, S, A)`; adaptive capacity enters
#'   with a negative sign.
#' @param sd_zeta Standard deviation of the structural disturbance.
#' @param n_genomic Number of genomic markers (coded 0/1/2) to add to the
#'   adaptive-capacity block, for testing penalised methods when predictors
#'   outnumber units. The first `n_causal` markers carry signal from
#'   \eqn{A^*}; the rest are noise.
#' @param n_causal Number of informative genomic markers.
#' @param seed Optional random seed.
#'
#' @return A `ccva_data` object (see [as_ccva()]) with an extra element
#'   `truth`, a data frame of the latent values `E_star`, `S_star`, `A_star`,
#'   `V_star`. Use `x$data` for the raw data frame.
#' @seealso [ccva_sim] for a bundled simulated dataset.
#' @export
#' @examples
#' sim <- ccva_simulate(n = 200, seed = 1)
#' sim
#' head(sim$truth)
#'
#' # High-dimensional adaptive-capacity block (p > n)
#' sim_g <- ccva_simulate(n = 150, n_genomic = 400, seed = 1)
#' length(sim_g$indicators$adaptive)
ccva_simulate <- function(n = 300, nonlinear = FALSE,
                          phi = matrix(c(1.00, 0.30, -0.20,
                                         0.30, 1.00, -0.25,
                                         -0.20, -0.25, 1.00), 3, 3),
                          b = c(E = 0.50, S = 0.35, A = 0.25),
                          sd_zeta = 0.30, n_genomic = 0, n_causal = 10,
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(n >= 10, all(dim(phi) == 3), length(b) == 3)

  ## Layer 1: correlated latent components
  eta <- matrix(stats::rnorm(n * 3), ncol = 3) %*% chol(phi)
  E <- eta[, 1]; S <- eta[, 2]; A <- eta[, 3]

  ## Layer 2: latent vulnerability
  if (nonlinear) {
    E01 <- .rescale01(E); S01 <- .rescale01(S); A01 <- .rescale01(A)
    V <- (1 - exp(-3 * E01)) * (1 - 0.5 * A01) +
      S01^0.6 * (1 - 0.4 * A01) + stats::rnorm(n, 0, 0.04)
    V <- as.numeric(scale(V))
  } else {
    V <- b[1] * E + b[2] * S - b[3] * A + stats::rnorm(n, 0, sd_zeta)
  }
  V01 <- .rescale01(V)

  ## Layer 3: measurement models
  gen <- function(latent, lambda, resid, common = 0) {
    out <- vapply(seq_along(lambda), function(j)
      lambda[j] * latent + stats::rnorm(n, 0, resid[j]) + common,
      numeric(n))
    colnames(out) <- names(lambda)
    out
  }
  lambda_E <- c(E1_risk = 0.85, E2_risk = 0.75, E3_risk = 0.60)
  lambda_S <- c(s1_thermal = 0.80, s2_range = 0.75, s3_climzones = 0.65,
                s4_habitat = 0.70, s5_diet = 0.60, s6_reprorate = 0.55,
                s7_repromode = 0.50, s8_clutch = 0.55, s9_survival = 0.60,
                s10_genlen = 0.70, s11_lifespan = 0.50)
  lambda_A <- c(a1_genetic = 0.70, a2_bodysize = 0.65, a3_lowfrag = 0.60,
                a4_microhab = 0.75)

  denom <- stats::rnorm(n, 0, 0.35)  # shared N_current nuisance
  X_E <- gen(E, lambda_E, c(0.45, 0.55, 0.70), common = denom)
  X_S <- gen(S, lambda_S, sqrt(1 - lambda_S^2))
  X_A <- gen(A, lambda_A, sqrt(1 - lambda_A^2))

  X_G <- NULL
  if (n_genomic > 0) {
    n_causal <- min(n_causal, n_genomic)
    X_G <- vapply(seq_len(n_genomic), function(j) {
      liab <- (if (j <= n_causal) 0.6 * A else 0) + stats::rnorm(n)
      maf <- stats::runif(1, 0.1, 0.5)
      as.numeric(cut(liab, stats::quantile(liab, c(0, (1 - maf)^2,
                                                    1 - maf^2, 1)),
                     include.lowest = TRUE, labels = FALSE)) - 1
    }, numeric(n))
    colnames(X_G) <- sprintf("g%04d", seq_len(n_genomic))
  }

  ## Vulnerability proxies
  iucn <- cut(pmin(pmax(V01 + stats::rnorm(n, 0, 0.15), 0), 1),
              c(-Inf, 0.25, 0.50, 0.75, Inf), labels = FALSE) - 1
  pop_trend <- ifelse(V01 > 0.65, -1,
                      ifelse(V01 > 0.35,
                             sample(c(-1, 0), n, TRUE, c(0.4, 0.6)), 1))
  range_change <- -0.4 * V01 + stats::rnorm(n, 0, 0.12)

  dat <- data.frame(species = sprintf("Species_%03d", seq_len(n)),
                    X_E, X_S, X_A,
                    iucn = iucn, pop_trend = pop_trend,
                    range_change = range_change)
  if (!is.null(X_G)) dat <- cbind(dat, X_G)

  x <- as_ccva(dat, exposure = names(lambda_E),
               sensitivity = names(lambda_S),
               adaptive = c(names(lambda_A), colnames(X_G)),
               outcomes = c("iucn", "pop_trend", "range_change"),
               direction = c(pop_trend = -1, range_change = -1),
               id = "species")
  x$truth <- data.frame(E_star = E, S_star = S, A_star = A, V_star = V)
  x
}

#' Simulated CCVA dataset
#'
#' 300 simulated species generated by [ccva_simulate()] with `seed = 99`
#' and the linear structural model. Contains the observed columns a modeller
#' would see plus the latent truth for validation. All values are simulated
#' and have no biological meaning.
#'
#' @format A data frame with 300 rows and 26 columns:
#' \describe{
#'   \item{species}{Unit identifier.}
#'   \item{E1_risk}{Exposure: habitat loss risk (SDM-derived, higher = more exposed).}
#'   \item{E2_risk}{Exposure: lack of overlap between current and future habitat.}
#'   \item{E3_risk}{Exposure: lack of protection of current habitat. The three exposure risks share a denominator, hence correlated residuals.}
#'   \item{s1_thermal}{Sensitivity: thermal tolerance (higher = more sensitive).}
#'   \item{s2_range}{Sensitivity: small geographic range.}
#'   \item{s3_climzones}{Sensitivity: few climate zones occupied.}
#'   \item{s4_habitat}{Sensitivity: habitat specialisation.}
#'   \item{s5_diet}{Sensitivity: diet specialisation.}
#'   \item{s6_reprorate}{Sensitivity: low reproductive rate.}
#'   \item{s7_repromode}{Sensitivity: reproductive mode.}
#'   \item{s8_clutch}{Sensitivity: small clutch size.}
#'   \item{s9_survival}{Sensitivity: low offspring survival.}
#'   \item{s10_genlen}{Sensitivity: long generation length.}
#'   \item{s11_lifespan}{Sensitivity: lifespan.}
#'   \item{a1_genetic}{Adaptive capacity: genetic diversity (higher = more capacity).}
#'   \item{a2_bodysize}{Adaptive capacity: body size as a dispersal proxy.}
#'   \item{a3_lowfrag}{Adaptive capacity: low habitat fragmentation.}
#'   \item{a4_microhab}{Adaptive capacity: microhabitat buffering.}
#'   \item{iucn}{IUCN-like status, 0 (least concern) to 3 (endangered).}
#'   \item{pop_trend}{Population trend: -1 declining, 0 stable, 1 increasing.}
#'   \item{range_change}{Proportional range change; negative = loss.}
#'   \item{E_star}{Latent exposure (simulation truth).}
#'   \item{S_star}{Latent sensitivity (simulation truth).}
#'   \item{A_star}{Latent adaptive capacity (simulation truth).}
#'   \item{V_star}{Latent vulnerability (simulation truth).}
#' }
#' @source Simulated with [ccva_simulate()]; see `data-raw/ccva_sim.R`.
#' @examples
#' head(ccva_sim[, 1:6])
"ccva_sim"

#' Specification of the bundled simulated data
#'
#' Convenience wrapper that turns [ccva_sim] into a `ccva_data` object with
#' the correct component tags and directions, and attaches the latent truth.
#'
#' @param data Defaults to [ccva_sim].
#' @return A `ccva_data` object with a `truth` element.
#' @export
#' @examples
#' d <- ccva_sim_data()
#' d
ccva_sim_data <- function(data = ccva::ccva_sim) {
  x <- as_ccva(data,
               exposure = c("E1_risk", "E2_risk", "E3_risk"),
               sensitivity = grep("^s[0-9]+_", names(data), value = TRUE),
               adaptive = grep("^a[0-9]+_", names(data), value = TRUE),
               outcomes = c("iucn", "pop_trend", "range_change"),
               direction = c(pop_trend = -1, range_change = -1),
               id = "species")
  tcols <- c("E_star", "S_star", "A_star", "V_star")
  if (all(tcols %in% names(data))) x$truth <- data[tcols]
  x
}
