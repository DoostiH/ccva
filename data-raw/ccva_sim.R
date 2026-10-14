# Generates the bundled `ccva_sim` dataset.
# Run from the package root after devtools::load_all():
#   source("data-raw/ccva_sim.R")
#
# All values are simulated and have no biological meaning. The design
# follows ccva_latent_data_generation.R (Doosti 2026): correlated latent
# exposure, sensitivity and adaptive capacity drive a latent vulnerability,
# each measured by observed indicators, with three vulnerability proxies.

sim <- ccva_simulate(n = 300, nonlinear = FALSE, seed = 99)
ccva_sim <- cbind(sim$data, sim$truth)
rownames(ccva_sim) <- NULL
usethis_available <- requireNamespace("usethis", quietly = TRUE)
if (usethis_available) {
  usethis::use_data(ccva_sim, overwrite = TRUE)
} else {
  save(ccva_sim, file = "data/ccva_sim.rda", compress = "xz")
}
