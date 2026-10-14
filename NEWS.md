# ccva 0.1.0

First release.

* Data contract `as_ccva()`: any dataset, columns tagged to exposure,
  sensitivity and adaptive capacity, with directions and outcome proxies.
* Simulator `ccva_simulate()` with a known latent structure (linear or
  nonlinear truth, optional genomic block) and bundled dataset `ccva_sim`.
* Component scoring: Ofori-style ordinal, continuous mean, PCA and
  sign-constrained penalised regression.
* Combination models: weighted index (Ofori et al. 2017, corrected
  normalisation), shape-constrained additive model (scam), constrained
  latent SEM (lavaan) with shared-error method factors, and penalised
  two-stage model (glmnet).
* Constraint system on predictors and components: sign, bounds, equality,
  fixed values and adaptive capacity as moderator.
* Diagnostics, validation, method comparison, bootstrap and weight
  uncertainty (with factor-score error for latent models), calibrated
  classification, plots for every step and an HTML report.
* Method registry `ccva_register()` for contributed methods.
