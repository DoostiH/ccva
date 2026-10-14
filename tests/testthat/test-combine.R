d <- ccva_sim_data()
se <- list(c("E1_risk", "E2_risk", "E3_risk"))

test_that("weighted index equals the corrected Ofori formula", {
  f <- ccva_combine(d, "weighted", component_method = "ordinal")
  C <- f$component_scores
  expect_equal(unname(f$score),
               (C$exposure + C$sensitivity - C$adaptive + 1) / 3)
  expect_true(all(f$score >= 0 & f$score <= 1))
  expect_s3_class(f$class, "factor")
  expect_output(print(f), "weighted index")
  # unequal weights stay on [0, 1]
  fw <- ccva_combine(d, "weighted", weights = c(exposure = 2, adaptive = 0.5))
  expect_true(all(fw$score >= 0 & fw$score <= 1))
})

test_that("weighted index without adaptive capacity is (E + S) / 2", {
  x <- as_ccva(d$data, d$indicators$exposure, d$indicators$sensitivity)
  f <- ccva_combine(x, "weighted")
  C <- f$component_scores
  expect_equal(unname(f$score), (C$exposure + C$sensitivity) / 2)
})

test_that("moderator form has no adaptive main effect", {
  f <- ccva_combine(d, "weighted",
                    constraints = ccva_constraints(constrain_moderator()))
  C <- f$component_scores
  expect_equal(unname(f$score),
               (C$exposure * (1 - 0.5 * C$adaptive) +
                  C$sensitivity * (1 - 0.4 * C$adaptive)) / 2)
  # with E = S = 0 the score is 0 whatever A is
})

test_that("fixing a component at zero removes it", {
  f <- ccva_combine(d, "weighted", constraints = ccva_constraints(
    constrain_fixed("adaptive", 0)))
  C <- f$component_scores
  expect_equal(unname(f$score), (C$exposure + C$sensitivity) / 2)
})

test_that("latent SEM recovers the simulated truth with correct signs", {
  f <- ccva_combine(d, "latent", shared_error = se)
  expect_true(f$details$converged)
  expect_gt(cor(f$score, d$truth$V_star, method = "spearman"), 0.85)
  e <- f$details$effects
  expect_gt(e$std[e$component == "exposure"], 0)
  expect_gt(e$std[e$component == "sensitivity"], 0)
  expect_lt(e$std[e$component == "adaptive"], 0)
  expect_false(f$details$sign_flipped)
  ld <- f$details$loadings
  expect_true(all(ld$estimate > 0))
  expect_true(is.numeric(f$details$score_se))
  ls <- f$details$latent_scores
  expect_gt(cor(ls$sensitivity, d$truth$S_star), 0.9)
  expect_equal(unname(predict(f, d)), unname(f$score), tolerance = 1e-6)
})

test_that("latent constraints bind when they contradict the data", {
  # force adaptive capacity to raise vulnerability: the path must sit at
  # the boundary instead of taking its (negative) free value
  f <- suppressWarnings(ccva_combine(d, "latent", shared_error = se,
    constraints = ccva_constraints(constrain_sign("adaptive", "+"))))
  e <- f$details$effects
  expect_gte(e$estimate[e$component == "adaptive"], -1e-6)
  # bounds
  fb <- ccva_combine(d, "latent", shared_error = se,
    constraints = ccva_constraints(constrain_bounds("exposure", 0, 0.5)))
  eb <- fb$details$effects
  expect_lte(eb$estimate[eb$component == "exposure"], 0.5 + 1e-6)
  # equal loadings
  fe <- ccva_combine(d, "latent", shared_error = se,
    constraints = ccva_constraints(constrain_equal(c("a1_genetic",
                                                     "a3_lowfrag"))))
  ld <- fe$details$loadings
  expect_equal(ld$estimate[ld$indicator == "a1_genetic"],
               ld$estimate[ld$indicator == "a3_lowfrag"], tolerance = 1e-5)
  # moderator is not expressible in the SEM
  expect_error(ccva_combine(d, "latent", constraints = ccva_constraints(
    constrain_moderator())), "moderator")
})

test_that("unconstrained latent fit still sign-aligns its scores", {
  f <- ccva_combine(d, "latent", shared_error = se,
                    constraints = ccva_unconstrained())
  expect_gt(cor(f$score, d$truth$V_star), 0.8)
})

test_that("large blocks enter the SEM as observed component scores", {
  g <- ccva_simulate(200, n_genomic = 40, seed = 3)
  f <- suppressWarnings(ccva_combine(g, "latent", max_indicators = 30))
  expect_true("adaptive" %in% f$details$observed_components)
  expect_true(any(grepl("observed component score", f$notes)))
})

test_that("shape-constrained model is monotone in the constrained direction", {
  skip_if_not_installed("scam")
  f <- ccva_combine(d, "shape")
  cv <- f$details$curves
  for (cmp in c("exposure", "sensitivity")) {
    y <- cv$y[cv$component == cmp]
    expect_true(all(diff(y) >= -1e-8))
  }
  ya <- cv$y[cv$component == "adaptive"]
  expect_true(all(diff(ya) <= 1e-8))
  expect_gt(cor(f$score, d$truth$V_star, method = "spearman"), 0.8)
  fm <- ccva_combine(d, "shape", constraints = ccva_constraints(
    constrain_moderator()))
  expect_true(fm$details$moderator)
  expect_equal(unname(predict(fm, d)), unname(fm$score), tolerance = 1e-6)
})

test_that("penalised two-stage respects component signs", {
  skip_if_not_installed("glmnet")
  set.seed(11)
  f <- ccva_combine(d, "penalised")
  cf <- f$details$coefficients
  expect_true(all(cf$coef[cf$component %in% c("exposure", "sensitivity")] >= 0))
  expect_true(all(cf$coef[cf$component == "adaptive"] <= 0))
  expect_gt(f$details$n_selected, 0)
  expect_gt(cor(f$score, d$truth$V_star, method = "spearman"), 0.7)
  expect_equal(unname(predict(f, d)), unname(f$score), tolerance = 1e-6)
})

test_that("ccva() wrapper accepts a data frame", {
  f <- ccva(ccva_sim, method = "weighted",
            exposure = c("E1_risk", "E2_risk"), sensitivity = "s1_thermal",
            adaptive = "a1_genetic", id = "species")
  expect_s3_class(f, "ccva_fit")
  expect_equal(nrow(ccva_rank(f)), 300)
  expect_true(all(diff(ccva_rank(f)$score) <= 0))
})
