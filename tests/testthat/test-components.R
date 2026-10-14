d <- ccva_sim_data()

test_that("ordinal scoring reproduces the Ofori scheme", {
  e <- ccva_exposure(d, method = "ordinal")
  expect_true(all(e$score >= 1 / 3 - 1e-12 & e$score <= 1 + 1e-12))
  items <- e$details$items
  expect_true(all(items %in% 1:3))
  expect_equal(unname(e$score), unname(rowSums(items) / 9))
  # explicit breaks and a reversed indicator
  df <- data.frame(t = c(0.5, 2, 3.5), r = c(10, 30, 60))
  x <- as_ccva(df, "t", "r", direction = c(t = -1))
  s <- ccva_exposure(x, method = "ordinal", breaks = list(t = c(1.5, 3)))
  # higher thermal ratio = less exposed: 0.5 -> 3, 2 -> 2, 3.5 -> 1
  expect_equal(unname(s$score), c(3, 2, 1) / 3)
})

test_that("Ofori weights act as maximum scores", {
  w <- ccva_sensitivity(d, method = "ordinal",
                        weights = c(s1_thermal = 11, s2_range = 10))
  expect_equal(sum(w$importance), 1)
  expect_gt(w$importance[["s1_thermal"]], w$importance[["s3_climzones"]])
})

test_that("continuous and PCA scores lie in [0, 1] and predict themselves", {
  for (m in c("mean", "pca")) {
    s <- ccva_sensitivity(d, method = m)
    expect_true(all(s$score >= -1e-12 & s$score <= 1 + 1e-12))
    expect_equal(unname(s$predict(d$data)), unname(s$score))
  }
  p <- ccva_sensitivity(d, method = "pca")
  expect_true(all(p$importance > 0))
  expect_gt(cor(p$score, d$truth$S_star), 0.9)
})

test_that("indicator direction is honoured", {
  a <- ccva_adaptive(d, method = "mean")
  expect_gt(cor(a$score, d$truth$A_star), 0.85)
})

test_that("fixed-at-zero indicators are dropped", {
  s <- ccva_sensitivity(d, constraints = ccva_constraints(
    constrain_fixed("s8_clutch", 0)))
  expect_false("s8_clutch" %in% s$indicators)
  expect_equal(s$dropped, "s8_clutch")
})

test_that("penalised scoring selects informative genomic markers with the right sign", {
  skip_if_not_installed("glmnet")
  set.seed(5)
  g <- ccva_simulate(150, n_genomic = 200, n_causal = 10, seed = 5)
  a <- suppressWarnings(ccva_adaptive(g, method = "penalised"))
  b <- a$importance
  expect_true(all(b >= 0))  # oriented: higher marker = more capacity
  sel <- names(b)[b != 0]
  causal <- sprintf("g%04d", 1:10)
  expect_gt(mean(sel %in% c(causal, g$indicators$adaptive[1:4])), 0.4)
  expect_gt(cor(a$score, g$truth$A_star), 0.5)
})

test_that("component scorers validate input", {
  expect_error(ccva_component(d, "exposure", method = "nope"), "Unknown")
  expect_error(ccva_exposure(d, weights = c(zzz = 1)), "Unknown weight")
  expect_error(ccva_component(list(), "exposure"), "ccva_data")
  expect_output(print(ccva_exposure(d)), "Exposure")
})
