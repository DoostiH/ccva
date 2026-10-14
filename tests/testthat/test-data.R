test_that("as_ccva builds a valid object", {
  d <- ccva_sim_data()
  expect_s3_class(d, "ccva_data")
  expect_length(d$id, 300)
  expect_length(d$indicators$sensitivity, 11)
  expect_equal(unname(d$direction["range_change"]), -1)
  expect_equal(unname(d$direction["E1_risk"]), 1)
  expect_output(print(d), "300 units")
  expect_s3_class(summary(d), "data.frame")
})

test_that("as_ccva validates its input", {
  df <- data.frame(id = letters[1:5], e = 1:5, s = c(2, 1, 4, 3, 5),
                   a = 5:1, y = c(1, 2, 2, 3, 1), txt = "x")
  expect_error(as_ccva(df, "nope", "s"), "not found")
  expect_error(as_ccva(df, "e", "txt"), "numeric")
  expect_error(as_ccva(df, "e", "s", outcomes = "e"), "indicator and outcome")
  expect_warning(as_ccva(df, "e", c("s", "e")), "double-counts")
  expect_error(as_ccva(df, "e", "s", direction = c(e = 2)), "direction")
  expect_error(as_ccva(df, character(0), "s"), "at least one")
  df$id[2] <- "a"
  expect_error(as_ccva(df, "e", "s", id = "id"), "unique")
})

test_that("missing values are dropped or rejected", {
  df <- data.frame(e = c(1, NA, 3, 4, 5), s = 1:5, a = 5:1)
  expect_message(x <- as_ccva(df, "e", "s", "a"), "Dropping 1")
  expect_length(x$id, 4)
  expect_error(as_ccva(df, "e", "s", "a", na_action = "fail"), "missing")
})

test_that("simulator produces truth and optional genomic block", {
  s <- ccva_simulate(80, n_genomic = 50, seed = 2)
  expect_s3_class(s, "ccva_data")
  expect_equal(nrow(s$truth), 80)
  expect_length(s$indicators$adaptive, 54)
  expect_true(all(s$data$g0001 %in% 0:2))
  sn <- ccva_simulate(60, nonlinear = TRUE, seed = 2)
  expect_equal(ncol(sn$truth), 4)
})

test_that("predictor helpers compute Ofori exposure and thermal margins", {
  e <- ccva_sdm_exposure(c(100, 50), c(60, 75), c(40, 50), c(20, 10))
  expect_equal(e$E1_change, c(-0.4, 0.5))
  expect_equal(e$E1_risk, c(0.4, -0.5))
  expect_equal(e$E2_risk, c(0.6, 0))
  expect_equal(e$E3_risk, c(0.8, 0.8))
  expect_error(ccva_sdm_exposure(0, 1, 0), "positive")
  expect_error(ccva_sdm_exposure(10, 5, 8), "cannot exceed")
  m <- ccva_thermal_margin(42, 30, 35)
  expect_equal(m$warming_tolerance, 12)
  expect_equal(m$safety_margin, 5)
})
