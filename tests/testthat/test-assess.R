d <- ccva_sim_data()
fw <- ccva_combine(d, "weighted", component_method = "ordinal")
fl <- ccva_combine(d, "latent",
                   shared_error = list(c("E1_risk", "E2_risk", "E3_risk")))

test_that("classification methods give three ordered classes", {
  for (m in c("uniform", "quantile", "mixture", "roc")) {
    cl <- ccva_classify(fl$score, m, data = d)
    expect_s3_class(cl, "ordered")
    expect_equal(levels(cl), c("Low", "Moderate", "High"))
    br <- attr(cl, "cutoffs")
    expect_true(br[1] < br[2])
  }
  q <- ccva_classify(fl$score, "quantile")
  expect_equal(as.numeric(table(q)), c(100, 100, 100))
  m <- ccva_classify(c(0.1, 0.5, 0.9), c(0.3, 0.7))
  expect_equal(as.character(m), c("Low", "Moderate", "High"))
  expect_error(ccva_classify(fl$score, "roc"), "needs `data`")
})

test_that("diagnostics report consistency and correlation", {
  chk <- ccva_check(d)
  expect_s3_class(chk, "ccva_check")
  expect_true(all(chk$alpha$alpha > 0.7))
  expect_equal(dim(chk$component_cor), c(3, 3))
  expect_output(print(chk), "Cronbach")
  g <- ccva_simulate(60, n_genomic = 100, seed = 1)
  cg <- ccva_check(g)
  expect_match(cg$alpha$verdict[cg$alpha$component == "adaptive"],
               "high-dimensional")
})

test_that("validation and comparison against outcomes and truth", {
  v <- ccva_validate(fl)
  expect_true("V* (truth)" %in% v$target)
  expect_true(all(v$spearman > 0))
  cmp <- ccva_compare(ofori = fw, latent = fl)
  expect_s3_class(cmp, "ccva_comparison")
  expect_gt(cmp$table[cmp$table$method == "latent", "V* (truth)"],
            cmp$table[cmp$table$method == "ofori", "V* (truth)"])
  expect_output(print(cmp), "Spearman")
  expect_error(ccva_compare(fw, fl), "name")
})

test_that("weight uncertainty is reproducible and well formed", {
  u1 <- ccva_uncertainty(fw, "weights", R = 200, seed = 4)
  u2 <- ccva_uncertainty(fw, "weights", R = 200, seed = 4)
  expect_equal(u1$summary, u2$summary)
  s <- u1$summary
  expect_true(all(s$lower <= s$upper))
  expect_true(u1$robust_share >= 0 && u1$robust_share <= 1)
  P <- u1$pairwise
  expect_true(all(diag(P) == 0))
  expect_true(all(P + t(P) <= 1 + 1e-12))
  expect_error(ccva_uncertainty(fl, "weights"), "weighted index")
  expect_output(print(u1), "random weights")
})

test_that("bootstrap refits the whole pipeline", {
  small <- ccva_simulate(80, seed = 9)
  f <- ccva_combine(small, "weighted")
  u <- ccva_uncertainty(f, "bootstrap", R = 15, seed = 2, progress = FALSE)
  expect_equal(nrow(u$summary), 80)
  expect_equal(u$R, 15)
  expect_true(all(u$summary$width >= 0))
})

test_that("every plot type renders", {
  chk <- ccva_check(d)
  cmp <- ccva_compare(ofori = fw, latent = fl)
  u <- ccva_uncertainty(fw, "weights", R = 100, seed = 1)
  plots <- list(
    ccva_plot(d, "indicators"), ccva_plot(d, "correlation"),
    ccva_plot(chk, "alpha"), ccva_plot(chk, "correlation"),
    ccva_plot(fl$components$exposure, "importance"),
    ccva_plot(fl$components$exposure, "distribution"),
    ccva_plot(fw, "structure"), ccva_plot(fl, "structure"),
    ccva_plot(fl, "components"), ccva_plot(fw, "scores"),
    ccva_plot(fl, "ranking", n = 5), ccva_plot(fw, "contributions", n = 5),
    ccva_plot(fl, "effects"), ccva_plot(cmp, "truth"),
    ccva_plot(cmp, "ranks"), ccva_plot(cmp, "validation"),
    ccva_plot(u, "intervals"), ccva_plot(u, "pairwise"))
  for (p in plots) {
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot2::ggplot_build(p))
  }
  expect_error(ccva_plot(fl, "map"), "coordinates")
  expect_error(ccva_plot(fl, "contributions"), "weighted")
  expect_error(ccva_plot(1), "No ccva_plot")
})
