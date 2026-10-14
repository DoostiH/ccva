d <- ccva_sim_data()
se <- list(c("E1_risk", "E2_risk", "E3_risk"))

test_that("directed links recover a simulated A -> E dependency", {
  s <- ccva_simulate(400, link_AE = -0.6, phi = diag(3), seed = 4)
  f <- ccva_combine(s, "latent", shared_error = se,
    constraints = ccva_constraints(constrain_link("adaptive", "exposure")))
  lk <- f$details$links
  expect_equal(nrow(lk), 1)
  expect_lt(lk$std, -0.4)
  te <- f$details$total_effects
  expect_equal(te$total, te$direct + te$indirect)
  expect_lt(te$indirect, 0)
  expect_output(print(f), "total effect of adaptive")
  expect_s3_class(ccva_plot(f, "structure"), "ggplot")
})

test_that("two links from adaptive capacity follow the causal diagram", {
  f <- ccva_combine(d, "latent", shared_error = se,
    constraints = ccva_constraints(
      constrain_link("adaptive", "exposure", "-"),
      constrain_link("adaptive", "sensitivity", "-")))
  expect_equal(nrow(f$details$links), 2)
  expect_true(all(f$details$links$estimate <= 1e-6))
  te <- f$details$total_effects
  expect_equal(te$indirect, sum(f$details$links$indirect))
})

test_that("links are resolved in two stages outside the SEM", {
  s <- ccva_simulate(300, link_AE = -0.6, seed = 4)
  w <- ccva_combine(s, "weighted", constraints = ccva_constraints(
    constrain_link("adaptive", "exposure")))
  expect_equal(w$details$links$method, "two-stage")
  expect_lt(abs(cor(w$component_scores$exposure,
                    w$component_scores$adaptive)), 1e-8)
  expect_equal(unname(predict(w, s)), unname(w$score), tolerance = 1e-8)
  expect_true(any(grepl("two-stage", w$notes)))
})

test_that("invalid links are rejected", {
  expect_error(constrain_link("adaptive", "adaptive"), "itself")
  expect_error(constrain_link("adaptive", "vulnerability"), "must each")
  expect_error(ccva_combine(d, "latent", constraints = ccva_constraints(
    constrain_link("adaptive", "exposure"),
    constrain_link("exposure", "adaptive"))), "cycle")
})

test_that("formative components are not judged by alpha and enter as composites", {
  x <- d; x$formative <- "adaptive"
  chk <- ccva_check(x)
  expect_match(chk$alpha$verdict[chk$alpha$component == "adaptive"],
               "formative")
  expect_true("mean_r" %in% names(chk$alpha))
  f <- ccva_combine(x, "latent", shared_error = se)
  expect_true("adaptive" %in% f$details$observed_components)
  xx <- as_ccva(ccva_sim, c("E1_risk", "E2_risk"), c("s1_thermal", "s2_range"),
                c("a1_genetic", "a2_bodysize"), formative = "adaptive")
  expect_equal(xx$formative, "adaptive")
  expect_output(print(xx), "formative")
})

test_that("curvature shapes map to concave / convex monotone splines", {
  skip_if_not_installed("scam")
  f <- ccva_combine(d, "shape", constraints = ccva_constraints(
    constrain_shape("exposure", "concave"),
    constrain_shape("sensitivity", "convex")))
  fml <- paste(deparse(f$details$formula), collapse = "")
  expect_match(fml, "micv")
  expect_match(fml, "micx")
  cv <- f$details$curves
  y <- cv$y[cv$component == "exposure"]
  expect_true(all(diff(diff(y)) <= 1e-6))
  w <- ccva_combine(d, "weighted", constraints = ccva_constraints(
    constrain_shape("exposure", "concave")))
  expect_true(any(grepl("constrain_shape", w$notes)))
  expect_error(constrain_shape("s1_thermal"), "component")
})

test_that("calibration gives probabilities with intervals", {
  f <- ccva_combine(d, "weighted", component_method = "ordinal")
  cal <- ccva_calibrate(f, "iucn", R = 100, seed = 1)
  expect_s3_class(cal, "ccva_calibration")
  u <- cal$units
  expect_true(all(u$probability >= 0 & u$probability <= 1))
  expect_true(all(u$lower <= u$upper))
  expect_equal(cal$threshold, 2)
  expect_gt(cal$metrics[["auc"]], 0.7)
  expect_gt(cal$metrics[["slope"]], 0)
  expect_false(cal$in_sample)
  fl <- ccva_combine(d, "latent", shared_error = se)
  un <- ccva_uncertainty(fl, R = 12, seed = 1, progress = FALSE)
  cal2 <- ccva_calibrate(fl, uncertainty = un, R = 50, seed = 1)
  expect_true(cal2$in_sample && cal2$propagated)
  for (t in c("curve", "risk", "reliability"))
    expect_s3_class(ccva_plot(cal2, t), "ggplot")
  expect_error(ccva_calibrate(f, "iucn", threshold = 99), "one class")
  cal3 <- ccva_calibrate(f, event = d$data$iucn >= 3, R = 50, seed = 1)
  expect_equal(cal3$outcome, "event")
})

test_that("method recommendations follow sample size and dimensionality", {
  expect_equal(ccva_recommend(d)$tier, "global")
  expect_equal(ccva_recommend(ccva_simulate(40, seed = 1))$tier, "small")
  r <- ccva_recommend(ccva_simulate(100, n_genomic = 100, seed = 1))
  expect_equal(r$tier, "regional")
  expect_true(any(grepl("penalised", r$steps$call)))
  expect_output(print(r), "regional")
  x <- as_ccva(ccva_sim, "E1_risk", "s1_thermal")
  expect_true("outcomes" %in% ccva_recommend(x)$steps$stage)
})

test_that("assay adjustment removes protocol effects", {
  set.seed(1)
  m <- data.frame(species = rep(paste0("sp", 1:30), each = 4),
                  ramp = rep(c(0.1, 0.25, 0.5, 1), 30),
                  assay = rep(c("dynamic", "dynamic", "dynamic", "static"), 30))
  truth <- rnorm(30, 42, 2)
  m$ctmax <- rep(truth, each = 4) + 3 * m$ramp +
    ifelse(m$assay == "static", -1.5, 0) + rnorm(120, 0, 0.2)
  a <- ccva_adjust_method(m, "ctmax", c("ramp", "assay"), unit = "species")
  e <- a$effects
  expect_equal(e$estimate[e$term == "ramp"], 3, tolerance = 0.1)
  expect_equal(e$estimate[e$term == "assaystatic"], -1.5, tolerance = 0.1)
  # adjusted values vary only by species (plus noise)
  within_sd <- tapply(a$data$ctmax_adj, a$data$species, stats::sd)
  expect_lt(mean(within_sd), 0.3)
  expect_gt(cor(a$by_unit$ctmax_adj, truth[order(paste0("sp", 1:30))]), 0.99)
  expect_output(print(a), "reference protocol")
  # confounded design is flagged
  m2 <- data.frame(sp = paste0("s", 1:10), ramp = rep(c(0.1, 1), 5),
                   v = rnorm(10))
  m2 <- rbind(m2, m2)
  expect_warning(ccva_adjust_method(m2, "v", "ramp", unit = "sp"),
                 "not identifiable")
})

test_that("site aggregation summarises communities", {
  f <- ccva_combine(d, "weighted")
  set.seed(3)
  occ <- unique(data.frame(site = rep(paste0("s", 1:6), each = 20),
                           species = sample(names(f$score), 120, TRUE)))
  st <- ccva_aggregate(f, occ)
  expect_s3_class(st, "ccva_sites")
  expect_true(all(diff(st$mean_vulnerability) <= 0))
  expect_false(anyNA(st$share_high))
  one <- occ$species[occ$site == st$site[1]]
  expect_equal(st$mean_vulnerability[1], mean(f$score[one]))
  M <- table(occ$site, occ$species) > 0
  st2 <- ccva_aggregate(f, unclass(M) * 1)
  expect_equal(sort(st2$mean_vulnerability), sort(st$mean_vulnerability))
  u <- ccva_uncertainty(f, "weights", R = 50, seed = 1)
  crd <- data.frame(site = paste0("s", 1:6), x = 1:6, y = 6:1)
  st3 <- ccva_aggregate(f, occ, uncertainty = u, coords = crd)
  expect_true(all(c("lower", "upper", "x") %in% names(st3)))
  expect_s3_class(ccva_plot(st3), "ggplot")
  expect_s3_class(ccva_plot(st3, "map"), "ggplot")
  expect_warning(ccva_aggregate(f, rbind(occ, data.frame(site = "s1",
                                                         species = "zzz"))),
                 "ignored")
})
