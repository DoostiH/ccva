d <- ccva_sim_data()

test_that("constraint objects print and resolve", {
  cc <- ccva_constraints(constrain_moderator(),
                         constrain_equal(c("s1_thermal", "s2_range")),
                         constrain_fixed("s10_genlen", 0),
                         constrain_bounds("exposure", 0, 2))
  expect_s3_class(cc, "ccva_constraints")
  expect_output(print(cc), "moderator")
  r <- ccva:::.resolve_constraints(cc, d)
  expect_equal(unname(r$comp_sign), c(1, 1, -1))
  expect_equal(r$ind_fixed$s10_genlen, 0)
  expect_equal(r$moderator, c(0.5, 0.4))
  expect_true(all(r$ind_sign == 1))
  expect_output(print(ccva_unconstrained()), "free")
})

test_that("predictor signs are stated on the raw scale", {
  r <- ccva:::.resolve_constraints(
    ccva_constraints(constrain_sign("range_change", "-")), d)
  # range_change has direction -1, so "-" on the raw scale is "+" oriented
  expect_equal(unname(r$ind_sign["range_change"]), 1)
})

test_that("invalid constraints are rejected", {
  expect_error(ccva_constraints("x"), "constrain_")
  expect_error(constrain_bounds("exposure", 2, 1), "below")
  expect_error(constrain_equal("a"), "two")
  expect_error(ccva:::.resolve_constraints(
    ccva_constraints(constrain_sign("nope", "+")), d), "not in the data")
  expect_error(ccva:::.resolve_constraints(
    ccva_constraints(constrain_equal(c("E1_risk", "s1_thermal"))), d),
    "one component")
})

test_that("methods can be registered and used", {
  ccva_register("first", function(X, anchor, signs, ...) {
    list(score = X[, 1], predict = function(M) M[, 1])
  }, type = "component", description = "first indicator")
  on.exit(ccva_unregister("first", "component"))
  expect_true("first" %in% ccva_methods()$method)
  s <- ccva_exposure(d, method = "first")
  expect_equal(unname(s$score), d$data$E1_risk)
  expect_error(ccva_register("mean", identity), "built-in")

  ccva_register("max_comp", function(x, components, constraints, ...) {
    list(score = apply(components, 1, max))
  }, type = "combine")
  on.exit(ccva_unregister("max_comp", "combine"), add = TRUE)
  f <- ccva_combine(d, "max_comp")
  expect_s3_class(f, "ccva_fit")
})
