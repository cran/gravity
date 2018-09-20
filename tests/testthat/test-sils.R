context("test-gravity.R")

test_that("SILS returns a valid output", {
  # fit model with example dataset
  data("gravity_no_zeros")
  countries_chosen <- names(sort(table(gravity_no_zeros$iso_o), decreasing = TRUE)[1:10])
  grav_small <- gravity_no_zeros[gravity_no_zeros$iso_o %in% countries_chosen, ]

  fit <- sils(
    dependent_variable = "flow", regressors = c("distw", "rta"),
    incomes = c("gdp_o", "gdp_d"),
    maxloop = 100, dec_places = 4, robust = TRUE, verbose = FALSE,
    data = grav_small
  )

  expect_is(fit, "summary.lm")
  expect_is(fit$coefficients, "matrix")
  expect_output(str(fit), "List of 11")
})
