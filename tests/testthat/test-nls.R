context("test-gravity.R")

test_that("NLS returns a valid output", {
  # fit model with example dataset
  data("gravity_no_zeros")
  countries_chosen <- names(sort(table(gravity_no_zeros$iso_o), decreasing = TRUE)[1:10])
  grav_small <- gravity_no_zeros[gravity_no_zeros$iso_o %in% countries_chosen, ]

  library(dplyr)

  grav_small <- grav_small %>%
    mutate(
      lgdp_o = log(gdp_o),
      lgdp_d = log(gdp_d)
    )

  fit <- nls(
    dependent_variable = "flow", regressors = c("distw", "rta", "lgdp_o", "lgdp_d"),
    robust = TRUE, data = grav_small
  )

  expect_is(fit, "summary.lm")
  expect_is(fit$coefficients, "matrix")
  expect_output(str(fit), "List of 9")
})
