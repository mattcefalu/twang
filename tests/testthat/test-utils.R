context("utils")

test_that("Testing weight calculation", {
  ps <- matrix(c(5, 6, 4, 3))
  indicator <- c(1, 1, 0, 0)
  expected <- matrix(c(1/5, 1/6, 1 / (1 - 4), 1 / (1 - 3)))
  wts <- calculate_weights(ps, indicator)
  expect_equal(wts, expected)
})

test_that("Testing weight calculation, opposite", {
  ps <- matrix(c(5, 6, 4, 3))
  indicator <- c(1, 1, 0, 0)
  expected <- matrix(c(1/(1 - 5), 1/ (1 - 6), 1/4, 1/3))
  wts <- calculate_weights(ps, indicator, use_opposite = T)
  expect_equal(wts, expected)
})

test_that("Expect error, missing data NA", {
  
  data <- matrix(c(5, NA, 4, 3))
  expect_error(check_missing(data))
  
})

test_that("Expect error, missing data NaN", {
  
  data <- matrix(c(5, NaN, 4, 3))
  expect_error(check_missing(data))
  
})