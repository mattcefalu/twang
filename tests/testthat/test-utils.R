context("utils")

test_that("Expect error, missing data NA", {
  
  data <- matrix(c(5, NA, 4, 3))
  expect_error(check_missing(data))
  
})

test_that("Expect error, missing data NaN", {
  
  data <- matrix(c(5, NaN, 4, 3))
  expect_error(check_missing(data))
  
})