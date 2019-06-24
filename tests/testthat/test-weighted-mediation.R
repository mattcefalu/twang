context("weighted_mediation")

make_fake_data <- function(n = 1000, med = 'binary') {

  if (med == 'binary') {
    mediator <- factor(sample(c(0, 1), n, replace=TRUE))
  } else if (med == 'categorical') {
    mediator <- factor(sample(c(1, 2, 3), n, replace=TRUE))
  } else {
    mediator <- rnorm(n)
  }

  data <- list(treatment = sample(c(0, 1), n, replace=TRUE),
               mediator = mediator,
               x_covariates = rnorm(n, 3),
               x_covariates_mediator = rnorm(n, 1),
               y_outcome = rnorm(n))
  return(data)
}

test_that("`weighted_mediation` fails with more than two treatment categories", {
  data <- make_fake_data()
  data$treatment <- sample(c(0, 1, 2), 1000, replace=TRUE)
  expect_error(suppressWarnings(do.call(weighted_mediation, list(data$treatment,
                                                                 data$mediator,
                                                                 data$x_covariates,
                                                                 y_outcome = data$y_outcome,
                                                                 estimate_ax = TRUE))),
               "The `a_treatment` argument must be .*",
               ignore.case = TRUE)

})

test_that("`weighted_mediation` fails with NaN in a_treatment", {
  data <- make_fake_data()
  data$treatment[3] <- NaN
  expect_error(suppressWarnings(do.call(weighted_mediation, list(data$treatment,
                                                                 data$mediator,
                                                                 data$x_covariates,
                                                                 y_outcome = data$y_outcome,
                                                                 estimate_ax = TRUE))),
               "There are NaN values in .*",
               ignore.case = TRUE)
  
})


test_that("`weighted_mediation` fails with NA in m_mediator", {
  data <- make_fake_data()
  data$mediator[3] <- NA
  expect_error(suppressWarnings(do.call(weighted_mediation, list(data$treatment,
                                                                 data$mediator,
                                                                 data$x_covariates,
                                                                 y_outcome = data$y_outcome,
                                                                 estimate_ax = TRUE))),
               "There are NA values in .*",
               ignore.case = TRUE)
  
})


test_that("`weighted_mediation` fails if no weights are provided and `estimate_ax=FALSE`", {
  data <- make_fake_data()
  expect_error(suppressWarnings(do.call(weighted_mediation, list(data$treatment,
                                                                 data$mediator,
                                                                 data$x_covariates,
                                                                 y_outcome = data$y_outcome,
                                                                 estimate_ax = FALSE))),
               "If no weights are provided via `ax_weights` .*",
               ignore.case = TRUE)
  
})


test_that("`weighted_mediation` fails with NA in m_mediator", {
  data <- make_fake_data()
  expect_error(suppressWarnings(do.call(weighted_mediation, list(data$treatment,
                                                                 data$mediator,
                                                                 data$x_covariates,
                                                                 y_outcome = data$y_outcome,
                                                                 ax_ps = c(1, 2, 3, 5),
                                                                 estimate_ax = TRUE))),
               "The `ax_ps` argument, if provided, must be a `ps` object, .*",
               ignore.case = TRUE)
  
})
