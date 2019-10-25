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
                                                                 estimate_total_effect_wts = TRUE))),
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
                                                                 estimate_total_effect_wts = TRUE))),
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
                                                                 estimate_total_effect_wts = TRUE))),
               "There are NA values in .*",
               ignore.case = TRUE)
  
})


test_that("`weighted_mediation` fails if no weights are provided and `estimate_total_effect_wts =FALSE`", {
  data <- make_fake_data()
  expect_error(suppressWarnings(do.call(weighted_mediation, list(data$treatment,
                                                                 data$mediator,
                                                                 data$x_covariates,
                                                                 y_outcome = data$y_outcome,
                                                                 estimate_total_effect_wts = FALSE))),
               "You have not provided `total_effect_wts` or `total_effect_ps`.*",
               ignore.case = TRUE)
  
})


test_that("`weighted_mediation` fails with NA in m_mediator", {
  data <- make_fake_data()
  expect_error(suppressWarnings(do.call(weighted_mediation, list(data$treatment,
                                                                 data$mediator,
                                                                 data$x_covariates,
                                                                 y_outcome = data$y_outcome,
                                                                 estimate_total_effect_wts = TRUE))),
               "At least two variables are needed in the right-hand side of the formula.*",
               ignore.case = TRUE)
  
})


test_that("`weighted_mediation` fails with only one X variable", {
  data <- read.csv(file.path(getwd() ,'data/test.csv'))
  expect_error(suppressWarnings(do.call(weighted_mediation, list(data[,'A'],
                                                                 data[,'M'],
                                                                 data[, c('X1')],
                                                                 y_outcome = data[,'Yobs'],
                                                                 estimate_total_effect_wts = TRUE))),
               "At least two variables are needed in the right-hand side of the formula.*",
               ignore.case = TRUE)
  
})

test_that("`weighted_mediation` works the same for estimate and passing PS object", {
  data <- read.csv(file.path(getwd() ,'data/test.csv'))
  res1 <- weighted_mediation(data[,'A'],
                             data[,'M'],
                             data[, c('X1', 'X2')],
                             y_outcome = data[,'Yobs'],
                             estimate_total_effect_wts = T,
                             ps_stop.method = c('ks.mean', 'ks.max'),
                             ps_version = 'fast')
  total_effects <- ps(A ~ X1 + X2,
                      data = data,
                      estimand = 'ATE',
                      version = 'fast',
                      booster = 'gbm',
                      n.trees = 10000,
                      interaction.depth = 3,
                      stop.method = c('ks.mean', 'ks.max'),
                      verbose = F)
  res2 <- weighted_mediation(data[,'A'],
                             data[,'M'],
                             data[, c('X1', 'X2')],
                             y_outcome = data[,'Yobs'],
                             total_effect_ps = total_effects,
                             ps_stop.method = c('ks.mean', 'ks.max'),
                             ps_version = 'fast')
  
  expect_equal(res1$natural_direct_wts,
               res2$natural_direct_wts,
               tolerance=1e-1)
  
})

test_that("`weighted_mediation` works the same for estimate and passing weight", {
  data <- read.csv(file.path(getwd() ,'data/test.csv'))
  res1 <- weighted_mediation(data[,'A'],
                             data[,'M'],
                             data[, c('X1', 'X2')],
                             y_outcome = data[,'Yobs'],
                             estimate_total_effect_wts = T,
                             ps_stop.method = c('ks.mean', 'ks.max'),
                             ps_version = 'legacy')

  total_effects <- ps(A ~ X1 + X2,
                      data = data,
                      estimand = 'ATE',
                      version = 'legacy',
                      n.trees = 10000,
                      interaction.depth = 3,
                      shrinkage = 0.01,
                      bag.fraction = 1.0,
                      perm.test.iters = 0,
                      stop.method = c('ks.mean', 'ks.max'),
                      verbose = F)

  res2 <- weighted_mediation(data[,'A'],
                             data[,'M'],
                             data[, c('X1', 'X2')],
                             y_outcome = data[,'Yobs'],
                             total_effect_wts = total_effects$w,
                             ps_stop.method = c('ks.mean', 'ks.max'),
                             ps_version = 'legacy')
  
  expect_equal(res1$ks.mean_effects$total_effect,
               res2$ks.mean_effects$total_effect,
               tolerance=1e-1)
  
  expect_equal(res1$ks.max_effects$total_effect,
               res2$ks.max_effects$total_effect,
               tolerance=1e-1)
  
})


test_that("`weighted_mediation` fails if stop methods are different", {
  data <- read.csv(file.path(getwd() ,'data/test.csv'))
  res1 <- weighted_mediation(data[,'A'],
                             data[,'M'],
                             data[, c('X1', 'X2')],
                             y_outcome = data[,'Yobs'],
                             estimate_total_effect_wts = T,
                             ps_stop.method = c('ks.mean', 'ks.max'),
                             ps_version = 'fast')
  total_effects <- ps(A ~ X1 + X2,
                      data = data,
                      estimand = 'ATE',
                      version = 'legacy',
                      n.trees = 10000,
                      interaction.depth = 3,
                      shrinkage = 0.01,
                      bag.fraction = 1.0,
                      perm.test.iters = 0,
                      stop.method = c('ks.mean'),
                      verbose = F)
  expect_error(weighted_mediation(data[,'A'],
                                  data[,'M'],
                                  data[, c('X1', 'X2')],
                                  y_outcome = data[,'Yobs'],
                                  total_effect_wts = total_effects$w,
                                  ps_stop.method = c('ks.mean', 'ks.max'),
                                  ps_version = 'fast'))
  
})
