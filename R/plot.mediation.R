#' Plot the `mediation` object.
#'
#' @param x weighted_mediation object
#' @param plots An indicator of which type of plot is desired. The options are
#'   * `"optimize" or 1` A plot of the balance criteria as a function of the GBM
#'     iteration.
#'   * `"boxplot" or 2` Boxplots of the propensity scores for the treatment and
#'     control cases
#'   * `"es" or 3` Plots of the standardized effect size of the pre-treatment
#'     variables before and after reweighing
#'   * `"t" or 4` Plots of the p-values from t-statistics comparing means of
#'     treated and control subjects for pretreatment variables, before and after
#'     weighting.
#'   * `"ks" or 5` Plots of the p-values from Kolmogorov-Smirnov statistics
#'     comparing distributions of pretreatment variables of treated and control
#'     subjects, before and after weighting.
#' @param subset Used to restrict which of the `stop.method`s will be used
#'   in the figure. For example `subset = c(1,3)` would indicate that the
#'   first and third `stop.method`s (in alphabetical order of those specified
#'   in the original call to the mediation function) should be included in the
#'   figure.
#' @param color If `color = FALSE`, figures will be gray scale. Default: `TRUE`.
#' @param model_subset integer
#'   Choose either model A or model M only.
#' @param ... Additional arguments.
#'
#' @method plot mediation
#' @export
plot.mediation <- function(x,
                           plots = "optimize",
                           subset = NULL,
                           color = TRUE,
                           model_subset = NULL,
                           ...) {

  # TODO : We'll probably want to fix / clean up this function ...
  if (!is.null(model_subset) && !(model_subset %in% c(1, 2))) {
    stop("The `model_subset` must be NULL, 1, or 2.")
  }

  # TODO : We'll probably want to change this ...
  if (plots == 'density'){

    mediator <- x$data[['M']]
    frames <- list()
    for (method in x$stopping_methods) {

      # we need to scale the weights to sum to 1 for the density plot
      # TODO : We'll want to change these to natural indirect weights, probably ...
      weights <- x[['natural_indirect_weights']][[paste(method, 'ATT', sep='.')]]

      # we replace any NA values with zero
      # TODO : Should we do this?
      weights[is.na(weights)] <- 0

      # we create normalized weights by dividing by the sum of weights
      weights_normed <- weights / sum(weights)

      # we create a new data frame with the mediator, weights, and method
      frames[[method]] <- data.frame(mediator = mediator,
                                     weights = weights,
                                     weights_normed = weights_normed,
                                     method = factor(method))
    }

    new_data <- do.call(rbind, frames)
    rownames(new_data) <- 1:nrow(new_data)

    if (is.factor(mediator)) {
      new_plot <- lattice::barchart(weights ~ mediator | method,
                                    data = new_data,
                                    main = "Weights Distribtion vs Mediator, Bar Chart")

    } else {
      new_plot <- lattice::densityplot(~mediator | method,
                                       data = new_data,
                                       weights = weights_normed,
                                       color = color,
                                       aspect = 1,
                                       main = "Weighted Mediator, Density Plot")
      
    }
    return(new_plot)
  }

  args <- list(plots = plots, subset = subset, color = color)

  model_a <- x$model_a
  model_m <- x$model_m
  model_names <- c('Model A', 'Model M')

  plot1 <- do.call(getS3method("plot", "ps"), c(list(model_a), args))
  plot2 <- do.call(getS3method("plot", "ps"), c(list(model_m), args))

  plot1 <- update(plot1, ylab.right = model_names)
  if (is.null(model_subset)) {
    new_plot <- suppressWarnings(c(plot1, plot2))
    new_plot <- update(new_plot, ylab.right = model_names)
  } else if (model_subset == 1) {
    new_plot <- update(plot1, ylab.right = model_names[[1]])
  } else {
    new_plot <- update(plot2, ylab.right = model_names[[2]])
  }
  return(new_plot)
}