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

  # TODO : We need a separate impelmentation for outcome mediation
  if (x$mediation_type != 'weighted') {
    stop('Plotting only implemented for weighted mediation objects.')
  }

  # TODO : We'll probably want to fix / clean up this function ...
  if (!is.null(model_subset) && !(model_subset %in% c(1, 2))) {
    stop("The `model_subset` must be NULL, 1, or 2.")
  }

  # TODO : We'll probably want to update this at some point
  if (plots == 'density'){

    mediator <- x$data[['M']]
    treatment <- x$data[['A']]

    # check if the mediator is a factor, because we'll be calculating
    # things differently if the mediator is a factor vs numeric
    mediator_as_factor <- is.factor(mediator)

    frames <- list()
    for (method in x$stopping_methods) {
  
      # extract the weights from the mediation object
      weights <- x[['natural_direct_weights']][, paste(method, 'ATT', sep='.')]
      
      # create normalized weights by dividing by the sum of weights;
      # we do this separately for treatment and control
      weights <- weights / sum(weights[which(treatment == 0)])
      mask <- which(treatment == 1)
      weights[mask] <- (weights[mask] / sum(weights[mask]))

      # if the mediator is a factor
      if (mediator_as_factor) {
  
        # calculate aggregates for the weighted mediator, treatment
        a1 <- aggregate(weights[which(treatment == 1)],
                        by = list(M = factor(mediator[which(treatment == 1)])),
                        FUN = sum)
        a1['A'] <- 1
  
        # calculate aggregates for the weighted mediator, control
        a0 <- aggregate(weights[which(treatment == 0)],
                        by = list(M = factor(mediator[which(treatment == 0)])),
                        FUN = sum)
        a0['A'] <- 0

        # combine the aggregates
        a1a0 <- rbind(a1, a0)
        a1a0['method'] <- method

        # add the combined data frame into the master frames list
        frames[[method]] <- a1a0
      } else {
  
        # create a new data frame an add it to the master frames list
        frames[[method]] <- data.frame(M = mediator,
                                       A = treatment,
                                       weights = weights,
                                       method = method)
        
      }

    }

    new_data <- do.call(rbind, frames)
    rownames(new_data) <- 1:nrow(new_data)

    if (mediator_as_factor) {
      new_plot <- lattice::barchart(x ~ factor(M) | factor(method),
                                    groups = factor(A, levels = c(0, 1), labels = c('control', 'treatment')),
                                    data = new_data,                                   
                                    origin = 0, 
                                    auto.key = TRUE,
                                    ylab = "Proportion",
                                    xlab = "Weighted Mediator",
                                    main = "Weighted Mediation Histogram, Control vs. Treatment",
                                    horiz = FALSE)
      

    } else {
      new_plot <- lattice::densityplot(~M | factor(method),
                                       groups = factor(A, levels = c(0, 1), labels = c('control', 'treatment')),
                                       data = new_data,
                                       weights = weights,
                                       plot.points = FALSE,
                                       origin = 0, 
                                       auto.key = TRUE,
                                       ylab = "Density",
                                       xlab = "Weighted Mediator",
                                       main = "Weighted Mediation Density, Control vs. Treatment")
      
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