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
  
  # return error if ask for plots="optimize" or 1 for method!=ps
  if (x$method!="ps" & (plots=="optimize" || plots==1)) { stop("The optimize plot is only available for method='ps'.")}
  # we want the mediator and the treatment variables
  mediators <- x$data[,x$mediator_names, drop = F]
  treatment <- x$data[[x$a_treatment]]
  
  if (plots != 'density') {
    args <- list(plots = plots, subset = subset, color = color)
    if(x$method=="logistic") {
       x$model_a$ps  <- data.frame(logistic=predict(x$model_a,type="response"))
       x$model_m0$ps  <- data.frame(logistic=predict(x$model_m0,type="response"))
       x$model_a$w  <- data.frame(logistic=ifelse(x$data[,x$a_treatment]==1,1/x$model_a$ps[,1],1/(1-x$model_a$ps[,1])))
       x$model_m0$w  <- data.frame(logistic=ifelse(x$data[,x$a_treatment]==1,1/x$model_m0$ps[,1],1/(1-x$model_m0$ps[,1])))
	 x$model_a$treat <- x$data[,x$a_treatment]
	 x$model_m0$treat <- x$data[,x$a_treatment]

	x$model_a$desc  <- vector("list",length=2)
	names(x$model_a$desc)  <- c("unw","logistic")
	bala	<- bal.table(x)
	x$model_a$desc[[1]]$bal.tab$results <- bala$balance_a[grep("unw",rownames(bala$balance_a)),1:8]
	x$model_a$desc[[2]]$bal.tab$results <- bala$balance_a[grep("logistic",rownames(bala$balance_a)),1:8]

	x$model_m0$desc  <- vector("list",length=2)
	names(x$model_m0$desc)  <- c("unw","logistic")
	x$model_m0$desc[[1]]$bal.tab$results <- bala$balance_m0[grep("unw",rownames(bala$balance_m0)),1:8]
	x$model_m0$desc[[2]]$bal.tab$results <- bala$balance_m0[grep("logistic",rownames(bala$balance_m0)),1:8]
    }
    if(x$method=="crossval") {
       best.iter.a 		<- gbm:::gbm.perf(x$model_a, method="cv",plot.it=FALSE)
       best.iter.m 		<- gbm:::gbm.perf(x$model_m0, method="cv",plot.it=FALSE)
       x$model_a$ps  <- data.frame(crossval=predict(x$model_a, n.trees=best.iter.a,newdata=x$data,type="response"))
       x$model_m0$ps  <- data.frame(crossval=predict(x$model_m0,n.trees=best.iter.m,newdata=x$data,type="response"))
       x$model_a$w  <- data.frame(crossval=ifelse(x$data[,x$a_treatment]==1,1/x$model_a$ps[,1],1/(1-x$model_a$ps[,1])))
       x$model_m0$w  <- data.frame(crossval=ifelse(x$data[,x$a_treatment]==1,1/x$model_m0$ps[,1],1/(1-x$model_m0$ps[,1])))
	 x$model_a$treat <- x$data[,x$a_treatment]
	 x$model_m0$treat <- x$data[,x$a_treatment]

	x$model_a$desc  <- vector("list",length=2)
	names(x$model_a$desc)  <- c("unw","crossval")
	bala	<- bal.table(x)
	x$model_a$desc[[1]]$bal.tab$results <- bala$balance_a[grep("unw",rownames(bala$balance_a)),1:8]
	x$model_a$desc[[2]]$bal.tab$results <- bala$balance_a[grep("crossval",rownames(bala$balance_a)),1:8]

	x$model_m0$desc  <- vector("list",length=2)
	names(x$model_m0$desc)  <- c("unw","crossval")
	x$model_m0$desc[[1]]$bal.tab$results <- bala$balance_m0[grep("unw",rownames(bala$balance_m0)),1:8]
	x$model_m0$desc[[2]]$bal.tab$results <- bala$balance_m0[grep("crossval",rownames(bala$balance_m0)),1:8]
    }
    model_a <- x$model_a
    model_m <- x$model_m0
    model_names <- c('Model A', 'Model M')
    
    plot1 <- do.call(getS3method("plot", "ps"), c(list(model_a), args))
    plot2 <- do.call(getS3method("plot", "ps"), c(list(model_m), args))
    
    plot1 <- update(plot1, ylab.right = model_names)
    if (is.null(model_subset)) {
      new_plot <- suppressWarnings(c(plot1, plot2))
      new_plot <- update(new_plot, ylab.right = rev(model_names),layout=c(length(new_plot$packet.sizes)/2,2))
    } else if (model_subset == 1) {
      new_plot <- update(plot1, ylab.right = model_names[[1]])
    } else {
      new_plot <- update(plot2, ylab.right = model_names[[2]])
    }
    new_plot$as.table	<- TRUE
    return(new_plot)
  }

  # we separate out the mediators that are factors and those that aren't
  mediator_is_factor <- x$mediator_names %in% names(Filter(is.factor, x$data[,x$mediator_names, drop=F]))
  
  mediators_factors <- x$mediator_names[mediator_is_factor == T, drop = T]
  mediators_numbers <- x$mediator_names[mediator_is_factor == F, drop = T]
  
  # get the actual weights from the mediation object
  w_00 <- attr(x, 'w_00')
  w_10 <- attr(x, 'w_10')
  
  # finally, create indicators for treatment and control
  treat <- which(treatment == 1)
  ctrl  <- which(treatment == 0)
  
  factor_plot <- NULL
  if (any((mediator_is_factor == T))) {
    factor_frames <- list()
    for (i in 1:length(x$stopping_methods)) {
      for (m in mediators_factors) {
        
        method <- x$stopping_methods[i]
        weights <- ifelse(!is.na(w_10[,i]), w_10[,i], w_00[,i])
        
        weights[treat] <- (weights[treat] / sum(weights[treat]))
        weights[ctrl] <- (weights[ctrl] / sum(weights[ctrl]))
        
        a1 <- aggregate(weights[treat], by = list(m = factor(mediators[treat, m])), FUN = sum)
        a0 <- aggregate(weights[ctrl], by = list(m = factor(mediators[ctrl, m])), FUN = sum)
        
        a1[x$a_treatment] <- 1; a0[x$a_treatment] <- 0 ##kat: check this works
        
        combined <- rbind(a1, a0)
        combined['mediator'] <- m
        combined['method'] <- method
        
        factor_frames[[paste(method, m, sep='')]] <- combined
      }
    }
    factor_data <- do.call(rbind, factor_frames)
    trt <- factor_data[,x$a_treatment]
    factor_plot <- lattice::barchart(x ~ factor(m) | factor(method) + factor(mediator),
                                     groups = factor(trt, levels = c(0, 1), labels = c('control', 'treatment')),
                                     data = factor_data,                                   
                                     origin = 0, 
                                     par.settings = list(superpose.polygon = list(col = c("#478BB8", "#B87447"))),
                                     auto.key = TRUE,
                                     ylab = "Proportion",
                                     xlab = "Weighted Mediator",
                                     main = "Weighted Mediation Histogram, Control vs. Treatment",
                                     horiz = FALSE)
  }
  
  number_plot <- NULL
  if (any((mediator_is_factor == F))) {
    number_frames <- list()
    for (i in 1:length(x$stopping_methods)) {
      for (m in mediators_numbers) {
        
        method <- x$stopping_methods[i]
        weights <- ifelse(!is.na(w_10[,i]), w_10[,i], w_00[,i])
        
        weights[treat] <- (weights[treat] / sum(weights[treat]))
        weights[ctrl] <- (weights[ctrl] / sum(weights[ctrl]))

        number_frames[[paste(method, m, sep='')]] <- data.frame('m' = mediators[,m],
                                                                'mediator' = m,
                                                                'A' = treatment,
                                                                'weights' = weights,
                                                                'method' = method)
      }
    }
    number_data <- do.call(rbind, number_frames)
    number_plot <- lattice::densityplot(~m | factor(method) + factor(mediator),
                                        groups = factor(A, levels = c(0, 1), labels = c('control', 'treatment')),
                                        data = number_data,
                                        weights = weights,
                                        plot.points = FALSE,
                                        origin = 0,
                                        par.settings = list(superpose.line = list(lwd = 3, col = c("#478BB8", "#B87447"))),
                                        auto.key = TRUE,
                                        ylab = "Density",
                                        xlab = "Weighted Mediator",
                                        main = "Weighted Mediation Density, Control vs. Treatment")
    
  }
  
  if (!is.null(factor_plot) && !is.null(number_plot)) {

    # NOTE :: Using `grid.arrange()`` is very slow, so we just plot the 
    #         images in the function here, rather than returning.
    #
    # new_plot <- grid.arrange(factor_plot, number_plot)
    # return(new_plot)
    plot(factor_plot, split = c(1, 1, 1, 2))
    plot(number_plot, split = c(1, 2, 1, 2), newpage = FALSE)
  } else if (!is.null(factor_plot)) {
    plot(factor_plot)
  } else {
    plot(number_plot)
  }
}

