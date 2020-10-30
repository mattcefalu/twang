#' Plot the `mediation` object.
#'
#' @param x weighted_mediation object
#' @param plots An indicator of which type of plot is desired. The options are
#'   * `"optimize"` A plot of the balance criteria as a function of the GBM
#'     iteration.
#'   * `"boxplot"` Boxplots of the propensity scores for the treatment and
#'     control cases
#'   * `"es"` Plots of the standardized effect size of the pre-treatment
#'     variables before and after reweighing
#'   * `"density"` Distriubtion plots of NIE1 (distribution of mediator for treatment
#'     sample weighted to match distribution of mediator under control for the population)
#'     and NIE0 (distribution of mediator for control sample weighted to match 
#'     distribution of mediator under treatment for the population) for each mediator.
#'     For continuous mediators, distributions are plotted with density curves and 
#'     for categorical (factor) mediators, distributions are plotted with barplots. 
#' @param subset Used to restrict which of the `stop.method`s will be used
#'   in the figure. For example `subset = c(1,3)` would indicate that the
#'   first and third `stop.method`s (in alphabetical order of those specified
#'   in the original call to the mediation function) should be included in the
#'   figure. If x$method = `logistic` or `crossval`, there is no need to subset 
#'   as there is only one method used. 
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
  
  # return error if model subset is not 1, 2 or NULL
  if (!is.null(model_subset) && !(model_subset %in% c(1, 2))) {
    stop("The `model_subset` must be NULL, 1, or 2.")
  }
  
  # return error if ask for any plots other than available 
  if (!plots %in% c("optimize","boxplot","es","density")) {
     stop("The `plots` options must be `optimize`,`boxplot`,`es`,or `density`.")
  }
  # return error if ask for plots="optimize" or 1 for method!=ps
  if (x$method!="ps" & (plots=="optimize" || plots==1)) { 
   stop("The optimize plot is only available for method='ps'.")}

  # we want the mediator and the treatment variables
  mediators <- x$data[,x$mediator_names, drop = F]
  treatment <- x$data[[x$a_treatment]]
  
  if (plots != 'density') {
    args <- list(plots = plots, subset = subset, color = color)
    if(x$method=="logistic") {
       x$model_a$ps  <- data.frame(logistic=predict(x$model_a,type="response"))
       x$model_m0$ps  <- data.frame(logistic=predict(x$model_m0,type="response"))
       x$model_a$w  <- data.frame(logistic=ifelse(x$data[,x$a_treatment]==1,1/x$model_a$ps[,1],1/(1-x$model_a$ps[,1])))
       x$model_m0$w  <- data.frame(logistic=ifelse(x$data[,x$a_treatment]==1,1/exp(predict(x$model_m0,type="link")),1))
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
       x$model_m0$w  <- data.frame(crossval=ifelse(x$data[,x$a_treatment]==1,1/exp(predict(x$model_m0, n.trees=best.iter.m, newdata=x$data, type="link")),1))
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
    model_names <- c('Model A', 'Model M0')
    
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
  
  ##check if any mediators binary and make factors 
  for(i in 1:length(x$mediator_names)) {
    if(length(unique(x$data[,x$mediator_names[i]]))<=5 & mediator_is_factor[i]==FALSE) {
      mediator_is_factor[i] <- TRUE
      warning(paste("Mediator",x$mediator_names[i],"is being treated like a factor to plot distributions
       with barplots instead of with density curves."))
	}
  }

  mediators_factors <- x$mediator_names[mediator_is_factor == T, drop = T]
  mediators_numbers <- x$mediator_names[mediator_is_factor == F, drop = T]
  
  ## Create a function for generating the plots so it can be replicated for NIE1 and NIE0

  genplot <- function(which_nie){

     if(which_nie == 1){
        # Plot to check that the weighted counterfactual density of mediator from the treatment sample 
        # matches the estimated population density of the mediator under control or M(0)  
        # The population sample is control and its weight is w_00
        # The counterfactual is M(0) from the treatment sample and its weight is w_10
        
        # get the actual weights from the mediation object
        w_pop <- attr(x, 'w_00')
        w_cfac <- attr(x, 'w_10')
      
        # finally, create indicators for treatment and control
        cfac <- which(treatment == 1)
        pop  <- which(treatment == 0)
     
     }else{
        # Plot to check that the weighted counterfactual density of mediator from the control sample 
        # matches the estimated population density of the mediator under treatment or M(1)  
        # The population sample is treatment and its weight is w_11
        # The counterfactual is M(1) from the treatment sample and its weight is w_01
  
        # get the actual weights from the mediation object
        # w_pop is the population wait -- now f
        w_pop <- attr(x, 'w_11')
        w_cfac <- attr(x, 'w_01')
   
        # finally, create indicators for treatment and control
        cfac <- which(treatment == 1)
        pop   <- which(treatment == 0)
        }
        
     if(which_nie==1){
        ptitle <- "NIE1: Distribution of Mediator for Treatment Sample Weighted to Match \n Distribution of Mediator under Control for the Population"
     }else{
        ptitle <- "NIE0: Distribution of Mediator for Control Sample Weighted to Match \n Distribution of Mediator under Treatment for the Population"
     }
  
  if(color) {
    cols <- c("#478BB8", "#B87447")
  } else { cols <- c("black","gray80") }
  stripBgCol <- ifelse(color, "#ffe5cc", "transparent")

  factor_plot <- NULL
  if (any((mediator_is_factor == T))) {
    factor_frames <- list()
    for (i in 1:length(x$stopping_methods)) {
      for (m in mediators_factors) {
        
        method <- x$stopping_methods[i]
        weights <- ifelse(!is.na(w_cfac[,i]), w_cfac[,i], w_pop[,i])
        
        weights[cfac] <- (weights[cfac] / sum(weights[cfac]))
        weights[pop] <- (weights[pop] / sum(weights[pop]))
        
        a1 <- aggregate(weights[cfac], by = list(m = factor(mediators[cfac, m])), FUN = sum)
        a0 <- aggregate(weights[pop], by = list(m = factor(mediators[pop, m])), FUN = sum)
        
        if(which_nie == 1){
              a1[x$a_treatment] <- 1; a0[x$a_treatment] <- 0 
           }else{        
              a1[x$a_treatment] <- 0; a0[x$a_treatment] <- 1 
           }
        
        combined <- rbind(a1, a0)
        combined['mediator'] <- m
        combined['method'] <- method
        
        factor_frames[[paste(method, m, sep='')]] <- combined
      }
    }
    factor_data <- do.call(rbind, factor_frames)
    trt <- factor_data[,x$a_treatment]
    factor_plot <- vector("list",length(mediators_factors))
    pos <- 0
    for(mm in mediators_factors) {
      pos <- pos+1   
      factor_plot[[pos]] <- lattice::barchart(x ~ factor(m) | factor(method) + factor(mediator),
                                     groups = factor(trt, levels = c(0, 1), labels = c('population', 'counterfactual')),
                                     data = factor_data[factor_data$mediator==mm,],                                   
                                     origin = 0, 
                                     par.settings = list(superpose.polygon = list(col = cols),strip.background = list(col=stripBgCol)),
                                     auto.key = TRUE,
                                     ylab = "Proportion",
                                     xlab = "Weighted Mediator",
                                     main = ptitle,
                                     horiz = FALSE,
                                     cex.main=0.85)
  }
  }
  
  number_plot <- NULL
  if (any((mediator_is_factor == F))) {
    number_frames <- list()
    for (i in 1:length(x$stopping_methods)) {
      for (m in mediators_numbers) {
        
        method <- x$stopping_methods[i]
        weights <- ifelse(!is.na(w_cfac[,i]), w_cfac[,i], w_pop[,i])
        
        weights[cfac] <- (weights[cfac] / sum(weights[cfac]))
        weights[pop] <- (weights[pop] / sum(weights[pop]))

        number_frames[[paste(method, m, sep='')]] <- data.frame('m' = mediators[,m],
                                                                'mediator' = m,
                                                                'A' = treatment,
                                                                'weights' = weights,
                                                                'method' = method)
      }
    }
    number_data <- do.call(rbind, number_frames)
    number_plot <- vector("list",length(mediators_numbers))
    pos <- 0
    for(mm in mediators_numbers) {
      pos <- pos+1
      number_plot[[pos]] <- lattice::densityplot(~m | factor(method) + factor(mediator),
                                        groups = factor(A, levels = c(0, 1), labels = c('population', 'counterfactual')),
                                        data = number_data[number_data$mediator==mm,],
                                        weights = weights,
                                        plot.points = FALSE,
                                        origin = 0,
                                        par.settings = list(superpose.line = list(lwd = 3, col = cols),strip.background = list(col=stripBgCol)),
                                        auto.key = TRUE,
                                        ylab = "Density",
                                        xlab = "Weighted Mediator",
                                        main = ptitle,
                                        cex.main=0.55)
    }
  }
  
  if (!is.null(factor_plot) && !is.null(number_plot)) {
    for(i in 1:length(factor_plot)) {
      plot(factor_plot[[i]])
      cc <- par()$ask
      par(ask=TRUE)
    }
      for(i in 1:length(number_plot)) {
      plot(number_plot[[i]])
      cc <- par()$ask
      par(ask=TRUE)
    }
       par(ask=cc)
  } else if (!is.null(factor_plot)) {
    for(i in 1:length(factor_plot)) {
      plot(factor_plot[[i]])
      cc <- par()$ask
      par(ask=TRUE)
    }
     par(ask=cc)
  } else {
      for(i in 1:length(number_plot)) {
      plot(number_plot[[i]])
      cc <- par()$ask
      par(ask=TRUE)
    }
     par(ask=cc)
  }
} # closes function genplot

  cask <- par()$ask

  # Run the density plots for NIE1
  genplot(which_nie=1)

  # Run the density plots for NIE0
  par(ask=TRUE)
  genplot(which_nie=0)
 
  par(ask=cask)
}

