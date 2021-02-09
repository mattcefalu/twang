#' Plot `dxwts`
#' 
#' @param x An `dxwts` object.
#' @param plots An indicator of which type of plot is desired. The options are
#' * `"optimize" or 1` A plot of the balance criteria as a function of the GBM 
#'   iteration.
#' * `"boxplot" or 2` Boxplots of the propensity scores for the treatment and 
#'   control cases
#' * `"es" or 3` Plots of the standardized effect size of the pre-treatment 
#'   variables before and after reweighing
#' * `"t" or 4` Plots of the p-values from t-statistics comparing means of 
#'   treated and control subjects for pretreatment variables, before and after 
#'   weighting.
#' * `"ks" or 5` Plots of the p-values from Kolmogorov-Smirnov statistics 
#'   comparing distributions of pretreatment variables of treated and control 
#'   subjects, before and after weighting.
#' @param ... Additional arguments.
#'
#' @method plot dxwts
#' @export
plot.dxwts <- function(x, plots="es", ...){
  pt1 <- diag.plot(x, plots, ...)
  return(pt1)	
}

