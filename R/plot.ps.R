#' Plots for `ps` objects
#'
#' This function produces a collection of diagnostic plots for `ps` objects.
#'
#' This function produces lattice-style graphics of diagnostic plots.
#'
#' @param x A `ps` object.
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
#' @param subset If multiple `stop.method` rules were used in the `ps()` call,
#'   `subset` restricts the plots of a subset of the stopping rules that were
#'   employed. This argument expects a subset of the integers from 1 to k,
#'   if k `stop.method`s were used. 
#' @param color If `color = FALSE`, figures will be gray scale. Default: `TRUE`.
#' @param ... Additional arguments.
#'
#' @references Dan McCaffrey, G. Ridgeway, Andrew Morral (2004). "Propensity
#'   Score Estimation with Boosted Regression for Evaluating Adolescent
#'   Substance Abuse Treatment", *Psychological Methods* 9(4):403-425.
#'
#' @seealso [ps]
#'
#' @method plot ps
#' @export
#' @md
plot.ps <- function(x, plots = "optimize", subset = NULL, color = TRUE, ...)
{
   # Creates diag.plot plots and sends to current device
   # x:     ps object 
   # label: Label added to the plot titles

   # extract the propensity scores and weights from the ps object
   pt1 <- diag.plot.color(x, plots, subset = subset, color = color, ...)


return(pt1)

}
