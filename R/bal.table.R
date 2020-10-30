#' Compute the balance table.
#'
#' Extract the balance table from [ps], [dx.wts], [mnps], and [mediation] objects
#'
#' `bal.table` is a generic function for extracting balance 
#' tables from [ps], [dx.wts], and [mediation] objects. These objects 
#' usually have several sets of candidate weights, one for an unweighted 
#' analysis and perhaps several `stop.methods`. `bal.table`
#' will return a table for each set of weights combined into a list. Each list 
#' component will be named as given in the `x`, usually the name of the 
#' `stop.method`. The balance table labeled \dQuote{unw} indicates the 
#' unweighted analysis.
#'
#' @param x A [ps], [dx.wts], or [mediation] object.
#' @param digits The number of digits that the numerical entries
#'   should be rounded to. Default: 3.
#' @param collapse.to For `mnps` ATE objects, the comparisons
#'   can be given for all pairs (default), summarized by pre-treatment
#'   covariate and stop.method, or as a single summary for each stop.method.
#' @param subset.var Eliminate all but a specified subset of covariates.
#' @param subset.treat Subset to either all pairs that include a specified
#'   treatment or a single pair of treatments.
#' @param subset.stop.method Subset to either all pairs that include a specified
#'   treatment or a single pair of treatments.
#' @param es.cutoff Subsets to comparisons with absolute ES values bigger than
#'   `es.cutoff`. Default: 0.
#' @param ks.cutoff Subsets to comparisons with KS values bigger than
#'   `ks.cutoff`. Default: 0.
#' @param p.cutoff Subsets to comparisons with t- or chi-squared p-values
#'   no bigger than `p.cutoff`. Default: 1.
#' @param ks.p.cutoff Subsets to comparisons with t- or chi-squared p-values
#'   no bigger than `p.cutoff`. Default: 1.
#' @param timePeriods Used to subset times for iptw fits.
#' @param ... Additional arugments.
#'
#' @return Returns a data frame containing the balance information.
#'   * `tx.mn` The mean of the treatment group.
#'   * `tx.sd` The standard deviation of the treatment group.
#'   * `ct.mn` The mean of the control group.
#'   * `ct.sd` The standard deviation of the control group.
#'   * `std.eff.sz` The standardized effect size, (tx.mn-ct.mn)/tx.sd.
#'     If tx.sd is small or 0, the standardized effect size can be large or INF. 
#'     Therefore, standardized effect sizes greater than 500 are set to NA.
#'   * `stat` The t-statistic for numeric variables and the chi-square 
#'     statistic for continuous variables.
#'   * `p` The p-value for the test associated with `stat`
#'     `ks` The KS statistic.
#'   * `ks.pval` The KS p-value computed using the analytic approximation,
#'     which does not necessarily work well with a lot of ties.
#'
#' @export
bal.table <- function(x,
                      digits = 3,
                      collapse.to = c("pair","covariate","stop.method")[1],
                      subset.var = NULL,
                      subset.treat = NULL,
                      subset.stop.method = NULL,
                      es.cutoff = 0,
                      ks.cutoff = 0,
                      p.cutoff = 1,
                      ks.p.cutoff = 1,
                      timePeriods = NULL,
                      ...){
	if(!(class(x) %in% c("mnps", "iptw", "mniptw","mediation"))){
   bal.tab <- bal.table.ps(x, digits = digits)
   return(bal.tab)
   }
   else if(class(x) == "mediation"){
   bal.tab <- bal.table.mediation(x,digits=digits)
   return(bal.tab)
	}
   else if(class(x) == "iptw"){
   	if(is.null(timePeriods)) timePeriods <- 1:length(x$psList)
   	for(i in timePeriods){
   		cat("Balance at time ", x$uniqueTimes[i], ":\n")
   		print(bal.table.ps(x$psList[[i]], digits = digits))
   		cat("\n")
   	}
   }
   else if(class(x) == "mnps"){
   	bal.table.mnps(x=x, digits = digits, collapse.to = collapse.to, subset.var = subset.var,
                     subset.treat = subset.treat, subset.stop.method = subset.stop.method,
                     es.cutoff = es.cutoff, p.cutoff = p.cutoff, ks.p.cutoff = ks.p.cutoff, ...)
   }
   else if(class(x) == "mniptw"){
   	if(is.null(timePeriods)) timePeriods <- 1:length(x$psList)
   	for(i in timePeriods){
   		cat("Balance at time ", x$uniqueTimes[i], ":\n")
   		print(bal.table.mnps(x$psList[[i]], digits = digits, collapse.to = collapse.to,
                              subset.var = subset.var, subset.treat = subset.treat,
                              subset.stop.method = subset.stop.method, es.cutoff = es.cutoff,
                              p.cutoff = p.cutoff, ks.p.cutoff = ks.p.cutoff, ...))
   		cat("\n")
   		
   	}
   }
   	
}


