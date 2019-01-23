\name{bal.table}
\alias{bal.table}


\title{Compute balance table}
\description{
 Extract the balance table from \code{\link{ps}}, 
 \code{\link{dx.wts}}, and \code{\link{mnps}} objects
}
\usage{
bal.table(x, digits = 3, collapse.to = c("pair","covariate","stop.method")[1], 
	subset.var = NULL, subset.treat = NULL, subset.stop.method = NULL, es.cutoff = 0, 
	ks.cutoff = 0, p.cutoff = 1, ks.p.cutoff = 1, timePeriods = NULL, \dots)
}

\arguments{
  \item{x}{a \code{\link{ps}} or \code{\link{dx.wts}} object}
  \item{digits}{The number of digits that the numerical entries should be rounded to.}
  \item{collapse.to}{For \code{mnps} ATE objects, the comparisons can be given for all pairs (default), summarized by 
  	pre-treatment covariate and stop.method, or as a single summary for each stop.method. }
  \item{subset.var}{Eliminate all but a specified subset of covariates.}
  \item{subset.treat}{Subset to either all pairs that include a specified treatment or a single pair of treatments.}
  \item{subset.stop.method}{Subset to a subset of stop.method\'s used to fit the ps object.}
  \item{es.cutoff}{Subsets to comparisons with absolute ES values bigger than \code{es.cutoff}.}
  \item{ks.cutoff}{Subsets to comparisons with KS values bigger than \code{ks.cutoff}.}
  \item{p.cutoff}{Subsets to comparisons with t- or chi-squared p-values no bigger than \code{p.cutoff}.}
  \item{ks.p.cutoff}{Subsets to comparisons with KS p-values no bigger than \code{ks.p.cutoff}.}
  \item{timePeriods}{Used to subset times for iptw fits. }
  \item{\dots}{Additional arugments.}
}

\details{ \code{bal.table} is a generic function for extracting balance 
tables from \code{\link{ps}} and \code{\link{dx.wts}} objects. These objects 
usually have several sets of candidate weights, one for an unweighted 
analysis and perhaps several \code{stop.methods}. \code{bal.table} 
will return a table for each set of weights combined into a list. Each list 
component will be named as given in the \code{x}, usually the name of the 
\code{stop.method}. The balance table labeled \dQuote{unw} indicates the 
unweighted analysis.}

\value{
    Returns a data frame containing the balance information.            
  \item{tx.mn}{The mean of the treatment group}
  \item{tx.sd}{The standard deviation of the treatment group}
  \item{ct.mn}{The mean of the control group}
  \item{ct.sd}{The standard deviation of the control group}
  \item{std.eff.sz}{The standardized effect size, (tx.mn-ct.mn)/tx.sd.  
If tx.sd is small or 0, the standardized effect size can be large or INF. 
Therefore standardized effect sizes greater than 500 are set to NA}
  \item{stat}{the t-statistic for numeric variables and the chi-square 
              statistic for continuous variables}
  \item{p}{the p-value for the test associated with \code{stat}}
  \item{ks}{the KS statistic}
  \item{ks.pval}{the KS p-value computed using the analytic approximation,
                 which does not necessarily work well with a lot of ties}  
}

\seealso{ The example for \code{\link{ps}} contains an example of the 
use of \code{bal.table}}

\keyword{ models }
