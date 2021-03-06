% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bal.table.R
\name{bal.table}
\alias{bal.table}
\title{Compute the balance table.}
\usage{
bal.table(
  x,
  digits = 3,
  collapse.to = c("pair", "covariate", "stop.method")[1],
  subset.var = NULL,
  subset.treat = NULL,
  subset.stop.method = NULL,
  es.cutoff = 0,
  ks.cutoff = 0,
  p.cutoff = 1,
  ks.p.cutoff = 1,
  timePeriods = NULL,
  ...
)
}
\arguments{
\item{x}{A \link{ps} or \link{dx.wts} object.}

\item{digits}{The number of digits that the numerical entries
should be rounded to. Default: 3.}

\item{collapse.to}{For \code{mnps} ATE objects, the comparisons
can be given for all pairs (default), summarized by pre-treatment
covariate and stop.method, or as a single summary for each stop.method.}

\item{subset.var}{Eliminate all but a specified subset of covariates.}

\item{subset.treat}{Subset to either all pairs that include a specified
treatment or a single pair of treatments.}

\item{subset.stop.method}{Subset to either all pairs that include a specified
treatment or a single pair of treatments.}

\item{es.cutoff}{Subsets to comparisons with absolute ES values bigger than
\code{es.cutoff}. Default: 0.}

\item{ks.cutoff}{Subsets to comparisons with KS values bigger than
\code{ks.cutoff}. Default: 0.}

\item{p.cutoff}{Subsets to comparisons with t- or chi-squared p-values
no bigger than \code{p.cutoff}. Default: 1.}

\item{ks.p.cutoff}{Subsets to comparisons with t- or chi-squared p-values
no bigger than \code{p.cutoff}. Default: 1.}

\item{timePeriods}{Used to subset times for iptw fits.}

\item{...}{Additional arugments.}
}
\value{
Returns a data frame containing the balance information.
\itemize{
\item \code{tx.mn} The mean of the treatment group.
\item \code{tx.sd} The standard deviation of the treatment group.
\item \code{ct.mn} The mean of the control group.
\item \code{ct.sd} The standard deviation of the control group.
\item \code{std.eff.sz} The standardized effect size, (tx.mn-ct.mn)/tx.sd.
If tx.sd is small or 0, the standardized effect size can be large or INF.
Therefore, standardized effect sizes greater than 500 are set to NA.
\item \code{stat} The t-statistic for numeric variables and the chi-square
statistic for continuous variables.
\item \code{p} The p-value for the test associated with \code{stat}
\code{ks} The KS statistic.
\item \code{ks.pval} The KS p-value computed using the analytic approximation,
which does not necessarily work well with a lot of ties.
}
}
\description{
Extract the balance table from \link{ps}, \link{dx.wts}, and \link{mnps} objects
}
\details{
\code{bal.table} is a generic function for extracting balance
tables from \link{ps} and \link{dx.wts} objects. These objects
usually have several sets of candidate weights, one for an unweighted
analysis and perhaps several \code{stop.methods}. \code{bal.table}
will return a table for each set of weights combined into a list. Each list
component will be named as given in the \code{x}, usually the name of the
\code{stop.method}. The balance table labeled \dQuote{unw} indicates the
unweighted analysis.
}
