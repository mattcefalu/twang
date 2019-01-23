\name{iptw}
\alias{iptw}

\title{Inverse probability of treatment weighting for marginal structural models}

\description{ \code{iptw} uses \code{\link[gbm]{gbm}} to estimate 
	 propensity scores for sequential treatments.}

\usage{ 
iptw(formula,
   data,
   timeInvariant = NULL,
   n.trees = 10000,
   stop.method = "es.max",
   cumulative = TRUE,
   timeIndicators = NULL,
   ID = NULL,
   priorTreatment = TRUE, ...)

}

\arguments{
  \item{formula}{Either a single formula (long format) or a list with formulas}
  \item{data}{The dataset, includes treatment assignment as well as covariates}
  \item{timeInvariant}{An optional formula (with no left-hand variable) specifying time-invariant chararacteristics.}
  \item{n.trees}{ number of gbm iterations passed on to \code{\link[gbm]{gbm}} }
\item{stop.method}{A method or methods of measuring and summarizing balance across 
pretreatment variables.  Current options are \code{ks.mean}, \code{ks.max}, \code{es.mean}, 
and \code{es.max}.  \code{ks} refers to the 
Kolmogorov-Smirnov statistic and \code{es} refers to standardized effect size.  These are summarized
across the pretreatment variables by either the maximum (\code{.max}) or the mean (\code{.mean}).}
\item{cumulative}{If \code{TRUE}, the time t model includes time-varying characteristics from times 1 through t-1.}  
  \item{timeIndicators}{For long format fits, a vector of times for each observation.}
  \item{ID}{For long format fits, a vector of numeric identifiers for unique analytic units.}
  \item{priorTreatment}{For long format fits, includes treatment levels from previous times if \code{TRUE}. This argument is ignored for wide format fits.}

   \item{...}{Additional arguments that are passed to \code{ps} function.  }
}


\details{This function uses generalized boosted models to estimate inverse probability of treatment weights for sequential treatments.}

\value{ 
Returns an object of class \code{iptw}, a list containing 
  \item{psList}{A list of ps objects with length equal to the number of time periods.}
\item{estimand}{The specified estimand.}
  \item{stop.methods}{The stopping rules used to optimize iptw balance.}      
  \item{nFits}{The number of ps objects (i.e., the number of distinct time points.)}
  \item{uniqueTimes}{The unique times in the specified model.}
  }



\seealso{ \code{\link{ps}} }


\keyword{models}
