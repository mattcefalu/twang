\name{sensitivity}
\alias{sensitivity}

\title{Sensitivity analyses for propensity score analyses}

\description{ \code{sensitivity} }

\usage{ 
sensitivity(ps1,
   data,
   outcome,
   order.by.importance = TRUE,
   verbose = TRUE)
}

\arguments{
  \item{ps1}{A ps object.}
  \item{data}{The dataset including the outcomes}
  \item{outcome}{The outcome of interest.}
  \item{order.by.importance}{Orders the output by relative importance of covariates.}
  \item{verbose}{If TRUE, extra information will be printed.}
}


\details{
Performs the sensitivity analyses described in Ridgeway (2006). This is a beta version of this functionality. Please let the developers know if you have problems with it. }

\value{ 
	\item{tx}{Summary for treated observations.}
	\item{ctrl}{Summary for control observations.}
 }

\references{
Ridgeway, G. (2006). \dQuote{The effect of race bias in post-traffic stop outcomes using propensity scores,}
\emph{Journal of Quantitative Criminology} 22(1):1-29.
}




