\name{get.weights}
\alias{get.weights}

\title{ Extract propensity score weights }
\description{
  Extracts propensity score weights from a ps or mnps object.
}
\usage{
get.weights(ps1,
			stop.method = NULL,
            estimand = NULL,
            withSampW = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ps1}{ a \code{ps} or \code{mnps} object }
  \item{stop.method}{ indicates which set of weights to retrieve from the \code{ps} object }  
  \item{estimand}{ indicates whether the weights are for the average treatment effect on the treated (ATT) or the average treatment effect on the population (ATE).  By default, \code{get.weights} will use the estimand used to fit the 
  	\code{ps} object. }
  	\item{withSampW}{Returns weights with sample weights multiplied in, if they were provided in the original \code{ps} or \code{mnps} call.}
}
\details{
  Weights for ATT are 1 for the treatment cases and p/(1-p) for the control cases.

  Weights for ATE are 1/p for the treatment cases and 1/(1-p) for the control cases.
}

\value{
a vector of weights
}

\seealso{ \code{\link{ps}} }

\keyword{ utilities }
