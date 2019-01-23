\name{get.weights.unstab}
\alias{get.weights.unstab}

\title{ Extract unstabilized propensity score weights for iptw fits. }
\description{
  Extracts propensity score weights from an iptw or mniptw object.
}
\usage{
get.weights.unstab(x,
			stop.method = NULL,
            withSampW = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A \code{iptw} or \code{mniptw} object.}
  \item{stop.method}{The stop method used for the fit of interest. }
  	\item{withSampW}{Returns weights with sample weights multiplied in, if they were provided in the original \code{iptw} call.}
}
\details{

  Weights are the reciprocal of the product of the probability of receiving the treatment received.
}

\value{
a data.frame of weights
}

\seealso{ \code{\link{iptw}} }

\keyword{ utilities }
