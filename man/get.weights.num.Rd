\name{get.weights.num}
\alias{get.weights.num}

\title{ Get numerators to stabilize propensity score weights for iptw fits. }
\description{
  Forms numerators to stabilize weights for an iptw object.
}
\usage{
get.weights.num(iptw,
			fitList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iptw}{An \code{iptw} object.}
  \item{fitList}{A list containing objects with an associated "fitted" function.}
}
\details{

  Returns numerator of stabilized weights to be used in conjunction with \code{get.weights.unstab}
}

\value{
A vector of stabilizing factors for weights.
}

\seealso{ \code{\link{iptw}} }

\keyword{ utilities }
