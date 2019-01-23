\name{iptwExLong}
\alias{iptwExLong}
\docType{data}
\title{Example data for iptw function (long version)}
\description{These data are simulated to demonstrate the \code{iptw} function in the "long" data format.}
\usage{data(lindner)}
\format{
  A list with a covariate matrix and outcomes.
  \describe{
    \item{covariates}{Time-invariant covariates are \code{gender} and \code{age}. The time-varying covariate is \code{use}. The reatment indicator is given by \code{tx}. An individual level identifier is given in \code{ID}, and the time period is \code{time}. }
    \item{outcome}{Vector of post-treatment outcomes.}   
  }
}


\keyword{datasets}

	