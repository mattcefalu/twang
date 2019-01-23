\name{iptwExWide}
\alias{iptwExWide}
\docType{data}
\title{Example data for iptw function (wide version)}
\description{These data are simulated to demonstrate the \code{iptw} function in the "wide" data format.}
\usage{data(lindner)}
\format{
  A list with a covariate matrix and outcomes.
  \describe{
  	\item{gender}{Gender.}
  	\item{age}{Age.}
  	\item{use0}{Baseline substance use}.
  	\item{use1}{Use following first time period treatment.}
  	\item{use2}{Use following second time period treatment.}
  	\item{tx1}{Treatment indicator (first time period).}
  	\item{tx2}{Treatment indicator (second time period).}
  	\item{tx3}{Treatment indicator (third time period).}
    \item{covariates}{Time-invariant covariates are \code{gender} and \code{age}. The time-varying covariate is \code{use}. The reatment indicator is given by \code{tx}. An individual level identifier is given in \code{ID}, and the time period is \code{time}. }
    \item{outcome}{Post-treatment outcomes.}   
  }
}


\keyword{datasets}

	