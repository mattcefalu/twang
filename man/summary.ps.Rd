% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.ps.R
\name{summary.ps}
\alias{summary.ps}
\title{Summarize a \code{ps} object}
\usage{
\method{summary}{ps}(object, ...)
}
\arguments{
\item{object}{An \code{ps} object.}

\item{...}{Additional arguments.}
}
\value{
See \link{ps} for details on the returned table.
}
\description{
Computes summary information about a stored \code{ps} object
}
\details{
Compresses the information in the \code{desc} component of the \code{ps} object
into a short summary table describing the size of the dataset and the quality of
the propensity score weights.
}
\seealso{
\link{ps}
}
\keyword{models}
