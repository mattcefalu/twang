\name{summary.mnps}
\alias{summary.mnps}

\title{ Summarize an mnps object}
\description{
Computes summary information about a stored \code{mnps} object
}
\usage{
\method{summary}{mnps}(object, ...)
}

\arguments{
  \item{object}{ a \code{\link{ps}} object }
  \item{\dots}{ additional arguments affecting the summary produced }
}
\details{
Compresses the information in the \code{desc} component of the \code{ps} object
into a short summary table describing the size of the dataset and the quality of
the propensity score weights.
}
\value{
See \code{\link{ps}} for details on the returned table
}

\seealso{\code{\link{ps}}, \code{\link{mnps}}}

\keyword{ models }
