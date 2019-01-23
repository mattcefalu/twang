\name{means.table}
\alias{means.table}

\title{ Extract table of means from an mnps object }
\description{
  Extracts table of means from an mnps object.
}
\usage{
means.table(mnps,
			stop.method = 1,
            includeSD = FALSE, digits = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mnps}{ An \code{mnps} object. }
  \item{stop.method}{Indicates which set of weights to retrieve from the \code{ps} object.  Either the name of the stop.method used, or a natural number with 1, for example, indicating the first stop.method specified. }  
  \item{includeSD}{Indicates whether standard deviations as well as means are to be displayed.  By default, they are not displayed. }
  \item{digits}{If not \code{NULL}, results will be rounded to the specified number of digits.}
}
\details{
  Displays a table with weighted and unweighted means and standardized effect sizes, and -- if requested -- standard deviations.
}

\value{
A table of means, standardized effect sizes, and perhaps standard deviations, by treatment group.
}

\seealso{ \code{\link{mnps}} }

\keyword{ utilities }
