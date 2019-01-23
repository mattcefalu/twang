\name{raceprofiling}
\alias{raceprofiling}
\docType{data}
\title{Traffic stop data}
\description{
 Simulated example data for assessing race bias in traffic stop outcomes
}
\usage{data(raceprofiling)}
\format{
  A data frame with 5000 observations on the following 10 variables.
  \describe{
    \item{\code{id}}{an ID for each traffic stop}
    \item{\code{nhood}}{a factor indicating the neighborhood in which the stop 
        occurred.}
    \item{\code{reason}}{The reason for the stop, mechanical/registration 
        violations, dangerous moving violation, non-dangerous moving violation}
    \item{\code{resident}}{an indicator whether the driver is a resident of the 
        city}
    \item{\code{age}}{driver's age}
    \item{\code{male}}{an indicator whether the driver was male}
    \item{\code{race}}{the race of the driver, with levels \code{A}, \code{B}, 
        \code{H}, \code{W}}
    \item{\code{hour}}{the hour of the stop (24-hour clock)}
    \item{\code{month}}{and ordered factor indicating in which month the stop 
        took place}
    \item{\code{citation}}{an indicator of whether the driver received a citation}
  }
}

\source{
This is simulated data to demonstrate how to use \code{twang} to adjust 
estimates of racial bias for important factors. This dataset does not represent
real data from any real law enforcement agency.
}
\references{
G. Ridgeway (2006). \dQuote{Assessing the effect of race bias in 
post-traffic stop outcomes using propensity scores,} \emph{Journal of 
Quantitative Criminology} 22(1). 

\url{http://www.i-pensieri.com/gregr/rp.shtml}
} 

\examples{
data(raceprofiling)

# the first five lines of the dataset
raceprofiling[1:5,]
}
\keyword{datasets}
