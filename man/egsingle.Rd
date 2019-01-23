\name{egsingle}
\alias{egsingle}
\docType{data}
\title{US Sustaining Effects study}
\description{
A subset of the mathematics scores from the U.S. Sustaining Effects Study. 
The subset consists of information on 1721 students from 60 schools. This 
dataset is available in the \code{mlmRev} package.
}
\usage{data(egsingle)}
\format{
  A data frame with 7230 observations on the following 12 variables.
  \describe{
    \item{\code{schoolid}}{a factor of school identifiers}
    \item{\code{childid}}{a factor of student identifiers}
    \item{\code{year}}{a numeric vector indicating the year of the test}
    \item{\code{grade}}{a numeric vector indicating the student's grade}
    \item{\code{math}}{a numeric vector of test scores on the IRT scale score 
                       metric}
    \item{\code{retained}}{a factor with levels \code{0} \code{1} indicating if
                           the student has been retained in a grade.}
    \item{\code{female}}{a factor with levels \code{Female} \code{Male}}
    \item{\code{black}}{a factor with levels \code{0} \code{1} indicating if the student is Black}
    \item{\code{hispanic}}{a factor with levels \code{0} \code{1} indicating if the student is Hispanic}
    \item{\code{size}}{a numeric vector indicating the number of students enrolled in the school}
    \item{\code{lowinc}}{a numeric vector giving the percentage of low-income students in the school}
    \item{\code{mobility}}{a numeric vector}
  }
}
\source{
Reproduced from the\code{mlmRev} package for use in the
section on nonresponse weighting in the \code{twang} package
vignette. These data are distributed with the HLM software package
(Bryk, Raudenbush, and Congdon, 1996). Conversion to the R format is
described in Doran and Lockwood (2006).
}
\references{
Doran, H.C. and J.R. Lockwood (2006). ``Fitting value-added models in R,''
\emph{Journal of Educational and Behavioral Statistics}, 31(1)
}
\keyword{datasets}
