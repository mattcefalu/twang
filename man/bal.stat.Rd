\name{bal.stat}
\alias{bal.stat}

\title{ Calculate weighted balance statistics }
\description{
\code{bal.stat} compares the treatment and control subjects by means, standard 
deviations, effect size, and KS statistics
}
\usage{
bal.stat(data, 
         vars = NULL, 
         treat.var, 
         w.all,
         sampw, 
         get.means = TRUE,
         get.ks = TRUE, 
         na.action = "level",
         estimand,
         multinom, fillNAs = FALSE)
}

\arguments{
  \item{data}{ a data frame containing the data }
  \item{vars}{ a vector of character strings with the names of the variables 
               on which the function will assess the balance }
  \item{treat.var}{ the name of the treatment variable }
  \item{w.all}{ observation weights (e.g. propensity score weights, sampling 
                weights, or both) }
  \item{sampw}{ sampling weights.  These are passed in addition to \code{w.all} because the "unweighted" results shoud be adjusted for sample weights (though not propensity score weights).}              
  \item{get.means}{ logical. If \code{TRUE} then \code{bal.stat} will compute means
                 and variances}
  \item{get.ks}{ logical. If \code{TRUE} then \code{bal.stat} will compute KS 
                 statistics}
  \item{na.action}{ a character string indicating how \code{bal.stat} should 
                    handle missing values. Current options are "level", 
                    "exclude", or "lowest" }
  \item{estimand}{ either "ATT" or "ATE"}
  \item{multinom}{\code{TRUE} if used for multinomial propensity scores.}
  \item{fillNAs}{If \code{TRUE}, fills in zeros for missing values.}
}
\details{
    \code{bal.stat} calls 
    auxiliary functions 
    for each variable and assembles the results in 
    a table
}
\value{
\code{get.means} and \code{get.ks} manipulate the inclusion of certain 
columns in the returned result.
}

\seealso{ The example for \code{\link{ps}} contains an example of the 
use of \code{bal.table} }

\keyword{ multivariate }
