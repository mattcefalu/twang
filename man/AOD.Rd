\name{AOD}
\alias{AOD}
\docType{data}
\title{Subset of Alcohol and Other Drug treatment data}
\description{
A small subset of the data from McCaffrey et al. (2013).
}
\usage{data(AOD)}
\format{
  A data frame with 600 observations on the following 10 variables.
  \describe{
    \item{\code{treat}}{Treatment that each study subject received. Either community, metcbt5, or scy. }
    \item{\code{suf12}}{outcome variable, substance use frequency at 12 month follow-up}
    \item{\code{illact}}{covariate, illicit activities scale}
    \item{\code{crimjust}}{covariate, criminal justice involvement}
    \item{\code{subprob}}{covariate, substance use problem scale}
    \item{\code{subdep}}{covariate, substance use dependence scale}
    \item{\code{white}}{1 if non-Hispanic white, 0 otherwise}

  }
}

\references{
McCaffrey, DF, BA Griffin, D Almirall, ME Slaughter, R Ramchand and LF Burgette (2013). A tutorial on 
propensity score estimation for multiple treatments using generalized boosted models.  Statistics in Medicine.

}

\keyword{datasets}
