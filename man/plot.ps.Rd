\name{plot.ps}
\alias{plot.ps}

\title{Plots for ps objects}

\description{This function produces a collection of diagnostic plots for ps objects. }

\usage{ 
\method{plot}{ps}(x, plots = "optimize", subset=NULL, color = TRUE, \dots)
}

\arguments{
  \item{x}{A ps object}
  \item{plots}{An indicator of which type of plot is desired.  The options are
\describe{
\item{"optimize" or 1}{A plot of the balance criteria as a function of the GBM iteration}
\item{"boxplot" or 2}{Boxplots of the propensity scores for the treatment and control cases}
\item{"es" or 3}{Plots of the standardized effect size of the pre-treatment variables before and after reweighing}
\item{"t" or 4}{Plots of the p-values from t-statistics comparing means of treated and control subjects for pretreatment variables, before and after weighting.}
\item{"ks" or 5}{Plots of the p-values from Kolmogorov-Smirnov statistics comparing distributions of pretreatment variables of treated and control subjects, before and after weighting.}
\item{"histogram" or 6}{Histogram of weights for treated and control subjects.}
}
}
\item{subset}{If multiple \code{stop.method} rules were used in the \code{ps()} call, \code{subset} restricts the plots of a subset of the stopping rules that were employed.  This argument expects a subset of the integers from 1 to k, if k \code{stop.method}s were used.}
\item{color}{If set to \code{FALSE}, grayscale figures will be produced}
  \item{...}{Additional arguments that may be passed to the underlying \code{lattice} package plotting functions}

}


\details{
This function produces lattice-style graphics of diagnostic plots. 
 }

\references{
Dan McCaffrey, G. Ridgeway, Andrew Morral (2004). \dQuote{Propensity Score Estimation
with Boosted Regression for Evaluating Adolescent Substance Abuse Treatment,}
\emph{Psychological Methods} 9(4):403-425.
}


\seealso{ \code{\link{ps}} }


\keyword{models}
\keyword{multivariate}
