\name{boxplot.ps}
\alias{boxplot.ps}

\title{Boxplots for ps objects}

\description{This function produces a collection of diagnostic plots for ps objects. }

\usage{ 
\method{boxplot}{ps}(x, subset=NULL, color = TRUE, time = NULL, \dots)
}

\arguments{
  \item{x}{A ps object}
\item{subset}{If multiple \code{stop.method} rules were used in the \code{ps()} call, \code{subset} restricts the plots of a subset of the stopping rules that were employed.  This argument expects a subset of the integers from 1 to k, if k \code{stop.method}s were used.}
\item{color}{If set to \code{FALSE}, grayscale figures will be produced}
\item{time}{Used to specify a subset of times for use with the iptw function. Ignored for standard ps fits.}
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

