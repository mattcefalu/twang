\name{plot.mnps}
\alias{plot.mnps}

\title{Plots for mnps objects}

\description{This function produces a collection of diagnostic plots for ps objects. }

\usage{ 
\method{plot}{mnps}(x, plots = "optimize", pairwiseMax = TRUE, figureRows = NULL,  
	color = TRUE, subset = NULL, treatments = NULL, singlePlot = NULL, 
	multiPage = FALSE, time = NULL, print = TRUE, \dots)
}

\arguments{
  \item{x}{An mnps object.}
  \item{plots}{An indicator of which type of plot is desired.  The options are
\describe{
\item{"optimize" or 1}{A plot of the balance criteria as a function of the GBM 
	iteration}
\item{"boxplot" or 2}{Boxplots of the propensity scores for the treatment and 
	control cases}
\item{"es" or 3}{Plots of the standardized effect size of the pre-treatment 
	variables before and after reweighing}
\item{"t" or 4}{Plots of the p-values from t-statistics comparing means of 
	treated and control subjects for pretreatment variables, before and after 
	weighting.}
\item{"ks" or 5}{Plots of the p-values from Kolmogorov-Smirnov statistics 
	comparing distributions of pretreatment variables of treated and control 
	subjects, before and after weighting.}
}
}
\item{pairwiseMax}{If \code{FALSE}, the plots for the underlying \code{ps} fits 
	will be returned.  Otherwise, pairwise maxima will be returned.}
\item{figureRows}{The number of rows of figures that should be used. If left as \code{NULL}, twang 
	tries to find a reasonable value.}
\item{color}{If \code{color = FALSE} figures will be gray scale.}
\item{subset}{Used to restrict which of the \code{stop.method}s will be used 
	in the figure. For example \code{subset = c(1,3)} would indicate that the 
	first and third \code{stop.method}s (in alphabetical order of those specified 
	in the original call to \code{mnps}) should be included in the figure.}
\item{treatments}{Only applicable when \code{pairwiseMax} is \code{FALSE} and \code{plots} 3, 4, and 5.  
	If left at \code{NULL}, panels for all treatment pairs are created.  If one level of the treatment 
	variable is specified, plots comparing that treatment to all others are produced.  If two levels are
	specified, a comparison for that single pair is produced.}
\item{singlePlot}{For plot calls that produce multiple plots, specifying an integer value of \code{singlePlot} will 
	return only the corresponding plot.  E.g., specifying \code{singlePlot = 2} will return the second plot.}
\item{multiPage}{When multiple frames of a figure are produced, \code{multiPage = TRUE} will print each frame on a different page. This is intended for situations where the graphical output is being saved to a file.}
\item{time}{For use with iptw.}
\item{print}{If FALSE, the figure is returned but not printed.}
\item{\dots}{Additional arguments.}
}


\details{
This function produces lattice-style graphics of diagnostic plots.}

\references{
Dan McCaffrey, G. Ridgeway, Andrew Morral (2004). \dQuote{Propensity Score Estimation
with Boosted Regression for Evaluating Adolescent Substance Abuse Treatment,}
\emph{Psychological Methods} 9(4):403-425.
}


\seealso{ \code{\link{mnps}} }


\keyword{models}
\keyword{multivariate}
