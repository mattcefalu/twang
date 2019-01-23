\name{mnps}
\alias{mnps}

\title{Propensity score estimation}

\description{ \code{mnps} calculates propensity scores and diagnoses them using 
a variety of methods, but centered on using boosted logistic regression as 
implemented in \code{\link[gbm]{gbm}} }

\usage{ 
mnps(formula = formula(data),
   data,
   n.trees = 10000,
   interaction.depth = 3,
   shrinkage = 0.01,
   bag.fraction = 1.0,
   perm.test.iters=0,
   print.level = 2,
   iterlim = 1000,
   verbose = TRUE, 
   estimand = "ATE",
   stop.method = "es.max", 
   sampw = NULL, 
   treatATT = NULL, ...)
}

\arguments{
  \item{formula}{A formula for the propensity score model with the treatment
                 indicator on the left side of the formula and the potential
                 confounding variables on the right side.}
  \item{data}{The dataset, includes treatment assignment as well as covariates}
  \item{n.trees}{ number of gbm iterations passed on to \code{\link[gbm]{gbm}} }
  \item{interaction.depth}{ \code{interaction.depth} passed on to
        \code{\link[gbm]{gbm}} }
  \item{shrinkage}{ \code{shrinkage} passed on to \code{\link[gbm]{gbm}} }
  \item{bag.fraction}{\code{bag.fraction} passed on to \code{\link[gbm]{gbm}} }
  \item{perm.test.iters}{a non-negative integer giving the number of iterations
        of the permutation test for the KS statistic. If \code{perm.test.iters=0}
        then the function returns an analytic approximation to the p-value. Setting
        \code{perm.test.iters=200} will yield precision to within 3\% if the true
        p-value is 0.05. Use \code{perm.test.iters=500} to be within 2\%}
  \item{print.level}{ the amount of detail to print to the screen }
  \item{iterlim}{ maximum number of iterations for the direct optimization }
  \item{verbose}{ if TRUE, lots of information will be printed to monitor the
                 the progress of the fitting }
\item{estimand}{The causal effect of interest.  Options are \code{"ATE"} (average treatment effect),
which attempts to estimate the change in the outcome if the treatment were applied to the entire population
versus if the control were applied to the entire population, or \code{"ATT"} (average treatment effect on 
the treated) which attempts to estimate the analogous effect, averaging only over the treated population.}
  \item{stop.method}{A method or methods of measuring and summarizing balance across 
pretreatment variables.  Current options are \code{ks.mean}, \code{ks.max}, \code{es.mean}, 
and \code{es.max}.  \code{ks} refers to the 
Kolmogorov-Smirnov statistic and \code{es} refers to standardized effect size.  These are summarized
across the pretreatment variables by either the maximum (\code{.max}) or the mean (\code{.mean}).}
   \item{sampw}{Optional sampling weights.} 
\item{treatATT}{If the estimand is specified to be ATT, this argument is 
used to specify which treatment condition is considered 'the treated'.  It must be one of the levels of the treatment variable.  It is ignored for ATE analyses.}
   \item{...}{Additional arguments.}
}


\details{
\code{formula} should be something like "treatment ~ X1 + X2 + X3". The
treatment variable should be a variable with three or more levels. There is no need to specify
interaction terms in the formula. \code{interaction.depth} controls the level
of interactions to allow in the propensity score model.

Note that --- unlike earlier versions of \code{twang} --- plotting functions 
are no longer included in the \code{ps()} function. See 
\code{\link{plot}} for details of the plots.}

\value{ 
Returns an object of class \code{mnps}, which consists of the following.  
  \item{psList}{A list of \code{ps} objects. }
\item{nFits}{The number of calls to \code{ps} that were used to form the \code{mnps} object.}
  \item{estimand}{The estimand -- either ATT or ATE -- that was specified in the call to \code{mnps}.}
\item{treatATT}{For ATT fits, the treatment category that is considered "the treated"}
\item{treatLev}{The levels of the treatment variable.}
\item{levExceptTreatAtt}{The levels of the treatment variable, excluding the \code{treatATT} level.}
\item{data}{The data used to fit the model.}
\item{treatVar}{The vector of treatment indicators}
\item{stopMethods}{The \code{stop.method} vector specified in the call to \code{mnps}.}
\item{sampw}{Sampling weights provided to \code{mnps}, if any.}

}      

\references{
Dan McCaffrey, G. Ridgeway, Andrew Morral (2004). \dQuote{Propensity Score Estimation
with Boosted Regression for Evaluating Adolescent Substance Abuse Treatment,}
\emph{Psychological Methods} 9(4):403-425.
}

\author{
Lane Burgette \email{burgette@rand.org},
Beth Ann Griffin \email{bethg@rand.org},
Dan McCaffrey \email{danielm@rand.org}
}

\seealso{ \code{ps} }


\keyword{models}
\keyword{multivariate}
