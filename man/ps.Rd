\name{ps}
\alias{ps}

\title{Propensity score estimation}

\description{ \code{ps} calculates propensity scores and diagnoses them using 
a variety of methods, but centered on using boosted logistic regression as 
implemented in \code{\link[gbm]{gbm}} }

\usage{ 
ps(formula = formula(data),
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
   stop.method = c("ks.mean", "es.mean"), 
   sampw = NULL, 
   multinom = FALSE, ...)
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
   \item{multinom}{Set to true only when called from \code{mnps} function.} 
   \item{...}{Additional arguments.  }
}


\details{
\code{formula} should be something like "treatment ~ X1 + X2 + X3". The
treatment variable should be a 0/1 indicator. There is no need to specify
interaction terms in the formula. \code{interaction.depth} controls the level
of interactions to allow in the propensity score model.

Note that --- unlike earlier versions of \code{twang} --- plotting functions 
are no longer included in the \code{ps()} function. See 
\code{\link{plot}} for details of the plots.}

\value{ 
Returns an object of class \code{ps}, a list containing 
  \item{gbm.obj}{The returned \code{\link[gbm]{gbm}} object}
\item{treat}{The treatment variable.}
  \item{desc}{a list containing balance tables for each method selected in
              \code{stop.methods}. Includes a component for the unweighted
              analysis names \dQuote{unw}. Each \code{desc} component includes
              a list with the following components
  \describe{
     \item{ess}{The effective sample size of the control group}
     \item{n.treat}{The number of subjects in the treatment group}
     \item{n.ctrl}{The number of subjects in the control group}
     \item{max.es}{The largest effect size across the covariates}
     \item{mean.es}{The mean absolute effect size}
     \item{max.ks}{The largest KS statistic across the covariates}
     \item{mean.ks}{The average KS statistic across the covariates}
     \item{bal.tab}{a (potentially large) table summarizing the quality of the 
                    weights for equalizing the distribution of features across 
                    the two groups. This table is best extracted using the
                    \code{\link{bal.table}} method. See the help for 
                    \code{\link{bal.table}} for details on the table's contents}
     \item{n.trees}{The estimated optimal number of \code{\link[gbm]{gbm}} 
                    iterations to optimize the loss function for the associated 
                    \code{stop.methods}}

  \item{ps}{a data frame containing the estimated propensity scores. Each
            column is associated with one of the methods selected in
            \code{stop.methods}}
  \item{w}{a data frame containing the propensity score weights. Each
            column is associated with one of the methods selected in
            \code{stop.methods}. If sampling weights are given then these are
            incorporated into these weights.}
\item{estimand}{The estimand of interest (ATT or ATE).}


     }      
  }
  \item{datestamp}{Records the date of the analysis}
  \item{parameters}{Saves the \code{ps} call}
  \item{alerts}{Text containing any warnings accumulated during the estimation}
\item{iters}{A sequence of iterations used in the GBM fits used by \code{plot} function.}
\item{balance}{The balance measures for the pretreatment covariates, with a column for each
\code{stop.method}.}
\item{n.trees}{Maximum number of trees considered in GBM fit.}
\item{data}{Data as specified in the \code{data} argument.}
 }

\references{
Dan McCaffrey, G. Ridgeway, Andrew Morral (2004). \dQuote{Propensity Score Estimation
with Boosted Regression for Evaluating Adolescent Substance Abuse Treatment,}
\emph{Psychological Methods} 9(4):403-425.
}

\author{
Greg Ridgeway \email{gregr@rand.org},
Dan McCaffrey \email{danielm@rand.org},
Andrew Morral \email{morral@rand.org},
Lane Burgette \email{burgette@rand.org}
}

\seealso{ \code{\link[gbm]{gbm}} }


\keyword{models}
\keyword{multivariate}
