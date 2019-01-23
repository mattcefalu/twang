\name{desc.wts}
\alias{desc.wts}

\title{Diagnosis of weights}
\description{
\code{desc.wts} assesses the quality of a set of weights on balancing a treatment
and control group.
}
\usage{
desc.wts(data, 
         w, 
         sampw = sampw,
         vars = NULL, 
         treat.var, 
         tp, 
         na.action = "level",
         perm.test.iters=0,
         verbose=TRUE,
         alerts.stack,
         estimand, multinom = FALSE, fillNAs = FALSE)
}

\arguments{
  \item{data}{ a data frame containing the dataset }
  \item{w}{ a vector of weights equal to \code{nrow(data)} }
  \item{sampw}{sampling weights, if provided}
  \item{vars}{ a vector of variable names corresponding to \code{data} }
  \item{treat.var}{ the name of the treatment variable }
  \item{tp}{a title for the method ``type" used to create the weights, used
            to label the results}
  \item{na.action}{a string indicating the method for handling missing data}
  \item{perm.test.iters}{an non-negative integer giving the number of iterations
        of the permutation test for the KS statistic. If \code{perm.test.iters=0}
        then the function returns an analytic approximation to the p-value. This
        argument is ignored is \code{x} is a \code{ps} object. Setting
        \code{perm.test.iters=200} will yield precision to within 3\% if the true
        p-value is 0.05. Use \code{perm.test.iters=500} to be within 2\%}
  \item{verbose}{ if TRUE, lots of information will be printed to monitor the
                 the progress of the fitting }
  \item{alerts.stack}{ an object for collecting warnings issued during the analyses }
  \item{estimand}{the estimand of interest: either "ATT" or "ATE"}
  \item{multinom}{Indicator that weights are from a propsensity score analysis with 3 or more treatment groups.}
  \item{fillNAs}{If \code{TRUE} fills NAs with zeros.}
}
\details{
\code{desc.wts} calls \code{\link{bal.stat}} to assess covariate balance.
If \code{perm.test.iters>0} it will call \code{\link{bal.stat}} multiple
times to compute Monte Carlo p-values for the KS statistics and the maximum KS
statistic. It assembles the results into a list object, which usually becomes 
the \code{desc} component of ps objects that \code{\link{ps}} returns.
}
\value{
See the description of the \code{desc} component of the \code{ps} object that 
\code{\link{ps}} returns
}

\seealso{ \code{\link{ps}} }

\keyword{ models }
