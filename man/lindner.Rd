\name{lindner}
\alias{lindner}
\docType{data}
\title{Lindner Center data on 996 PCI patients analyzed by Kereiakes et al. (2000)}
\description{These data are adapted from the \code{lindner} dataset in the \code{USPS} package.
	The description comes from that package, except for the variable \code{sixMonthSurvive}, which 
	is a recode of \code{lifepres}
	
	Data from an observational study of 996 patients receiving an initial
  Percutaneous Coronary Intervention (PCI) at Ohio Heart Health, Christ Hospital,
  Cincinnati in 1997 and followed for at least 6 months by the staff of the Lindner Center.
  The patients thought to be more severely diseased were assigned to treatment with abciximab
  (an expensive, high-molecular-weight IIb/IIIa cascade blocker); in fact, only 298 (29.9
  percent) of patients received usual-care-alone with their initial PCI.}
\usage{data(lindner)}
\format{
  A data frame of 10 variables collected on 996 patients; no NAs.
  \describe{
    \item{lifepres}{Mean life years preserved due to survival for at least 6 months following
       PCI; numeric value of either 11.4 or 0.}
    \item{cardbill}{Cardiac related costs incurred within 6 months of patient's initial PCI;
       numeric value in 1998 dollars; costs were truncated by death for the 26 patients with
       lifepres == 0.}
    \item{abcix}{Numeric treatment selection indicator; 0 implies usual PCI care alone; 1
       implies usual PCI care deliberately augmented by either planned or rescue treatment
       with abciximab.}
    \item{stent}{Coronary stent deployment; numeric, with 1 meaning YES and 0 meaning NO.}
    \item{height}{Height in centimeters; numeric integer from 108 to 196.}
    \item{female}{Female gender; numeric, with 1 meaning YES and 0 meaning NO.}
    \item{diabetic}{Diabetes mellitus diagnosis; numeric, with 1 meaning YES and 0 meaning NO.}
    \item{acutemi}{Acute myocardial infarction within the previous 7 days; numeric, with 1
       meaning YES and 0 meaning NO.}
    \item{ejecfrac}{Left ejection fraction; numeric value from 0 percent to 90 percent.}
    \item{ves1proc}{Number of vessels involved in the patient's initial PCI procedure; numeric
       integer from 0 to 5.}
    \item{sixMonthSurvive}{Survival at six months --- a recoded version of \code{lifepres}.}   
  }
}
\references{
  Kereiakes DJ, Obenchain RL, Barber BL, et al.  Abciximab provides cost
  effective survival advantage in high volume interventional practice.
  \emph{Am Heart J} 2000; \bold{140}: 603-610.

  Obenchain RL. (2009) \bold{USPSinR.pdf}  ../R\_HOME/library/USPS  40 pages.
}

\keyword{datasets}

	