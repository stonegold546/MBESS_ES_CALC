\name{conf.limits.nct.M3}
\alias{conf.limits.nct.M3}

\title{Confidence limits for a noncentrality parameter from a t-distribution (Method 3 of 3)}

\description{
Function to determine the noncentrality parameters necessary to form a confidence interval around the population noncentrality parameter and related parameters.
Method 3 uses an iterative scheme in order to determine the critical values.
}

\usage{
conf.limits.nct.M3(ncp, df, conf.level = NULL, alpha.lower = NULL, alpha.upper = NULL, tol = 1e-09, sup.int.warns = TRUE, max.steps = 2500, ...)
}
\arguments{
 \item{ncp}{the noncentrality parameter (e.g., observed t-value) of interest}
  \item{df}{the degrees of freedom}
  \item{conf.level}{the level of confidence for a symmetric confidence interval}
  \item{alpha.lower}{the proportion of values beyond the lower limit of the confidence interval (cannot be used with \code{conf.level})}
  \item{alpha.upper}{the proportion of values beyond the upper limit of the confidence interval (cannot be used with \code{conf.level})}
  \item{tol}{is the tolerance of the iterative method for determining the critical values}
  \item{sup.int.warns}{Suppress internal warnings (from internal functions): TRUE or FALSE}
  \item{max.steps}{maximum number of iterations when finding the lower and upper confidence limit}
  \item{\dots}{ allows one to potentially include parameter values for inner functions}
  }


\value{ 
  \item{Lower.Limit}{Value of the distribution with \code{Lower.Limit} noncentral value that has at its specified quantile \code{F.value}}
  \item{Prob.Less.Lower}{Proportion of the distribution beyond (i.e., less than) \code{Lower.Limit}}
  \item{Upper.Limit}{Value of the distribution with \code{Upper.Limit} noncentral value that has at its specified quantile \code{F.value}}
  \item{Prob.Greater.Upper}{Proportion of the distribution beyond (i.e., larger than) \code{Upper.Limit}}
  }


\note{This code was adapted from code written by Ken Kelley and Joseph R. Rausch (University of Notre Dame; \email{JRausch@ND.Edu}), which was adapted from
code written by Michael Smithson (Australian National University; \email{Michael.Smithson@ANU.Edu.AU}).}

\author{Ken Kelley (Indiana University; \email{KKIII@Indiana.Edu})}

\seealso{\code{conf.limits.nct}}

\keyword{models}
\keyword{htest}

