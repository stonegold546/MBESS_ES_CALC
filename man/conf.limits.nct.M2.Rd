\name{conf.limits.nct.M2}
\alias{conf.limits.nct.M2}

\title{Confidence limits for a noncentrality parameter from a t-distribution (Method 2 of 3)}

\description{
Largely internal function to determine the noncentrality parameters necessary to form a confidence interval around the population noncentrality parameter and related parameters.
Method 2 uses the \code{nlm} function to determine the critical values. This function requires
the \code{ncp} to be positive, but the function that should be used \code{conf.limits.nct}
does not.
}

\usage{
conf.limits.nct.M2(ncp = ncp, df = df, conf.level = NULL, alpha.lower = NULL, 
alpha.upper = NULL, tol = 1e-09, sup.int.warns = TRUE, ...)
}

\arguments{
 \item{ncp}{the noncentrality parameter (e.g., observed t-value) of interest}
  \item{df}{the degrees of freedom}
  \item{conf.level}{the level of confidence for a symmetric confidence interval}
  \item{alpha.lower}{the proportion of values beyond the lower limit of the confidence interval (cannot be used with \code{conf.level})}
  \item{alpha.upper}{the proportion of values beyond the upper limit of the confidence interval (cannot be used with \code{conf.level})}
   \item{tol}{is the tolerance of the iterative method for determining the critical values}
  \item{sup.int.warns}{Suppress internal warnings (from internal functions): TRUE or FALSE}
  \item{\dots}{ allows one to potentially include parameter values for inner functions}
}

\value{ 
  \item{Lower.Limit}{Value of the distribution with \code{Lower.Limit} noncentral value that has at its specified quantile \code{F.value}}
  \item{Prob.Less.Lower}{Proportion of the distribution beyond (i.e., less than) \code{Lower.Limit}}
  \item{Upper.Limit}{Value of the distribution with \code{Upper.Limit} noncentral value that has at its specified quantile \code{F.value}}
  \item{Prob.Greater.Upper}{Proportion of the distribution beyond (i.e., larger than) \code{Upper.Limit}}
  }
  
\author{Ken Kelley (Indiana University; \email{KKIII@Indiana.Edu})}

\seealso{\code{conf.limits.nct}, \code{nlm}}

\keyword{models}
\keyword{htest}
