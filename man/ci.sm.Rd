\name{ci.sm}
\alias{ci.sm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Confidence Interval for the Standardized Mean }
\description{
  Function to obtain the exact confidence interval for the standardized mean.
}
\usage{


ci.sm(sm = NULL, Mean = NULL, SD = NULL, ncp = NULL, N = NULL, conf.level = 0.95,
 alpha.lower = NULL, alpha.upper = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sm}{ standardized mean}
  \item{Mean}{ mean }
  \item{SD}{ standard deviation}
  \item{ncp}{ noncentral parameter }
  \item{N}{ sample size }
  \item{conf.level}{ confidence interval coverage (i.e., 1 - Type I error rate); default is .95 }
  \item{alpha.lower}{ Type I error for the lower confidence limit }
  \item{alpha.upper}{ Type I error for the upper confidence limit }
  \item{\dots}{ allows one to potentially include parameter values for inner functions }
}
\details{

The user must specify the standardized mean in one and only one of the three ways: a) mean and standard deviation ('Mean' and 'SD'), b) standardized mean ('sm'), and c) noncentral parameter ('ncp'). The confidence level must be specified in one of following two ways: using confidence interval coverage ('conf.level'), or lower and upper confidence limits ('alpha.lower' and 'alpha.upper').

This function uses the exact confidence interval method based on noncentral t distribution. The confidence interval for noncentral t parameter can be abtained from function 'conf.limits.nct' in MBESS.
}
\value{
  
  \item{Lower.Conf.Limit.Standardized.Mean }{lower confidence limit of the standardized mean}
  \item{Standardized.Mean }{standardized mean}
  \item{Upper.Conf.Limit.Standardized.Mean }{upper confidence limit of the standardized mean}

}
\references{
Kelley, K. (2006). Constructing confidence intervals for standardized effect sizes: Theory, 
application, and implementation. \emph{Manuscript submitted for publication}.
}
\author{ Ken Kelley (Indiana University;<KKIII@Indiana.Edu>) }
\note{ 
The standardized mean is the mean divided by the standard deviation.
}
\seealso{ 'conf.limits.nct' }
\examples{
ci.sm(sm=2.037905, N=13, conf.level=.95)
ci.sm(Mean=30, SD=14.721, N=13, conf.level=.95)
ci.sm(ncp=7.347771, N=13, conf.level=.95)
ci.sm(sm=2.037905, N=13, alpha.lower=.05, alpha.upper=0)
ci.sm(Mean=50, SD=10, N=25, conf.level=.95)
}
\keyword{ design }% at least one, from doc/KEYWORDS
