\name{ci.sc}
\alias{ci.sc}

\title{ Confidence Interval for a Standardized Contrast in a Fixed Effets ANOVA}
\description{
Function to obtain the confidence interval for a standardized contrast in a fixed effets analysis of variance context.
}
\usage{
ci.sc(means = NULL, error.variance = NULL, c.weights = NULL, n = NULL, 
N = NULL, Psi = NULL, ncp = NULL, conf.level = 0.95, 
alpha.lower = NULL, alpha.upper = NULL, df.error = NULL, ...)
}

\arguments{
  \item{means}{ a vector of the group means or the means of the particular
   level of the effect (for fixed effect designs) }
  
  \item{error.variance}{ the common variance of the error (i.e., the mean square error) }
  \item{c.weights}{ the contrast weights (the sum of the contrast weights should be zero) }
  \item{n}{ sample sizes per group or level of the particular factor (if length 1 it is 
  assumed that the per group/level sample sizes are equal) }
  
  \item{N}{ total sample size }
  \item{Psi}{ the (unstandardized) contrast effect, obtained by multiplying the \emph{j}th mean 
  by the \emph{j}th contrast weight (this is the unstandardized effect) }
  
  \item{ncp}{ the noncentrality paramter from the t-distribution }
  \item{conf.level}{ desired level of confidence for the computed interval (i.e., 1 - the Type I error rate) }
  \item{alpha.lower}{ the Type I error rate for the lower confidence interval limit }
  \item{alpha.upper}{ the Type I error rate for the upper confidence interval limit }
  \item{df.error}{ the degrees of freedom for the error. In one-way designs, this is simply 
  \emph{N}-length (means) and need not be specified; it must be specified if the design has multiple factors. }
  
  \item{\dots}{ optional additional specifications for nested functions }
}
%\details{
}
\value{
  \item{Lower.Conf.Limit.Standardized.Contrast }{the lower confidence limit for the standardized contrast}
  \item{Standardized.contrast }{standardized contrast}
  \item{Upper.Conf.Limit.Standardized.Contrast}{the upper confidence limit for the standardized contrast}
}

\references{ 
Kelley, K. (2007). Constructing confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20} (8), 1-24.

Lai, K., \& Kelley, K. (2007). Sample size planning for standardized ANCOVA and ANOVA 
contrasts: Obtaining narrow confidence intervals. \emph{Manuscript submitted for publication}.

Steiger, J. H. (2004). Beyond the \emph{F} Test: Effect size confidence intervals and tests of close fit in the 
Analysis of Variance and Contrast Analysis. \emph{Psychological Methods, 9}, 164--182.  

}

\author{ Ken Kelley (Indiana University, <KKIII@Indiana.edu>) }

\note{ Be sure to use the error varaince and not its square root (i.e., the standard deviation of the errors).
}
\seealso{ \code{conf.limits.nct}, \code{ci.src}, \code{ci.smd}, \code{ci.smd.c}, \code{ci.sm}, \code{ci.c} }
\examples{
ci.sc(means=c(2, 4, 9, 13), error.variance=1, c.weights=c(1, -1, -1, 1), 
n=c(3, 3, 3, 3), N=12, conf.level=.95)

ci.sc(means=c(2, 4, 9, 13), error.variance=1, c.weights=c(1, -1, -1, 1), 
n=c(3, 3, 3, 3), N=12, conf.level=.95)

ci.sc(means=c(1.6, 0), error.variance=1, c.weights=c(1, -1), n=c(10, 10), 
N=20, conf.level=.95)
}
\keyword{ design}
