\name{ci.c}
\alias{ci.c}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Confidence interval for a contrast in a fixed effets ANOVA }
\description{
Function to calculate the exact confidence interval for a contrast in a fixed effects analysis of variance context.
}
\usage{
ci.c(means = NULL, error.variance = NULL, c.weights = NULL, n = NULL, 
N = NULL, Psi = NULL, conf.level = 0.95, alpha.lower = NULL, 
alpha.upper = NULL, df.error = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{means}{a vector of the group means or the means of the particular level of the effect (for fixed effect designs) }
  \item{error.variance}{ the common variance of the error (i.e., the mean square error) }
  \item{c.weights}{ the contrast weights (the sum of the contrast weights should be zero) }
  \item{n}{ sample sizes per group or level of the particular factor (if length 1 it is assumed that the per group/level sample sizes are equal) }
  \item{N}{ total sample size }
  \item{Psi}{ the (unstandardized) contrast effect, obtained by multiplying the \emph{j}th mean by the jth contrast weight (this is the unstandardized effect) }
  \item{conf.level}{ confidence interval coverage (i.e., 1- Type I error rate); default is .95 }
  \item{alpha.lower}{ Type I error for the lower confidence limit }
  \item{alpha.upper}{ Type I error for the upper confidence limit}
  \item{df.error}{ the degrees of freedom for the error. In one-way designs, this is simply \emph{N}-length (means) and need not be specified; it must be specified if the design has multiple factors. }
  \item{\dots}{ allows one to potentially include parameter values for inner functions }
}
%\details{}
\value{
Returns the confidence limits for the contrast: 
  
  \item{Lower.Conf.Limit.Contrast }{The lower confidence limit for the contrast effect}
  \item{Contrast }{the value of the estimated unstandardized contrast effect}
  \item{Upper.Conf.Limit.Contrast}{The upper confidence limit for the contrast effect}
}
\references{
Kelley, K. (2007). Constructing confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20} (8), 1-24.

Steiger, J. H. (2004). Beyond the \emph{F} Test: Effect size confidence intervals and tests of close fit in the 
Analysis of Variance and Contrast Analysis. \emph{Psychological Methods, 9}, 164--182.
}

\author{ Ken Kelley (Indiana University, <KKIII@Indiana.edu>) }

\note{ Be sure to use the error varaince and not its square root (i.e., the standard deviation of the errors).
}

\seealso{ \code{conf.limits.nct}, \code{ci.sc}, \code{ci.src}, \code{ci.smd}, \code{ci.smd.c}, \code{ci.sm}}

\examples{
ci.c(means=c(2, 4, 9, 13), error.variance=1, c.weights=c(1, -1, -1, 1), 
n=c(3, 3, 3, 3), N=12, conf.level=.95)

ci.c(means=c(2, 4, 9, 13), error.variance=1, c.weights=c(1, -1, -1, 1), 
n=c(3, 3, 3, 3), N=12, conf.level=.95)

ci.c(means=c(1.6, 0), error.variance=1, c.weights=c(1, -1), n=c(10, 10), 
N=20, conf.level=.95)

# An example given by Maxwell and Delaney (2004, pp. 155--171) :
# 24 subjects of mild hypertensives are assigned to one of four treatments: drug 
# therapy, biofeedback, dietary modification, and a treatment combining all the 
# three previous treatments. Subjects' blood pressure is measured two weeks
# after the termination of treatment. Now we want to form a 95% level
# confidence interval for the difference in blood pressure between subjects
# received drug treatment and those received biofeedback treatment 

## Drug group's mean = 94; group size=4
## Biofeedback group's mean = 91; group size=6 
## Diet group's mean = 92; group size=5
## Combination group's mean = 83; group size=5
## Mean Square Within (i.e., 'error.variance') = 67.375

ci.c(means=c(94, 91, 92, 83), error.variance=67.375, c.weights=c(1, -1, 0, 0), 
n=c(4, 6, 5, 5), N=20, conf.level=.95)

}
\keyword{design}% at least one, from doc/KEYWORDS

