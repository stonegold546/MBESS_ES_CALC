\name{conf.limits.nct}
\alias{conf.limits.nct}

\title{Confidence limits for a noncentrality parameter from a t-distribution}

\description{
Function to determine the noncentrality parameters necessary to form a confidence interval around the population noncentrality parameter and related parameters.
Due to the difficulties in estimating the necessary values, three different methods are implemented within the present function (\code{conf.limits.nct.M1}, \code{conf.limits.nct.M2}, and \code{conf.limits.nct.M3}) and the best set of results are taken.
}

\usage{
conf.limits.nct(ncp, df, conf.level = 0.95, alpha.lower = NULL, 
alpha.upper = NULL, t.value, tol = 1e-09, sup.int.warns = TRUE, method = "all", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{ncp}{the noncentrality parameter (e.g., observed t-value) of interest.}
  \item{df}{the degrees of freedom.}
  \item{conf.level}{the level of confidence for a symmetric confidence interval.}
  \item{alpha.lower}{the proportion of values beyond the lower limit of the confidence interval (cannot be used with \code{conf.level}).}
  \item{alpha.upper}{the proportion of values beyond the upper limit of the confidence interval (cannot be used with \code{conf.level}).}
  \item{t.value}{alias for \code{ncp}}
  \item{tol}{is the tolerance of the iterative method for determining the critical values.}
  \item{sup.int.warns}{Suppress internal warnings (from internal functions): TRUE or FALSE}
  \item{method}{which of the three methods should be used: "all", "1", "2", "3" ("all" is default).}
  \item{\dots}{ allows one to potentially include parameter values for inner functions}
}
\details{
Function for finding the upper and lower confidence limits for a noncentral parameter from a noncentral \emph{t}-distribution with \code{df} degrees of freedom. 
This function is especially helpful when forming confidence intervals around standardized mean differences (i.e., Cohen's d; Glass's g; Hedges g'), standardized regression coefficients, and 
coefficients of variations. The \code{Lower.Limit} and the \code{Upper.Limit} values correspond to the noncentral parameters of a \emph{t}-distribution with \code{df} degrees of
freedom whose upper and lower tails contain the desired proportion of the curves, respectively.
When \code{ncp} is zero, the \code{Lower.Limit} and \code{Upper.Limit} are simply the desired quantiles of the
central \emph{t}-distribution with \code{df} degrees of freedom.

See the documentation and code for each of the three methods (if interested). Each of the three methods 
should reach the same values for the confidence limits. However, due to the iterative nature of the functions, 
one function may arrive at a more optimal solution. Furthermore, in some situations one (or more) 
functions could fail to find the optimum values, which necessitates the use of multiple methods to (essentially) ensure
that optimal values are found. 
}

\value{ 
  \item{Lower.Limit}{Value of the distribution with \code{Lower.Limit} noncentral value that has at its specified quantile \code{F.value}}
  \item{Prob.Less.Lower}{Proportion of the distribution beyond (i.e., less than) \code{Lower.Limit}}
  \item{Upper.Limit}{Value of the distribution with \code{Upper.Limit} noncentral value that has at its specified quantile \code{F.value}}
  \item{Prob.Greater.Upper}{Proportion of the distribution beyond (i.e., larger than) \code{Upper.Limit}}
  }
\references{
Cumming, G. \& Finch, S. (2001) A primer on the understanding, use, and calculation of confidence intervals that are
based on central and noncentral distributions, \emph{Educational and Psychological Measurement, 61}, 532--574.

Kelley, K. (2005) The effects of nonnormal distributions on confidence intervals around the standardized mean
difference: Bootstrap and parametric confidence intervals, \emph{Educational and Psychological Measurement, 65}, 51--69.

Steiger, J. \& Fouladi, T. (1997) Noncentrality interval estimation and the evaluation of statistical models. In L. Harlow, 
S. Muliak, \& J. Steiger (Eds.), \emph{What if there were no significance tests?}. Mahwah, NJ: Lawrence Erlbaum.

}
\author{Ken Kelley (Indiana University; \email{KKIII@Indiana.Edu})}

\section{Warning}{At the present time, the largest \code{ncp} that R can accurately handle is 37.62.
}

\seealso{'pt', 'qt', 'ci.smd', 'ci.smd.c', 'ss.aipe', 'conf.limits.ncf', 'conf.limits.nc.chisq'}
\examples{

# Suppose observed t-value based on 'df'=126 is 2.83. Finding the lower 
# and upper critical values for the population noncentrality parameter 
# with a symmetric confidence interval with 95\% confidence is given as:
conf.limits.nct(ncp=2.83, df=126, conf.level=.95)

# Modifying the above example so that nonsymmetric confidence intervals 
# can be formed:
conf.limits.nct(ncp=2.83, df=126, alpha.lower=.01, alpha.upper=.04,
conf.level=NULL)
}
\keyword{models}
\keyword{htest}
