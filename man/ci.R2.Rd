\name{ci.R2}
\alias{ci.R2}
\title{Confidence intervals for the multiple correlation coefficient}
\description{A function to calculate the exact confidence interval for the multiple correlation coefficient.}
\usage{
ci.R2(R2 = NULL, df.1 = NULL, df.2 = NULL, conf.level = NULL, 
F.value = NULL, N = NULL, p = NULL, alpha.lower = NULL, 
alpha.upper = NULL, tol = 1e-09)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R2}{multiple correlation coefficient}
  \item{df.1}{numerator degrees of freedom}
  \item{df.2}{denominator degrees of freedom}
  \item{conf.level}{confidence interval coverage; 1-Type I error rate}
  \item{F.value}{obtained \code{F}-value}
  \item{N}{sample size}
  \item{p}{number of predictors}
  \item{alpha.lower}{Type I error for the lower confidence limit}
  \item{alpha.upper}{Type I error for the upper confidence limit}
  \item{tol}{tolerance for iterative convergence}}

\details{
Specify either \code{R2} or \code{F.value} and \code{N} and \code{p} or \code{R2} and \code{df.1} 
and \code{df.2}. Note that this function uses the noncentral \code{F}-distribution to determine 
critical values.}

\value{
  \item{Lower.Conf.Limit.R2 }{upper limit of the confidence interval around the population multiple correlation coefficient}
  \item{Prob.Less.Lower }{proportion of the distribution less than \code{Lower.Conf.Limit.R2}}
  \item{Upper.Conf.Limit.R2 }{upper limit of the confidence interval around the population multiple correlation coefficient}
  \item{Prob.Greater.Upper }{proportion of the distribution greater than \code{Lower.Conf.Limit.R2}}
}
\references{
Algina, J. \& Olejnik, S. (2000) Determining Sample Size for Accurate Estimation of 
the Squared Multiple Correlation Coefficient. \emph{Multivariate Behavioral Research, 35}, 
119--136.

Steiger, J. H. \& Fouladi, R. T. (1992) R2: {A} computer program for interval estimation, power calculation,
and hypothesis testing for the squared multiple correlation. \emph{Behavior research methods, instruments and computers, 4}, 581--582.

}
\author{Ken Kelley (Indiana University; \email{KKIII@Indiana.Edu})}

\seealso{'ss.aipe.R2', 'conf.limits.nct'}
\examples{
ci.R2(R2=.25, N=100, p=5, conf.level=.95)
ci.R2(F.value=6.266667, N=100, p=5, conf.level=.95)

# Using conversion functions
ci.R2(F.value=Rsquare2F(R2=.25, p = 5, N = 100), N=100, p=5, 
conf.level=.95)
ci.R2(R2=F2Rsquare(F.value=6.266667, df.1 = 5, df.2 = 94), N=100, p=5,
conf.level=.95)

# One sided confidence intervals.
ci.R2(R2=.25, N=100, p=5, alpha.lower=.05, alpha.upper=0, conf.level=NULL)
ci.R2(R2=.25, N=100, p=5, alpha.lower=0, alpha.upper=.05, conf.level=NULL)
}
\keyword{multivariate}
\keyword{htest}
\keyword{regression}
