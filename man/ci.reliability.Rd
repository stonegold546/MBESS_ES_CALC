\name{ci.reliability}
\alias{ci.reliability}
\title{ Confidence Interval for a Reliability Coefficient }
\description{
A function to calculate the confidence interval for a reliability coefficient using 
the factor analytic approach or the normal theory approach. 
}
\usage{
ci.reliability(S = NULL, data = NULL, N = NULL, 
model = "True-Score Equivalent", type = "Factor Analytic", conf.level = 0.95, 
interval = TRUE, Bootstrap = FALSE, B = 10000, BootstrapCI = "BCa")
}

\arguments{
  \item{S}{symmetric covariance or correlation matrix }
  \item{data}{the data set that the reliability coefficient is obtained from}
  \item{N}{ the total sample size }
  \item{model}{ the type of measurement model (e.g., \code{"parallel items"}, \code{"true-score equivalent"}, or 
    \code{"congeneric model"}) for a homogeneous single common factor test }
  \item{type}{the type of method to base the formation of the confidence interval on, either the \code{"Factor Analytic"} (McDonald, 1999) or \code{"Normal Theory"} (van Zyl, Neudecker, & Nel, 2000) }
  \item{conf.level}{ the confidence level (i.e., 1-Type I error rate)}
  \item{interval}{whether or not to compute a confidence interval. If \code{FALSE}, then only a point estimate of the reliability will be returned.}
  \item{Bootstrap}{ whether or not to use the bootstrap to form a confidence interval. The default is to use analytic methods.}
  \item{B}{the number of bootstrap replications}
  \item{BootstrapCI}{the type of bootstrap confidence interval. It can be \code{BCa} or \code{percentile}.}
}
\details{
This function calculates a reliability coefficient for a set of scores based on the measurement model chosen. 
Thus, coefficient alpha is calculated for parallel items and true-score equivalent; and coefficient omega 
is calculated for congeneric items. See McDonald (1999) for the assumptions of each of these models. 
Under the Normal Theory method the asymptotic distribution of the maximum likelihood estimator of coefficient 
alpha is used for the true-score equivalent model (van Zyl et al., 2000). Note that this model is not 
optimal for small samples (Yuan & Bentler, 2002), for instance with N<20. 
}

\value{
  \item{ci.lower }{the lower bound of the computed confidence interval}
  \item{ci.upper}{the upper bound of the computed confidence interval}
  \item{Estimated.reliability}{the estimated reliability coefficient}
  \item{SE.reliability}{the standard error of the reliability coefficient}
  \item{Conf.Level}{the confidence level (i.e., 1 - Type I error rate)}
}

\references{ 
McDonald, R. P. (1999). \emph{Test theory: A unified approach}. Mahwah, New Jersey: Lawrence 
Erlbaum Associates, Publishers.

van Zyl, J. M., Neudecker, H., & Nel, D. G. (2000) On the distribution of the maximum likelihood estimator 
of Cronbach's alpha. \emph{Psychometrika, 65} (3), 271-280.

Yuan, K. & Bentler, P. M. (2002) On robustness of the normal-theory based asymptotic distributions of three
reliability coefficient estimates. \emph{Psychometrika, 67} (2), 251-259.
}

\author{Keke Lai (University of Notre Dame); Leann J. Terry (Indiana University; \email{ljterry@Indiana.Edu});
Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}); Sunthud Pornprasertmanit (University of Kansas; \email{PSunthud@KU.Edu})}


\seealso{ \code{\link{CFA.1}}; \code{\link[sem]{sem}} }
\examples{
# library(sem)

## Forming a hypothetical population covariance matrix
# Pop.Cov.Mat <- matrix(.3, 9, 9)
# diag(Pop.Cov.Mat) <- 1

# ci.reliability(S=Pop.Cov.Mat, N=50, model="True-Score", type="Normal Theory")
}

\keyword{ htest }
