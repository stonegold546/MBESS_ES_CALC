\name{ci.reliability}
\alias{ci.reliability}
\title{ Confidence Interval for a Reliability Coefficient }
\description{
A function to calculate the confidence interval for a reliability coefficient: coefficient alpha or coefficient omega. 
}
\usage{
ci.reliability(data = NULL, S = NULL, N = NULL, type = "omega", 
    analysis.type = "default", interval.type = "normal-theory", B = 1000, conf.level = 0.95)
}
\arguments{
  \item{data}{The dataset that the reliability coefficient is obtained from.}
  \item{S}{Symmetric covariance or correlation matrix }
  \item{N}{The total sample size }
  \item{type}{The type of reliability coefficient to be calculated: \code{"alpha"} for coefficient alpha or \code{"omega"} for coefficient omega}
  \item{analysis.type}{the type of method used to calculate reliability coefficient: \code{"analytic"} which uses formula or \code{"cfa"} which uses confirmatory factor analysis to calculate reliability coefficient. The default of this argument for coefficient alpha is the analytic approach. The default of this argument for coefficient omega is the confirmatory factor analytic approach.}
  \item{interval.type}{There are 13 options for the methods. See \code{details} below. The default is to use \code{"wald"} approach for both reliability coefficients.}
  \item{B}{the number of bootstrap replications}
  \item{conf.level}{ the confidence level (i.e., 1-Type I error rate)}
}
\details{
When the coefficient alpha is used, the measurement model is assumed to be true-score equivalent (or tau equivalent) model such that factor loadings are equal across items. When the coefficient omega is used, the measurement model is assumed to be congeneric model (similar to one-factor confirmatory factor analysis model). If researchers wish to make the measurement model with all parallel items (equal factor loadings and equal error variances), users can specify it by setting \code{interval.type="parallel"} and \code{type="alpha"}. See McDonald (1999) for the assumptions of each of these models. Note that the \code{"wald"} approach to estimate confidence interval of a reliability coefficient is not optimal for small samples (Yuan & Bentler, 2002), for instance with \emph{N} < 20. 

The list below shows all methods to find the confidence interval of reliability. 

\enumerate{
	\item \code{"none"} to not find any confidence interval
	\item \code{"fisher"} for the Fisher's \eqn{z} transformation on the correlation coefficient approach, \eqn{z = 0.5 \times \log{\frac{1 + \alpha}{1 - \alpha}}}, directly on the coefficient alpha and find confidence interval of transformed scale (Fisher, 1950). The variance of the \eqn{z} is \eqn{\frac{1}{N - 3}} where \eqn{N} is the total sample size.
	\item \code{"icc"} for the Fisher's \eqn{z} transformation on the intraclass correlation approach, \eqn{z = \log{1 - \alpha}}. The variance of the \eqn{z} is \eqn{\frac{2p}{N(p - 1)}} where \eqn{p} is the number of items (Fisher, 1991, p. 221; van Zyl et al., 2000, p. 277).
	\item \code{"bonett"} for the Fisher's \eqn{z} transformation on the intraclass correlation approach with the variance of \eqn{\frac{2p}{(N - 2)(p - 1)}} (Bonett, 2002, Equation 6).
	\item \code{"feldt"} is based on that \eqn{\frac{1 - \alpha}{1 - \hat{alpha}}} is distributed as \eqn{F} distribution with the degree of freedoms of \eqn{N - 1} and \eqn{(N - 1) \times (p - 1)} (Feldt, 1965).
	\item \code{"siotani"} is the same as the \code{"feldt"} method but using the degree of freedoms of \eqn{N} and \eqn{N \times (p - 1)} (Siotani, Hayakawa, & Fujikoshi, 1985; van Zyl et al., 2000, Equations 7 and 8; also referred as the exact method of Koning & Franses, 2003).
	\item \code{"parallel"} to assume that the items are parallel and analyze confidence interval based on normal-theory approach (see van Zyl, Neudecker, & Nel, 2000, Equation 22; also referred as the asymptotic method of Koning & Franses, 2003).
	\item \code{"hakstian"} uses the cube root transformation and assumes normal distribution on the cube root transformation (Hakstian & Whalen, 1976). The variance of the transformed reliability is based on the degrees of freedom in the \code{"feldt"} method.
	\item \code{"hakstianbarchard"} uses a correction of the violation of compound symmetry of covariance matrix by adjusting the degrees of freedom in the \code{"hakstian"}. This correction is used for the inference in type 12 sampling (both persons and items are sampled from the population of persons and items) See Hakstian and Barchard (2000) for further details.
	\item \code{"wald"} to analyze the confidence interval based on normal-theory approach (or multivariate delta method). See van Zyl, Neudecker, & Nel (2000, Equation 21) for the confidence interval of coefficient alpha (also be referred as Iacobucci & Duhachek's, 2003, method). See Raykov (2002) for details for coefficient omega. If users use \code{analytic.type="cfa"}, the \code{sem} package will be used to obtain parameter estimates and standard errors used for the formula proposed by Raykov (2002)
	\item \code{"adf"} for asymtotic distribution-free method (see Maydeu-Olivares, Coffman, & Hartman, 2007 for further details for coefficient omega; we use phantom variable approach, Cheung, 2009, and \code{"WLS"} estimator for coefficient omega, Browne, 1984, in the \code{lavaan} package, Rosseel, in press).
	\item \code{"ll"} for profile likelihood-based confidence interval of both reliability coefficients (Cheung, 2009) analyzed by the \code{OpenMx} package (Boker et al., 2011)
	\item \code{"bsi"} for standard bootstrap confidence interval which find the standard deviation across the bootstrap estimates, multiply the standard deviation by critical value, and add and subtract from the reliability estimate.
	\item \code{"perc"} for percentile bootstrap confidence interval.
	\item \code{"bca"} for bias-corrected and accelerated bootstrap confidence interval.
}
All method can be used with \code{type="alpha"}. Only \code{"wald"}, \code{"adf"}, \code{"ll"}, \code{"perc"}, and \code{"bca"} can be used for \code{type="omega"}. The \code{"fisher"}, \code{"icc"}, \code{"bonett"}, \code{"feldt"}, \code{"siotani"}, and \code{"parallel"} assumes that the items are paralleled, which means the true scores and error variances are equal across items.
}
\value{
  \item{est}{the estimated reliability coefficient}
  \item{se}{the standard error of the reliability coefficient}
  \item{ci.lower }{the lower bound of the computed confidence interval}
  \item{ci.upper}{the upper bound of the computed confidence interval}
  \item{conf.Level}{the confidence level (i.e., 1 - Type I error rate)}
  \item{est.type}{the type of estimated reliability coefficient (alpha or omega)}
  \item{analysis.type}{the type of method used to estimate reliability coefficient (analytic or cfa)}
  \item{interval.type}{the method used to find confidence interval}
  \item{call}{the input the users used to call this function}
}
\references{ 
Boker, S., M., N., Maes, H., Wilde, M., Spiegel, M., Brick, T., et al. (2011).
OpenMx: An open source extended structural equation modeling framework.
\emph{Psychometrika, 76,} 306-317.

Bonett, D. G. (2002). Sample size requirements for testing and estimating coefficient
alpha. \emph{Journal of Educational and Behavioral Statistics, 27,} 335-340.

Browne, M. W. (1984). Asymptotic distribution free methods in the analysis of
covariance structures. \emph{British Journal of Mathematical and Statistical
Psychology, 24,} 445-455.

Cheung, M. W.-L. (2009). Constructing approximate confidence intervals for
parameters with structural constructing approximate confidence intervals for
parameters with structural equation models. \emph{Structural Equation Modeling,
16,} 267-294.

Feldt, L.S. (1965). The approximate sampling distribution of Kuder-Richardson reliability coefficient twenty. \emph{Psychometrika,
30}, 357-370.

Fisher, R. A. (1950). \emph{Statistical methods for research workers.} Edinburgh, UK:
Oliver & Boyd.

Fisher, R.A. (1991). Statistical methods for research workers. In J.H. Bennett (Ed.), \emph{Statistical methods, experimental
design, and scientific inference.} Oxford: Oxford University Press.

Hakstian, A. R., & Whalen, T. E. (1976). A k-sample significance test for independent alpha
coefficients. \emph{Psychometrika, 41,} 219-231.

Iacobucci, D., & Duhachek, A. (2003). Advancing alpha: measuring reliability with confidence.
\emph{Journal of Consumer Psychology, 13,} 478-487.

Koning, A. J., & Franses, P. H. (2003). \emph{Confidence intervals for Cronbach&grave;s coefficient
alpha values} (ERIM Report Series Ref. No. ERS-2003-041-MKT). Rotterdam, The Netherlands:
Erasmus Research Institute of Management.

Maydeu-Olivares, A., Coffman, D. L., & Hartmann, W. M. (2007). Asymptotically
distribution-free (ADF) interval estimation of coefficient alpha. \emph{Psychological
Methods, 12,} 157-176.

McDonald, R. P. (1999). \emph{Test theory: A unified approach}. Mahwah, New Jersey: Lawrence 
Erlbaum Associates, Publishers.

Raykov, T. (2002). Analytic estimation of standard error and confidence interval
for scale reliability. \emph{Multivariate Behavioral Research, 37,} 89-103.

Rosseel, Y. (in press). lavaan: An R package for structural equation modeling. \emph{Journal of Statistical Software}.

Siotani, M., Hayakawa, T., & Fujikoshi, Y. (1985). \emph{Modem multivariate statistical analysis: A graduate course and
handbook}. Columbus, Ohio: American Sciences Press.

van Zyl, J. M., Neudecker, H., & Nel, D. G. (2000) On the distribution of the maximum likelihood estimator 
of Cronbach's alpha. \emph{Psychometrika, 65} (3), 271-280.

Yuan, K. & Bentler, P. M. (2002) On robustness of the normal-theory based asymptotic distributions of three
reliability coefficient estimates. \emph{Psychometrika, 67} (2), 251-259.
}
\author{Sunthud Pornprasertmanit (University of Kansas; \email{PSunthud@KU.Edu})}
\seealso{ \code{\link{CFA.1}}; \code{\link[sem]{sem}}; \code{\link[lavaan]{lavaan}} }
\examples{
# Use this function for the attitude dataset (ignoring the overall rating variable)
ci.reliability(data=attitude[,-1])

## Forming a hypothetical population covariance matrix
Pop.Cov.Mat <- matrix(.3, 9, 9)
diag(Pop.Cov.Mat) <- 1

ci.reliability(S=Pop.Cov.Mat, N=50, type="alpha")
}
\keyword{ htest }
