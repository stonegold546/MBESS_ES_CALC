\name{ss.aipe.R2.sensitivity}
\alias{ss.aipe.R2.sensitivity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sensitivity analysis for sample size planning with the goal of Accuracy in Parameter Estimation (i.e., a narrow observed confidence interval)}

\description{
Given \code{Estimated.R2} and \code{True.R2}, one perform a sensitivity analysis to determine the effect of a misspecified population 
squared multiple correlation coefficient using the Accuracy in Parameter Estimation (AIPE) approach to sample size planning. The function
evaluates the effect of a misspecified \code{True.R2} on the width of obtained confidence intervals. 
}
\usage{
ss.aipe.R2.sensitivity(True.R2 = NULL, Estimated.R2 = NULL, w = NULL, 
p = NULL, Selected.N=NULL, degree.of.certainty = NULL, 
conf.level = 0.95, rho.yx = 0.3, rho.xx = 0.3,  G = 10000, 
print.iter = TRUE, ...)
}
\arguments{
  \item{True.R2}{value of the population squared multiple correlation coefficient}
  \item{Estimated.R2}{value of the estimated (for sample size planning) squared multiple correlation coefficient}
  \item{w}{full confidence interval width of interest}
  \item{p}{number of predictors}
  \item{Selected.N}{selected sample size to use in order to determine distributional properties of at a given value of sample size}
  \item{degree.of.certainty}{parameter to ensure confidence interval width with a specified degree of certainty}
  \item{conf.level}{confidence interval coverage (symmetric coverage)}
  \item{rho.yx}{value of the correlation between \code{y} (dependent variable) and each of the \code{x} variables (independent variables)}
  \item{rho.xx}{value of the correlation among the \code{x} variables (independent variables)}
  \item{G}{number of generations (i.e., replications) of the simulation}
  \item{print.iter}{should the iteration number (between 1 and \code{G}) during the run of the function}
  \item{...}{for modifying parameters of functions this function calls upon}
}
\details{

When \code{Estimated.R2}=\code{True.R2}, the results are that of a simulation study when all assumptions
are satisfied. Rather than specifying \code{Estimated.R2}, one can specify \code{Selected.N} to determine the results of a particular sample size (when doing this \'Estimated.R2\' cannot be specified).

The sample size estimation procedure technically assumes multivariate normal variables (\code{p}+1) with fixed predictors (\code{x}/indepdent variables), 
yet the function assumes random multivariate normal predictors (having a \code{p}+1 multivariate distribution). As Gatsonis and Sampson (1989) note in the context of statistical 
power analysis (recall this function is used in the context of precision), there is little difference in the outcome.

}
\value{
  \item{Results}{a list containing vectors of the empirical results}
  \item{Specifications}{outputs the input specifications and required sample size}
  \item{Summary}{summary values for the results of the sensitivity analysis (simulation study)}
}
\references{
Gatsonis, C. \& Sampson, A. R. (1989). Multiple Correlation: Exact power and sample size calculations. \emph{Psychological Bulletin, 106}, 516--524.

Steiger, J. H. \& Fouladi, R. T. (1992). R2: A computer program for interval estimation, power calculation,
and hypothesis testing for the squared multiple correlation. \emph{Behavior research methods, instruments and computers, 4}, 581--582.

}
\author{Ken Kelley (Indiana University; \email{KKIII@Indiana.Edu})}

\seealso{'ci.R2', 'conf.limits.nct', 'ss.aipe.R2'}
\examples{
# Change 'G' to some large number (e.g., G=5,000)
ss.aipe.R2.sensitivity(True.R2=.5, Estimated.R2=.4, w=.10, p=5, conf.level=0.95,
G=100)
}
\keyword{design}
