\name{ss.aipe.reg.coef.sensitivity}
\alias{ss.aipe.reg.coef.sensitivity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sensitivity analysis for sample size planning from the Accuracy in Parameter Estimation Perspective for the (standardized and unstandardized) regression coefficient}
\description{
Performs a sensitivity analysis when planning sample size from the Accuracy in Parameter Estimation Perspective for the standardized or unstandardized regression coefficient.}
\usage{
ss.aipe.reg.coef.sensitivity(True.Var.Y = NULL, True.Cov.YX = NULL,
True.Cov.XX = NULL, Estimated.Var.Y = NULL, Estimated.Cov.YX = NULL,
Estimated.Cov.XX = NULL, Specified.N = NULL, which.predictor = 1, 
w = NULL, Noncentral = FALSE, Standardize = FALSE, conf.level = 0.95, 
degree.of.certainty = NULL, assurance=NULL, certainty=NULL,
 G = 1000, print.iter = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{True.Var.Y}{Population variance of the dependent variable (\emph{Y})}
  \item{True.Cov.YX}{Population covariances vector between the \code{p} predictor variables and the dependent variable (\emph{Y})}
  \item{True.Cov.XX}{Population covariance matrix of the \code{p} predictor variables}
  \item{Estimated.Var.Y}{Estimated variance of the dependent variable (\emph{Y})}
  \item{Estimated.Cov.YX}{Estimated covariances vector between the \code{p} predictor variables and the dependent variable (\code{Y})}
  \item{Estimated.Cov.XX}{Estimated Population covariance matrix of the \code{p} predictor variables}
  \item{Specified.N}{Directly specified sample size (instead of using \code{Estimated.Rho.YX} and \code{Estimated.RHO.XX})}
  \item{which.predictor}{identifies which of the \emph{p} predictors is of interest}
  \item{w}{desired confidence interval width for the regression coefficient of interest}
  \item{Noncentral}{specify with a \code{TRUE} or \code{FALSE} statement whether or not the noncentral approach to sample size planning should be used}
  \item{Standardize}{specify with a \code{TRUE} or \code{FALSE} statement whether or not the regression coefficient will be standardized}
  \item{conf.level}{desired level of confidence for the computed interval (i.e., 1 - the Type I error rate)}
  \item{degree.of.certainty}{degree of certainty that the obtained confidence interval will be sufficiently narrow}
  \item{assurance}{ an alias for \code{degree.of.certainty} }
\item{certainty}{an alias for \code{degree.of.certainty}}
  \item{G}{the number of generations/replication of the simulation student within the function}
  \item{print.iter}{specify with a \code{TRUE}/\code{FALSE} statement if the iteration number should be printed as the simulation within the function runts}
}
  
  \details{ Direct specification of \code{True.Rho.YX} and \code{True.RHO.XX} is necessary, even if one is interested in a single regression
coefficient, so that the covariance/correlation structure can be specified when the simulation student within the function runs.} 
\value{\item{Results}{a matrix containing the empirical results from each of the \code{G} replications of the simulation}
  \item{Specifications}{a list of the input specifications and the required sample size}
  \item{Summary.of.Results}{summary values for the results of the sensitivity analysis (simulation study) given the input specification}}

\references{Kelley, K. & Maxwell, S. E. (2003). Sample size for Multiple Regression: Obtaining regression coefficients that are accuracy, not simply significant. \emph{Psychological Methods, 8}, 305--321.}

\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}

\note{Note that when \code{True.Rho.YX}=\code{Estimated.Rho.YX} and \code{True.RHO.XX}=\code{Estimated.RHO.XX}, the results are not
literally from a sensitivity analysis, rather the function performs a standard simulation study. A simulation study 
can be helpful in order to determine if the sample size procedure under or overestimates necessary sample size.} 
\seealso{\code{ss.aipe.reg.coef}, \code{ci.reg.coef}}
\keyword{design}
