\name{ss.aipe.R2}
\alias{ss.aipe.R2}
\title{Sample Size Planning for Accuracy in Parameter Estimation (i.e., precision) 
for the multiple correlation coefficient.}
\description{Determines necessary sample size for the multiple correlation coefficient so that the 
confidence interval around the population multiple correlation coefficient is 
sufficiently narrow. Optionally, there is a certainty parameter that allows one to be 
a specified percent certain that the observed interval will be no wider than specified.}
\usage{ss.aipe.R2(Population.R2 = NULL, conf.level = .95, width = NULL, 
Random.Predictors=TRUE, which.width = "Full", p = NULL, 
degree.of.certainty = NULL, 
Tol = 1e-09, ...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Population.R2}{value of the population multiple correlation coefficient}
  \item{conf.level}{confidence interval level (e.g., .95, .99, .90); 1-Type I error rate}
  \item{width}{width of the confidence interval (see \code{which.width})}
  \item{Random.Predictors}{whether or not the predictor variables are random (set to TRUE) or are fixed (set to FALSE)}
  \item{which.width}{defines the width that \code{width} refers to}
  \item{p}{the number of predictor variables}
  \item{degree.of.certainty}{value with which confidence can be placed that describes the likelihood of obtaining a confidence interval less than the value specified (e.g., .80, .90, .95)}
  \item{Tol}{the tolerance of the iterative function \code{conf.limits.nct} for convergence}
  \item{\dots}{for modifying parameters of functions this function calls }
}

%\details
% {
%  ~~ If necessary, more details than the description above ~~
% }
\value{\item{Required.Sample.Size}{sample size that should be used given conditions specified.}
  \item{Expected.Width}{given the sample size, the expected (half)-width of interest.}
  }

\references{
Algina, J. \& Olejnik, S. (2000) Determining Sample Size for Accurate Estimation of the Squared
Multiple Correlation Coefficient. \emph{Multivariate Behavioral Research, 35}, 119--136.

Steiger, J. H. \& Fouladi, R. T. (1992) R2: {A} computer program for interval estimation, power calculation,
  and hypothesis testing for the squared multiple correlation. \emph{Behavior research methods, instruments and computers, 4}, 581--582.

Kelley, K. Sample size planning for the squared multiple correlation coefficient:
Accuracy in parameter estimation via narrow confidence intervals, \emph{manuscripted submitted for publication}.

}

\author{Ken Kelley (Indiana University; \email{KKIII@Indiana.Edu})}
\note{This function can be slow to converge (e.g., 1 minute convergence is not uncommon). This is 
because the function is written so that it handles (essentially) all values and is thus very general.}

\section{Warning}{Only \code{Full} in \code{which.width} should be used at the present time. Sample size returned for \code{Lower} and \code{Upper} widths are only approxiamate. If these widths are
of interest, try using the \code{ss.aipe.R2.sensitivity} function.}


\seealso{'ci.R2', 'conf.limits.nct', 'ss.aipe.R2.sensitivity'}

\examples{
ss.aipe.R2(Population.R2=.50, conf.level=.95, width=.10, which.width="Full",
p=5, Random.Predictors=TRUE)

# Same as above, except the predictor variables are considered fixed.
# Fixed predictors do not occur often in applied behavioral, educational,
# and social science research. 
ss.aipe.R2(Population.R2=.50, conf.level=.95, width=.10, which.width="Full",
p=5, Random.Predictors=FALSE)

ss.aipe.R2(Population.R2=.50, conf.level=.95, width=.10, which.width="Full",
p=5, degree.of.certainty=.85, Random.Predictors=TRUE)

ss.aipe.R2(Population.R2=.50, conf.level=.95, width=.10, which.width="Full",
p=5, degree.of.certainty=.85, Random.Predictors=FALSE)
}
\keyword{design}
