\name{ss.aipe.R2}
\alias{ss.aipe.R2}
\title{Sample Size Planning for Accuracy in Parameter Estimation (i.e., precision) 
for the multiple correlation coefficient.}
\description{Determines necessary sample size for the multiple correlation coefficient so that the 
confidence interval around the population multiple correlation coefficient is 
sufficiently narrow. Optionally, there is a certainty parameter that allows one to be 
a specified percent certain that the observed interval will be no wider than specified.}

\usage{ss.aipe.R2(Population.R2 = NULL, conf.level = .95, width = NULL, 
Random.Predictors=TRUE, Random.Regressors, which.width = "Full", p = NULL, 
degree.of.certainty = NULL, verify.ss=FALSE, Tol = 1e-09, ...)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Population.R2}{value of the population multiple correlation coefficient}
  \item{conf.level}{confidence interval level (e.g., .95, .99, .90); 1-Type I error rate}
  \item{width}{width of the confidence interval (see \code{which.width})}
  \item{Random.Predictors}{whether or not the predictor variables are random (set to TRUE) or are fixed (set to FALSE)}
  \item{Random.Regressors}{an alias for \code{Random.Predictors}; \code{Random.Regressors} overrides \code{Random.Predictors}}
  \item{which.width}{defines the width that \code{width} refers to}
  \item{p}{the number of predictor variables}
  \item{degree.of.certainty}{value with which confidence can be placed that describes the likelihood of obtaining a confidence interval less than the value specified (e.g., .80, .90, .95)}
  \item{verify.ss}{evaluates numerically via an internal Monte Carlo simualtion the exact sample size given the specifications}
  \item{Tol}{the tolerance of the iterative function \code{conf.limits.nct} for convergence}
  \item{\dots}{for modifying parameters of functions this function calls }
}

\details{
This function determines a necessary sample size so that the expected confidence interval width for
the squared multiple correlation coefficient is sufficiently narrow (when \code{degree.of.certainty=NULL})
or so that the obtained confidence interval is no larger than the value specified with some desired degree 
of certainty (i.e., a probability that the obtained width is less than the specified width). The 
method depends on whether or not the regressors are regarded as fixed or random. This is the case because the distribution theory for the two cases is different and thus 
the confidence interval procedure is conditional the type of regressors. Kelley (2006) and 
Kelley \& Maxwell (In press) detail the methods used in the function, with the former focusing on
random regressors and the latter on fixed regressors.

It is recommended that the option \code{verify.ss} should always be used! Doing so uses the method impled sample size as an estimate
and then evaluates with an internal Monte Carlo simulation (i.e., via "brute-force" methods) the exact 
sample size given the goals specified. When \code{verify.ss=TRUE}, the default number of iterations is 10,000
but this can be changed by specifying G=5000 (or some other value; 10000 is the recommended)
When \code{verify.ss=TRUE} is specified, an internal function \code{verify.ss.aipe.R2} calls upon the 
\code{ss.aipe.R2.sensitivity} function for purposes of the internal Monte Carlo simulation study. See 
the \code{verify.ss.aipe.R2} function for arguments that can be passed from \code{ss.aipe.R2} to 
\code{verify.ss.aipe.R2}.
}


\value{
\item{Required.Sample.Size}{sample size that should be used given the conditions specified.}
%\item{Expected.Width}{given the sample size, the expected (half)-width of interest.}
}

\references{
Algina, J. \& Olejnik, S. (2000) Determining sample size for accurate estimation of the squared
multiple correlation coefficient. \emph{Multivariate Behavioral Research, 35}, 119--136.

Steiger, J. H. \& Fouladi, R. T. (1992) R2: {A} computer program for interval estimation, power calculation,
and hypothesis testing for the squared multiple correlation. \emph{Behavior research methods, instruments and computers, 4}, 581--582.

Kelley, K. Sample size planning for the squared multiple correlation coefficient:
Accuracy in parameter estimation via narrow confidence intervals, \emph{manuscripted submitted for publication}.

Kelley, K. \& Maxwell, S. E. (In press) Power and accuracy for omnibus and targeted effects: 
Issues of sample size planning with applications to multiple regression. In P. Alasuuta, J. Brannen, \& L. Bickman (Eds.),
\textit{Handbook of Social Research Methods}. Newbury Park, CA: Sage.

}

\author{Ken Kelley (Indiana University; \email{KKIII@Indiana.Edu})}
\note{This function can be slow to converge (e.g., 1 minute convergence in some situations). This is 
because the function is written so that it handles (essentially) all values and is thus very general.}

\section{Warning}{Only \code{Full} in \code{which.width} should be used at the present time. Sample size returned for \code{Lower} and \code{Upper} widths are only approxiamate. If these widths are
of interest, try using the \code{ss.aipe.R2.sensitivity} function to determine sample size through 
brute force (trial and error) procedures.}

\seealso{'ci.R2', 'conf.limits.nct', 'ss.aipe.R2.sensitivity'}

\examples{
ss.aipe.R2(Population.R2=.50, conf.level=.95, width=.10, which.width="Full",
p=5, Random.Predictors=TRUE)

# Same as above, except the predictor variables are considered fixed.
ss.aipe.R2(Population.R2=.50, conf.level=.95, width=.10, which.width="Full",
p=5, Random.Predictors=FALSE)

ss.aipe.R2(Population.R2=.50, conf.level=.95, width=.10, which.width="Full",
p=5, degree.of.certainty=.85, Random.Predictors=TRUE)

ss.aipe.R2(Population.R2=.50, conf.level=.95, width=.10, which.width="Full",
p=5, degree.of.certainty=.85, Random.Predictors=FALSE)
}
\keyword{design}
