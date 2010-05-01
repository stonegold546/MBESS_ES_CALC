\name{ss.aipe.rmsea.sensitivity}
\Rdversion{1.1}
\alias{ss.aipe.rmsea.sensitivity}
\title{a priori Monte Carlo simulation for sample size planning for RMSEA in SEM
}
\description{Conduct a priori Monte Carlo simulation to empirically study the effects of (mis)specifications of input information on the calculated sample size. The sample size is planned so that the expected width of a confidence interval for the population RMSEA is no larger than desired. Random data are generated from the true covariance matrix but fit to the proposed model, whereas sample size is calculated based on the input covariance matrix and proposed model.
}

\usage{
ss.aipe.rmsea.sensitivity(width, model.est, Sigma.est, 
Sigma.true = Sigma.est, N = NULL, conf.level = 0.95, G = 1000, ...)
}

\arguments{
  \item{width}{desired confidence interval width for the model parameter of interest }
  \item{model.est}{ the model the researcher proposes, may or may not be the true model. This argument should be an RAM (reticular action model; e.g., McArdle & McDonald, 1984) specification of a structural equation model, and should be of class \code{mod}. The model is specified in the same manner as does the \code{\link[sem]{sem}} package; see \code{\link[sem]{sem}} and \code{\link[sem]{specify.model}} for detailed documentation about model specifications in the RAM notation.}

  \item{Sigma.est}{ the covariance matrix used to calculate sample size, may or may not be the true covariance matrix. The row names and column names of \code{Sigma.est} should be the same as the manifest variables in \code{model.est}.}
  
  \item{Sigma.true}{the true population covariance matrix, which will be used to generate random data for the simulation study. The row names and column names of \code{Sigma.true} should be the same as the manifest variables in \code{model.est}.
}

  \item{N}{if \code{N} is specified, random sample of the specified \code{N} size will be generated. Otherwise the sample size is calculated with the sample size planning method with the goal that the expected width of a confidence interval for population RMSEA is no larger than desired.
}
  \item{conf.level}{ confidence level (i.e., 1- Type I error rate)}
  \item{G}{number of replications in the Monte Carlo simulation }
  \item{\dots}{  allows one to potentially include parameter values for inner functions}

}
\details{This function implements the sample size planning methods proposed in Kelley and Lai (2010). It depends on the 
function \code{\link[sem]{sem}} in the \code{sem} package to fit the proposed model to random data, and uses the same notation to specify SEM 
models as does \code{\link[sem]{sem}}. Please refer to \code{\link[sem]{sem}} for more detailed documentation 
about model specifications, the RAM notation, and model fitting techniques. For technical discussion 
on how to obtain the model implied covariance matrix in the RAM notation given model parameters, see McArdle and McDonald (1984)
}
\value{
\item{w}{the \code{G} random confidence interval widths }
\item{RMSEA.hat}{the \code{G} estimated RMSEA values based on the \code{G} random samples }
\item{sample.size}{the sample size calculated }
\item{df}{degrees of freedom of the proposed model}
\item{RMSEA.pop}{the input RMSEA value that is used to calculated the necessary sample size}
\item{desired.width}{desired confidence interval width}
\item{mean.width}{mean of the \code{G} random confidence interval widths}
\item{median.width}{median of the \code{G} random confidence interval widths}
\item{assurance}{the proportion of confidence interval widths narrower than desired} 
\item{quantile.width}{99, 97, 95, 90, 80, 70, and 60 percentiles of the \code{G} random confidence interval widths}
\item{alpha.upper}{the upper empirical Type I error rate} 
\item{alpha.lower}{the lower empirical Type I error rate}
\item{alpha}{total empirical Type I error rate}  
\item{conf.level}{confidence level}
\item{suc.rep}{successful replications} 
}

\references{
Fox, J. (2006). Structural equation modeling with the sem package in R. \emph{Structural Equation Modeling, 13}, 465-486.

Kelley, K., & Lai, K. (2010). Accuracy in parameter estimation for the root mean square of approximation: Sample size planning for narrow confidence intervals. \emph{Manuscript under review}.

McArdle, J. J., & McDonald, R. P. (1984). Some algebraic properties of the reticular action model. \emph{British Journal of Mathematical and Statistical Psychology, 37}, 234-251.
}
\author{
Keke Lai (University of Notre Dame; \email{KLai1@ND.Edu})
}
\note{Sometimes the simulation stops in the middle of fitting the model to the random data. The reason is because \code{\link{nlm}}, the 
function \code{\link[sem]{sem}} calls to fit the model, fails to converge. We suggest using the \code{\link[base]{try}} function in simulation so that
the simulation can proceed with unsuccessful iterations.
}


\seealso{
\code{\link[sem]{sem}}; \code{\link[sem]{specify.model}}; \code{\link{ss.aipe.rmsea}}; \code{\link{theta.2.Sigma.theta}}
}
%\examples{}

\keyword{ design }
\keyword{ multivariate }