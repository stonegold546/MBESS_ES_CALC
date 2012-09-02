\name{ci.reliability.bs}
\alias{ci.reliability.bs}

\title{ Bootstrap the confidence interval for reliability coefficient }
\description{
Obtain the confidence interval for reliability coefficient with bootstrap
}
\usage{
ci.reliability.bs(data, type = "Congeneric", analysis.type = "Factor Analytic", conf.level = 0.95, B = 100, ...)
}

\arguments{
  \item{data}{ the data set that the reliability coefficient is obtained from}
  \item{type}{ the type of measurement model (e.g., \code{"parallel items"}, \code{"true-score equivalent"}, or 
    \code{"congeneric model"}) for a homogeneous single common factor test }
  \item{analysis.type}{the type of method to base the formation of the confidence interval on, either the \code{"Factor Analytic"} (McDonald, 1999) or \code{"Normal Theory"} (van Zyl, Neudecker, & Nel, 2000) }
  \item{conf.level}{ the confidence level (i.e., 1-Type I error rate)}
  \item{B}{ The number of bootstrap replicates  }
  \item{\dots}{ Any other arguments for the \code{boot} function}
}
%\details{}

\value{
  \item{Percentile.Method }{Confidence interval obtained with the percentile bootstrap method}
  \item{BCa.Method }{confidence interval obtained with the bias-correct and accelerated bootstrap method}
}

\references{ Efron, B. & Tibshirani, R. J. (1993). \emph{An introduction to the bootstrap}. New York: Chapman & Hall }
\author{ Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}); Keke Lai (University of Notre Dame) }

%\note{ }

\seealso{ \code{\link{ci.reliability}}, \code{\link{CFA.1}}, \code{\link[boot]{boot}} }

%\examples{}

\keyword{ htest }
