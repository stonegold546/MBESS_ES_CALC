\name{CFA.1}
\alias{CFA.1}

\title{ One-factor confirmatory factor analysis model}

\description{Returns the MLE estimates and the estimated asymptotic covariance matrix of parameter estimates for one-factor confirmatory factor analysis model
}

\usage{
CFA.1(S, N, equal.loading = FALSE, equal.error = FALSE)
}

\arguments{
  \item{S}{ covariance matrix of the indicators }
  \item{N}{ total sample size }
  \item{equal.loading}{logical statement indicating whether the path coefficients are the same}
  \item{equal.error}{ logical statement indicating whether the manifest variables have the same error variances}
}

\value{
  \item{Model }{the factor analysis model specified by the user}
  \item{Factor.Loadings}{factor loadings}
  \item{Indicator.var}{the error variances of the indicator variables}
  \item{Parameter.cov}{the covariance matrix of the parameters}
}
%\references{}

\author{ Keke Lai (University of Notre Dame, \email{Lai.15@ND.Edu}) }

\seealso{ \code{\link[sem]{sem}}, \code{covmat.from.cfm} }
\examples{
# Construct covariance matrix
library(sem)
cov.mat<-read.moments()
1.384
1.484 2.756
1.988 2.874 4.845
2.429 3.588 4.894 6.951
3.031 4.390 6.080 7.476 10.313

tran.cov.mat <- t(cov.mat)
cov.mat[upper.tri(cov.mat)] <- tran.cov.mat[upper.tri(tran.cov.mat)]
CFA.1(N=300, S=cov.mat)
}

\keyword{ multivariate}
