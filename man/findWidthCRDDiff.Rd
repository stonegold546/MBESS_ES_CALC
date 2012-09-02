\name{findWidthCRDDiff}
\alias{findWidthCRDDiff}
\title{
	Find the width of CI of unstandardized condition means difference
}
\description{
	Find the width of CI of unstandardized condition means difference. Users may use the model with or without covariates. See further details at Pornprasertmanit and Schneider (2010, submitted).
}
\usage{
findWidthCRDDiff(nclus, nindiv, prtreat, tauy=NULL, sigma2y=NULL, 
	totalvar=NULL, iccy=NULL, r2between = 0, r2within = 0, numpredictor = 0, 
	assurance=NULL, conflevel = 0.95)
}
\arguments{
  \item{nclus}{
	The desired number of clusters
  }
  \item{nindiv}{
	The number of individuals in each cluster (cluster size)
  }
  \item{prtreat}{
	The proportion of treatment clusters
  }
  \item{tauy}{
	The residual variance in the between level before accounting for the covariates
  }
  \item{sigma2y}{
	The residual varaince in the within level before accounting for the covariate
  }
  \item{totalvar}{
	The total resiudal variance before accounting for the covariate
  }
  \item{iccy}{
	The intraclass correlation of the dependent variable
  }
  \item{r2within}{
	The proportion of variance explained in the within level (used when \code{covariate = TRUE})
  }
  \item{r2between}{
	The proportion of variance explained in the between level (used when \code{covariate = TRUE})
  }
  \item{numpredictor}{
	The number of predictors used in the between level
  }
  \item{assurance}{
	The degree of assurance, which is the value with which confidence can be placed that describes the likelihood of obtaining a confidence interval less than the value specified (e.g, .80, .90, .95)
  }
  \item{conflevel}{
	The desired level of confidence for the confidence interval
  }
}
\value{
	The width of CI of unstandardized condition means difference. If \code{assurance = NULL}, the value represents the expected width. If \code{assurance} is specified as a number, the width value will have the proprotion of the specified assurance that the the likelihood of obtaining a confidence interval less than the outcome width 
}
\references{ 
Pornprasertmanit, S., & Schneider, W. J. (2010). \emph{Efficient sample size for power and desired accuracy in Cohen's d estimation in two-group cluster randomized design} (Master Thesis). Illinois State University, Normal, IL. 

Pornprasertmanit, S., & Schneider, W. J. (submitted). \emph{Accuracy in parameter estimation in two-condition cluster randomized design}.
}
\author{
    Sunthud Pornprasertmanit (University of Kansas; \email{psunthud@ku.edu})
}
\examples{
\dontrun{
findWidthCRDDiff(80, 10, 0.25, tauy=0.25, sigma2y=0.75)
}
}
