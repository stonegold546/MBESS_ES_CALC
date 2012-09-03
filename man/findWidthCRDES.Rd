\name{findWidthCRDES}
\alias{findWidthCRDES}
\title{
	Find the width of CI of standardized condition means difference
}
\description{
	Find the width of CI of standardized condition means difference. Users may use the model with one covariate or without any covariates. See further details at Pornprasertmanit and Schneider (2010, submitted).
}
\usage{
findWidthCRDES(nrep, nclus, ntreatclus, nindiv, iccy, es, estype = 1, 
	totalvar=1, covariate=FALSE, iccz=NULL, r2within=NULL, r2between=NULL, 
	totalvarz = 1, assurance=NULL, seed=123321, multicore=FALSE, numProc=NULL,
	conflevel = 0.95)
}
\arguments{
  \item{nrep}{
	The number of replications used in a priori Monte Carlo simulation
  }
  \item{nclus}{
	The total number of clusters
  }
  \item{ntreatclus}{
	The number of treatment clusters
  }
  \item{nindiv}{
	The number of individuals in each cluster (cluster size)
  }
  \item{iccy}{
	The intraclass correlation of the dependent variable
  }
  \item{es}{
	The amount of effect size
  }
  \item{estype}{
	The type of effect size. There are only three possible options: 0 = the effect size using total standard deviation, 1 = the effect size using the individual-level standard deviation (level 1), 2 = the effect size using the cluster-level standard deviation (level 2)
  }
  \item{totalvar}{
	The amount of total variance in the dependent variable
  }
  \item{covariate}{
	If \code{TRUE}, a single covariate is included into the model. If \code{FALSE}, the no-covariate model is used.
  }
  \item{iccz}{
	The intraclass correlation of the covariate (used when \code{covariate = TRUE}). If \code{iccz = 0}, the within-level covariate will be only used. If \code{iccz = 1}, the between-level covariate will be only used.
  }
  \item{r2within}{
	The proportion of variance explained in the within level (used when \code{covariate = TRUE})
  }
  \item{r2between}{
	The proportion of variance explained in the between level (used when \code{covariate = TRUE})
  }
  \item{totalvarz}{
	The total variance of the covariate (used when \code{covariate = TRUE})
  }
  \item{assurance}{
	The degree of assurance, which is the value with which confidence can be placed that describes the likelihood of obtaining a confidence interval less than the value specified (e.g, .80, .90, .95)
  }
  \item{seed}{
	A desired seed number
  }
  \item{multicore}{
	Use multiple processors within a computer. Specify as \code{TRUE} to use it. 
  }
  \item{numProc}{
	The number of processors to be used when \code{multicore=TRUE}. If it is not specified, the package will use the maximum number of processors in a machine. 
  }
  \item{conflevel}{
	The desired level of confidence for the confidence interval
  }
}
\value{
	The width of CI of standardized condition means difference. If \code{assurance = NULL}, the value represents the expected width. If \code{assurance} is specified as a number, the width value will have the proprotion of the specified assurance that the the likelihood of obtaining a confidence interval less than the outcome width 
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
# The number of replications should be much higher than 20.
findWidthCRDES(20, 70, 35, 20, 0.05, 0.2, estype = 1)
}
}
