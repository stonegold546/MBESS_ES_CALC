\name{ss.aipe.crd}
\alias{ss.aipe.crd.nclus.fixedwidth}
\alias{ss.aipe.crd.nindiv.fixedwidth}
\alias{ss.aipe.crd.nclus.fixedbudget}
\alias{ss.aipe.crd.nindiv.fixedbudget}
\alias{ss.aipe.crd.both.fixedbudget}
\alias{ss.aipe.crd.both.fixedwidth}
\title{
	Find target sample sizes for the accuracy in unstandardized conditions means estimation in CRD
}
\description{
	Find target sample sizes (the number of clusters, cluster size, or both) for the accuracy in unstandardized conditions means estimation in CRD. If users wish to seek for both types of sample sizes simultaneously, an additional constraint is required, such as a desired width or a desired budget.
}
\usage{
ss.aipe.crd.nclus.fixedwidth(width, nindiv, prtreat, tauy=NULL, sigma2y=NULL, 
	totalvar=NULL, iccy=NULL, r2between = 0, r2within = 0, numpredictor = 0, 
	assurance=NULL, conflevel = 0.95, cluscost=NULL, indivcost=NULL)
ss.aipe.crd.nindiv.fixedwidth(width, nclus, prtreat, tauy=NULL, sigma2y=NULL, 
	totalvar=NULL, iccy=NULL, r2between = 0, r2within = 0, numpredictor = 0, 
	assurance=NULL, conflevel = 0.95, cluscost=NULL, indivcost=NULL)
ss.aipe.crd.nclus.fixedbudget(budget, nindiv, cluscost = 0, indivcost = 1, 
	prtreat = NULL, tauy=NULL, sigma2y=NULL, totalvar=NULL, iccy=NULL, r2between = 0, 
	r2within = 0, numpredictor = 0, assurance=NULL, conflevel = 0.95)
ss.aipe.crd.nindiv.fixedbudget(budget, nclus, cluscost = 0, indivcost = 1, 
	prtreat = NULL, tauy=NULL, sigma2y=NULL, totalvar=NULL, iccy=NULL, r2between = 0, 
	r2within = 0, numpredictor = 0, assurance=NULL, conflevel = 0.95)
ss.aipe.crd.both.fixedbudget(budget, cluscost=0, indivcost=1, prtreat, tauy=NULL, 
	sigma2y=NULL, totalvar=NULL, iccy=NULL, r2between = 0, r2within = 0, 
	numpredictor = 0, assurance=NULL, conflevel = 0.95)
ss.aipe.crd.both.fixedwidth(width, cluscost=0, indivcost=1, prtreat, tauy=NULL, 
	sigma2y=NULL, totalvar=NULL, iccy=NULL, r2between = 0, r2within = 0, 
	numpredictor = 0, assurance=NULL, conflevel = 0.95)
}
\arguments{
  \item{width}{
	The desired width of the confidence interval of the unstandardized means difference
  }
  \item{budget}{
	The desired amount of budget
  }
  \item{nclus}{
	The desired number of clusters
  }
  \item{nindiv}{
	The number of individuals in each cluster (cluster size)
  }
  \item{prtreat}{
	The proportion of treatment clusters
  }
  \item{cluscost}{
	The cost of collecting a new cluster regardless of the number of individuals collected in each cluster
  }
  \item{indivcost}{
	The cost of collecting a new individual
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
\details{
	Here are the functions description:
	\itemize{
		\item \code{ss.aipe.crd.nclus.fixedwidth} Find the number of clusters given a specified width of the confidence interval and the cluster size
		\item \code{ss.aipe.crd.nindiv.fixedwidth} Find the cluster size given a specified width of the confidence interval and the number of clusters
		\item \code{ss.aipe.crd.nclus.fixedbudget} Find the number of clusters given a budget and the cluster size
		\item \code{ss.aipe.crd.nindiv.fixedbudget} Find the cluster size given a budget and the number of clusters
		\item \code{ss.aipe.crd.both.fixedbudget} Find the sample size combinations (the number of clusters and that cluster size) providing the narrowest confidence interval given the fixed budget
		\item \code{ss.aipe.crd.both.fixedwidth} Find the sample size combinations (the number of clusters and that cluster size) providing the lowest cost given the specified width of the confidence interval
	}
}
\value{
	The \code{ss.aipe.crd.nclus.fixedwidth} and \code{ss.aipe.crd.nclus.fixedbudget} functions provide the number of clusters. The \code{ss.aipe.crd.nindiv.fixedwidth} and \code{ss.aipe.crd.nindiv.fixedbudget} functions provide the cluster size. The \code{ss.aipe.crd.both.fixedbudget} and \code{ss.aipe.crd.both.fixedwidth} provide the number of clusters and the cluster size, respectively.
}
\author{
    Sunthud Pornprasertmanit (University of Kansas; \email{psunthud@ku.edu})
}
\examples{
\dontrun{
ss.aipe.crd.nclus.fixedwidth(0.3, 30, 0.5, 0.25, 0.75)
ss.aipe.crd.nindiv.fixedwidth(0.3, 250, 0.5, 0.25, 0.75)
ss.aipe.crd.nclus.fixedbudget(10000, 20, 20, 1)
ss.aipe.crd.nindiv.fixedbudget(10000, 30, 20, 1, prtreat=0.5, tauy=0.05, sigma2y=0.95, assurance=0.8)
ss.aipe.crd.both.fixedbudget(10000, 30, 1, 0.5, 0.25, 0.75)
ss.aipe.crd.both.fixedwidth(0.3, 0, 1, 0.5, 0.25, 0.75)
}
}