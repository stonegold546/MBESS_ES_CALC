\name{createDataCRD}
\alias{createDataCRD}
\alias{createDataCRDWide}
\title{
	Create the cluster-randomized design data with two conditions
}
\description{
	Create the cluster-randomized design data with two conditions. If a covaraite is specified, the covariate will not correlate with the grouping variable.
}
\usage{
createDataCRD(nclus, ntreatclus, nindiv, iccy, es, estype = 1, totalvar=1, 
	covariate=FALSE, iccz=NULL, r2within=NULL, r2between=NULL, totalvarz = 1)
createDataCRDWide(nclus, ntreatclus, nindiv, iccy, es, estype = 1, totalvar=1, 
	covariate=FALSE, iccz=NULL, r2within=NULL, r2between=NULL, totalvarz = 1)
}
\arguments{
  \item{nclus}{
	The number of clusters
  }
  \item{ntreatclus}{
	The number of clusters in the treatment condition
  }
  \item{nindiv}{
	The number of individuals within each group (cluster size)
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
	The intraclass correlation of the covariate (used when \code{covariate = TRUE}). If \code{iccz = 0}, the within-level covariate will be only created. If \code{iccz = 1}, the between-level covariate will be only created.
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
}
\value{
	A dataset in the long format if the \code{createDataCRD} function is used and in the wide format if the \code{createDataCRDWide} function is used.
}
\author{
    Sunthud Pornprasertmanit (University of Kansas; \email{psunthud@ku.edu})
}
\examples{
\dontrun{
dat <- createDataCRD(100, 50, 4, iccy=0.25, es=0.5)
datWide <- createDataCRDWide(100, 50, 5, iccy=0.25, es=0.5, estype = 0, totalvar=1, 
	covariate=TRUE, iccz=0.25, r2within=0.5, r2between=0.5, totalvarz = 1)
}
}