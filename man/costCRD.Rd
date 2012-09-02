\name{costCRD}
\alias{costCRD}
\title{
	Calculate the total cost used in the cluster randomized design
}
\description{
	Calculate the total cost used in the cluster randomized design. See further details at Pornprasertmanit and Schneider (submitted).
}
\usage{
costCRD(nclus, nindiv, cluscost, indivcost)
}
\arguments{
  \item{nclus}{
	The total number of clusters
  }
  \item{nindiv}{
	The number of individuals in each cluster (cluster size)
  }
  \item{cluscost}{
	The cost of collecting a new cluster regardless of the number of individuals collected in each cluster
  }
  \item{indivcost}{
	The cost of collecting a new individual
  }
}
\value{
	The total cost used in the cluster randomized design
}
\references{ 
Pornprasertmanit, S., & Schneider, W. J. (submitted). \emph{Accuracy in Parameter Estimation in Two-Condition Cluster Randomized Design}.
}
\author{
    Sunthud Pornprasertmanit (University of Kansas; \email{psunthud@ku.edu})
}
\examples{
\dontrun{
costCRD(500, 10, 5 ,1)
}
}
