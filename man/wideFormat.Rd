\name{wideFormat}
\alias{wideFormat}
\title{
	Transform nested data from long format into wide format
}
\description{
	Transform nested data from long format into wide format to run multilevel model in structural equation modeling framework (Bauer, 2003).
}
\usage{
wideFormat(data, betweencol, withincol, idcol)
}
\arguments{
  \item{data}{
	A long-format dataset
}
  \item{betweencol}{
	The columns of the between-level variables
  }
  \item{withincol}{
	The columns of the within-level variables
  }
  \item{idcol}{
	The variable specifying cluster ID
  }
}
\value{
	A dataset in wide format
}
\references{
Bauer, D. J. (2003). Estimating multilevel linear models as structural equation models. 
\emph{Journal of Educational and Behavioral Statistics, 28}, 135-167.
}
\author{
    Sunthud Pornprasertmanit (University of Kansas; \email{psunthud@ku.edu})
}
\examples{
\dontrun{
dat <- createDataCRD(100, 50, 4, iccy=0.25, es=0.5)
datwide <- wideFormat(dat, 3, 2, 1)
}
}