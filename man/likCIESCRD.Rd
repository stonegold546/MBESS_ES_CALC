\name{likCIESCRD}
\alias{likCIESCRD}
\title{
	Compute the likelihood-based confidence interval of the effect size in cluster randomized design
}
\description{
	Compute the likelihood-based confidence interval (Cheung, 2009) by running the multilevel model using the structural equation modeling with the data in wide format (Bauer, 2003). This program uses the \code{OpenMx} package (Boker et al., 2011) to run the likelihood-based confidence interval of the effect size. See further details at Pornprasertmanit and Schneider (2010, submitted).
}
\usage{
likCIESCRD(datawide, ylab, xlab, zwlab=NULL, zblab=NULL, estype=1, 
	iccy=0.25, es=0.5, totalvar=1, covariate=FALSE, iccz=0.25, r2within=0.5, 
	r2between=0.5, totalvarz = 1, conflevel=0.95)
}
\arguments{
  \item{datawide}{
	The nested data in wide format
  }
  \item{ylab}{
	The label of the dependent variable
  }
  \item{xlab}{
	The label of the independent variable
  }
  \item{zwlab}{
	The label of the within-level covariate variables
  }
  \item{zblab}{
	The label of the between-level covariate variable
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
  \item{conflevel}{
	The desired level of confidence for the confidence interval
  }
}
\value{
	The likelihood-based confidence interval of effect size
}
\references{ 
Bauer, D. J. (2003). Estimating multilevel linear models as structural equation models. 
\emph{Journal of Educational and Behavioral Statistics, 28}, 135-167.

Boker, S., M., N., Maes, H., Wilde, M., Spiegel, M., Brick, T., et al. (2011).
OpenMx: An open source extended structural equation modeling framework.
\emph{Psychometrika, 76,} 306-317.

Cheung, M. W.-L. (2009). Constructing approximate confidence intervals for
parameters with structural constructing approximate confidence intervals for
parameters with structural equation models. \emph{Structural Equation Modeling,
16,} 267-294.

Pornprasertmanit, S., & Schneider, W. J. (2010). \emph{Efficient sample size for power and desired accuracy in Cohen's d estimation in two-group cluster randomized design} (Master Thesis). Illinois State University, Normal, IL. 

Pornprasertmanit, S., & Schneider, W. J. (submitted). \emph{Accuracy in parameter estimation in two-condition cluster randomized design}.
}
\author{
    Sunthud Pornprasertmanit (University of Kansas; \email{psunthud@ku.edu})
}
\examples{
\dontrun{
datWide <- createDataCRDWide(100, 50, 5, iccy=0.25, es=0.5, estype = 0, 
	totalvar=1, covariate=FALSE)
likCIESCRD(datWide, ylab=c("y1", "y2", "y3", "y4", "y5"), xlab="x", 
	estype=0, iccy=0.25, es=0.5, totalvar=1, covariate=FALSE)
}
}