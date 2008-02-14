\name{ss.aipe.rmsea}
\alias{ss.aipe.rmsea}
\title{ Sample siza planning for population root mean square 
error of approximation }
\description{
Sample siza planning for population root mean square error of approximation 
from the Accuracy in Parameter Estimation (AIPE) perspective
}
\usage{
ss.aipe.rmsea(RMSEA, df, width, conf.level = 0.95)
}
\arguments{
  \item{RMSEA}{population root mean square error of approximation }
  \item{df}{degrees of freedom of the model}
  \item{width}{desired confidence interval width}
  \item{conf.level}{desired confidence level (e.g., .90, .95, .99)}
}
%\details{}
\value{
Returns the necessary total sample size in order to achieve the desired 
degree of accuracy (i.e., the sufficiently narrow confidence interval). 
}
%\references{  }
\author{Ken Kelley (Indiana University; \email{KKIII@Indiana.Edu})}
%\note{}
\seealso{ \code{\link{ci.rmsea}}}
\examples{
ss.aipe.rmsea(RMSEA=.035, df=50, width=.05, conf.level=.95)
}

\keyword{ design }
