\name{upsilon}
\alias{upsilon}
\title{
A function for estimating the mediation effect size as discussed in Lachowicz, Preacher, & Kelley (in press). 
}
\description{
This function implements the upsilon effect size statistic as described in Lachowicz, Preacher, & Kelley (in press).  
}
\usage{
upsilon(x, mediator, dv, conf.level = 0.95, bootstrap = TRUE, B = 1000,
boot.type = "ordinary", ...)
}
\arguments{
  \item{x}{\code{x} is the independent variable. }
  \item{mediator}{\code{mediator} is the mediator variable. }
  \item{dv}{\code{dv} is the outcome or dependent variable. }
  \item{conf.level}{\code{conf.level} is the desired confidence coefficient (i.e., the complement of the Type I error rate).}
  \item{bootstrap}{TRUE or FALSE regarding if a bootstrap confidence interval should be performed.}
  \item{B}{The number of bootstrap replications}
\item{boot.type}{The type of bootstrap confidence interval.}
\item{\dots}{Allows specifictions for functions that are used within this function.}
}
%\details{
%See the examples below for applications of the function. 
%}
\value{
Returns the value of the effect size upsilon for a simple mediation model. 
}
\references{
Lachowicz, M. J., Preacher, K. J., & Kelley, K. (in press). A novel measure of effect size for mediation analysis. \emph{Psychological Methods}, \emph{X}, X--X.  

Preacher, K. J., & Kelley, K. (2011). Effect size measures for mediation models: quantitative strategies for communicating indirect effects. \emph{Psychological Methods}, \emph{16}, 93--115. 

Wen, Z., & Fan, X. (2015). Monotonicity of effect sizes: Questioning kappa-squared as mediation effect size measure. 
\emph{Psychological Methods}, \emph{20}, 193--203. 
}

\author{
Lachowicz
Mark J. Lachowicz (Vanderbilt University; \email{Mark.J.Lachowicz@Vanderbilt.edu})
}
\note{
Note that this function overcomes some limitations of other effects for mediation models, such as those discussed in Preacher and Kelley (2012) and Wen and Fan (2015). This function can only be used for simple mediation models at this time. 
}

\seealso{
\link{mediation}
}
\keyword{design}
