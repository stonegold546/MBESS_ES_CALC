\name{F.and.R2.Noncentral.Conversion}
\alias{Rsquare2F}
\alias{Rsquare2Lambda}
\alias{F2Rsquare}
\alias{Lambda2Rsquare}

\title{Conversion functions from noncentral noncentral values to their corresponding
and vice versa, for those related to the F-test and R Square.}
\description{
Given values of test statistics (and the appropriate additional information) the value of the noncentral
values can be obtained. Likewise, given noncentral values (and the appropriate additional information)
the value of the test statistic can be obtained. 
}
\usage{
Rsquare2F(R2 = NULL, df.1 = NULL, df.2 = NULL, p = NULL, N = NULL)
Lambda2Rsquare(Lambda = NULL, N = NULL)
Rsquare2Lambda(R2 = NULL, N = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R2}{squared multiple correlation coefficient (population or observed)}
  \item{df.1}{degrees of freedom for the numerator of the \code{F}-distribution}
  \item{df.2}{degrees of freedom for the denominator of the \code{F}-distribution}
  \item{p}{number of predictor variables for \code{R2}}
  \item{N}{sample size}
  \item{Lambda}{The noncentral parameter from an \code{F}-distribution}
}
\details{
These functions are especially helpful for forming confidence intervals around noncentral parameters, as they 
convert to and from related quantities.}
\value{
Returns the converted value from the specified function.
}
%\references{}
\author{Ken Kelley (Indiana University; \email{KKIII@Indiana.Edu})}

%\note{ ~~further notes~~ 
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{'ss.aipe.R2', 'ci.R2', 'conf.limits.nct', 'conf.limits.ncf'}
\examples{Rsquare2Lambda(R2=.5, N=100)}
\keyword{multivariate}
\keyword{design}
