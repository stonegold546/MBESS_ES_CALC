\name{mediation.effect.bar.plot}
\alias{mediation.effect.bar.plot}

\title{ Bar plots of mediation effects}
\description{
Provides an effect bar plot in the context of simple mediation. }
\usage{
mediation.effect.bar.plot(x, mediator, dv, 
main = "Mediation Effect Bar Plot", width = 1, left.text.adj = 0, 
right.text.adj = 0, rounding = 3, file = "", save.pdf = FALSE, 
save.eps = FALSE, save.jpg = FALSE, ...)
}

\arguments{
\item{x}{vector of the predictor/independent variable}
\item{mediator}{vector of the mediator variable}
\item{dv}{vector of the dependent/outcome variable}
  \item{main}{ main title }
  \item{width}{ width of bar, default 1}
  \item{left.text.adj}{ for fine tuning left side text adjustment }
  \item{right.text.adj}{ for fine tuning right side text adjustment }
  \item{rounding}{ how to round so that the values displayed in the plot do not have too few or too many significant digits }
  \item{file}{ file name of the plot to be saved (not necessary) }
  \item{save.pdf}{TRUE or FALSE if the produced figure should be saved as a PDF file}
  \item{save.eps}{TRUE or FALSE if the produced figure should be saved as an EPS file}
  \item{save.jpg}{TRUE or FALSE if the produced figure should be saved as a JPG file}
  \item{\dots}{optional additional specifications for nested functions }
}
\details{
Provides an effect bar for mediation (Bauer, Preacher, & Gil, 2006) may be used to plot the results of a mediation analysis compactly. Effect bars represent, in a single metric, the relative magnitudes of several values that are important for interpreting indirect effects. Preacher and Kelley (2010) discuss this plotting method also.

}
\value{
Only a figure is returned
}
\references{
Bauer, D. J., Preacher, K. J., & Gil, K. M. (2006). Conceptualizing and testing random indirect effects and moderated mediation in multilevel models: New procedures and recommendations. \emph{Psychological Methods, 11}, 142-163.

Preacher, K. J., & Kelley, K. (2010). Effect size measures for mediation models: Quantitative and graphical strategies for communicating indirect effects. \emph{Manuscript in preparation}.
}


\author{Ken Kelley (University of Notre Dame; KKelley@nd.edu)}
%\note{}
\seealso{\code{\link{mediation}}, \code{\link{mediation.effect.bar.plot}}}
%\examples{}

\keyword{multivariate}