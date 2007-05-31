\name{vit.fitted}
\alias{vit.fitted}
\title{ Visualize individual trajectories with fitted curve and quality of fit }
\description{
  A function to help visualize individual trajectories in a longitudinal (i.e., analysis of change) context with fitted curve
  and quality of fit after analyzing the data with lme, lmer, or nlme function.
}
\usage{
vit.fitted(fit.Model, layout = c(3, 3), ylab = "", xlab = "", 
pct.rand = NULL, number.rand = NULL, subset.ids = NULL, 
same.scales = TRUE, save.pdf = FALSE, save.eps = FALSE, 
save.jpg = FALSE, file = "", ...)
}
\arguments{
  \item{fit.Model}{ lme, nlme object produced by nlme package or lmer object produced by lme4 package }
  \item{layout}{ define the per-page layout when All.in.One==FALSE }
  \item{ylab}{ label for the ordinate (i.e., y-axis; see par) }
  \item{xlab}{ label for the abscissa (i.e., x-axis; see par) }
  \item{pct.rand}{ percentage of random trajectories to be plotted }
  \item{number.rand}{ number of random trajectories to be plotted }
  \item{subset.ids}{ id values for a selected subset of individuals to be plotted }
  \item{same.scales}{ should the y-axes have the same scales }
  \item{save.pdf}{ save a pdf file }
  \item{save.eps}{ save a postscript file }
  \item{save.jpg}{ save a jpg file }
  \item{file}{ file name and file path for the graph(s) to save, if file="" a file would be saved in the current working directory }
  \item{\dots}{ optional plotting specifications }
}
\details{
  This function use the fitted model from nlme and lme functions in nlme package, and lmer function in lme4 package. 
  It returns a set of plots of individual observed data, the fitted curves and the quality of fit.
}
%\value{}

%\references{ }

\author{ Ken Kelley (Indiana University; \email{KKIII@indiana.edu}) and Po-Ju Wu (Indiana University; \email{pojwu@indiana.edu}) }

%\note{ }

\seealso{ par, nlme, lme4, lme, lmer, vit.fitted }

\examples{

data(Gardner.LD)
library(nlme)
grouped.Gardner.LD <- groupedData(Score ~ Trial|ID, data=Gardner.LD, order.groups=FALSE)    
G.L.nlsList<- nlsList(SSlogis,grouped.Gardner.LD)
G.L.nlme <- nlme(G.L.nlsList)

vit.fitted(G.L.nlme)
vit.fitted(G.L.nlme, pct.rand = 50)
}

\keyword{hplot}
\keyword{device}
\keyword{dynamic}
