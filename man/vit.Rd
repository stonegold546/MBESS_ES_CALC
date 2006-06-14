\name{vit}
\alias{vit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Visualize individual trajectories}
\description{
A function to help visualize individual trajectories in a longitudinal (i.e., analysis of change) context.
}
\usage{
vit(id = "", occasion = "", score = "", Data = NULL, group = NULL, subset.ids = NULL, 
pct.rand = NULL, All.in.One = TRUE,  ylim = NULL, xlim = NULL, ylab = "Score", 
xlab = "Occasion", main = "", plot.points = TRUE, draw.which.lines = "Observed", 
predicted.id = NULL, predicted.occasion = NULL, predicted.score = NULL, 
predicted.data = NULL, pch = 16, points.cex = 0.7, lty = 1, lwd = 1, 
per.page.layout = rbind(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9), c(10, 11, 12)), 
points.col = "Black", lines.col = "Black", mar = c(2.2, 2, 1, 1), 
mgp = c(2.2, 1, 0), pty = "s", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{string variable of the column name of id}
  \item{occasion}{string variable of the column name of id}
  \item{score}{string variable of the column name where the score (i.e., dependent variable) is located}
  \item{Data}{data set with named column variables (see above)}
  \item{group}{if plotting parameters should be conditional on group membership}
  \item{subset.ids}{\code{id} values for a selected subset of individuals}
  \item{pct.rand}{percentage of random trajectories to be plotted}
  \item{All.in.One}{should trajectories be in a single or multiple plots}
  \item{ylim}{optional limit of the ordinate (i.e., \code{y}-axis; see \code{par})}
  \item{xlim}{optional limit for the abscissa (i.e., \code{x}-axis; see \code{par})}
  \item{ylab}{label for the ordinate (i.e., \code{y}-axis; see \code{par})}
  \item{xlab}{label for the abscissa (i.e., \code{x}-axis; see \code{par})}
  \item{main}{main title of the plot}
  \item{plot.points}{should the points be plotted}
  \item{draw.which.lines}{draw the observed or predicted lines}
  \item{predicted.id}{like \code{id} but for the predicted data}
  \item{predicted.occasion}{like \code{occasion} but for the predicted data}
  \item{predicted.score}{like \code{score} but for the predicted scores}
  \item{predicted.data}{like \code{Data} but for the predicted data}
  \item{pch}{plotting character(s); see \code{par}}
  \item{points.cex}{size of the points (1 is the R default; see \code{par})}
  \item{lty}{type of line to be plotted}
  \item{lwd}{width of line to be plotted}
  \item{per.page.layout}{define the per-page layout when \code{All.in.One==FALSE}}
  \item{points.col}{color(s) of the points}
  \item{lines.col}{color(s) of the line(s)}
  \item{mar}{adjusts the margin spacing (see \code{par})}
  \item{mgp}{adjusts the margin lines (see \code{mgp})}
  \item{pty}{plot region type (see \code{par}; "s" is square)}
  \item{\dots}{optional plotting specifications}
}
\details{
This function makes visualizing individual trajectories simple. If the plot will be saved for a publication, 
it is best to call upon a device so that the file can be saved. See \code{Devices} (e.g., use \code{pdf} to 
save a PDF file of the figure). Data should be in the "univariate format" (i.e., the same format as 
\code{lmer} and \code{nlme} data.

}
\value{
Returns a plot of individual trajectories with the specifications provided.
}
% \references{}
\author{Ken Kelley (Indiana University; \email{KKIII@Indiana.Edu})}
\note{
At present, the use of predicted statements is limited (and not recommended).   
}
\seealso{\code{par}, \code{Devices}, \code{pdf}, \code{jpeg}, \code{bmg}}
\examples{
data(Gardner.LD)

# Although many options are possible, a simple call to \code{vit} is of the form:
vit(id="ID", occasion="Trial", score="Score", Data=Gardner.LD)

# Now line color is conditional on group membership.
vit(id="ID", occasion="Trial", score="Score", Data=Gardner.LD, 
group="Group", draw.which.lines="Observed", pch=16, 
points.col=c("Black"), lines.col=c("Black", "Blue"))

# Now line and plotting character are conditional on group membership
vit(id="ID", occasion="Trial", score="Score", Data=Gardner.LD, 
group="Group", draw.which.lines="Observed",  pch=16, 
points.col=c("Black", "Blue"), lines.col=c("Black", "Blue"))

# Now plotting character is different with same color lines.
vit(id="ID", occasion="Trial", score="Score", Data=Gardner.LD, 
group="Group", draw.which.lines="Observed",  pch=c(1, 16), 
points.col=c("Black"), lines.col=c("Black"))

# Point color and line color can be mixed (this feature may not 
# be useful very often, but the function is general nonetheless).
vit(id="ID", occasion="Trial", score="Score", Data=Gardner.LD, 
group="Group", draw.which.lines="Observed",  pch=16, 
points.col=c("Red", "Orange"), lines.col=c("Black", "Blue"))

# Now randomly selects 50\% of the sample to plot
vit(id="ID", occasion="Trial", score="Score", Data=Gardner.LD, 
pct.rand=50, group="Group", draw.which.lines="Observed",  pch=16)

# Specified individuals are plotted (by group)
vit(id="ID", occasion="Trial", score="Score", Data=Gardner.LD, 
subset.ids=c(1, 4, 8, 13, 17, 21), group="Group", 
lines.col=c("Black", "Blue"), pch=16)

}
\keyword{hplot}
\keyword{device}
\keyword{dynamic}
