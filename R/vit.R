"vit" <-
function(id="", occasion="", score="", Data=NULL, group=NULL, subset.ids=NULL, 
pct.rand=NULL, All.in.One=TRUE,  ylim=NULL, xlim=NULL, ylab="Score", xlab="Occasion", main="",
plot.points=TRUE, draw.which.lines="Observed", predicted.id=NULL, predicted.occasion=NULL, 
predicted.score=NULL, predicted.data=NULL, pch=16, points.cex=.7, lty=1, lwd=1, 
per.page.layout=rbind(c(1,2,3),c(4,5,6),c(7,8,9),c(10,11,12)), points.col="Black", 
lines.col="Black", mar=c(2.2,2,1,1), mgp=c(2.2, 1, 0), pty="s", ...)
{

if(is.null(Data)) stop("You need to specify the data set with the \'Data\' argument.")

if(is.data.frame(Data))
{
id <- Data[,which(names(Data)==id)]
occasion <- Data[,which(names(Data)==occasion)]
score <- Data[,which(names(Data)==score)]
if(!is.null(group)) group <- Data[,which(names(Data)==group)]

}

if(!is.data.frame(Data))
{
id <- Data[,(colnames(Data)==id)]
occasion <- Data[,(colnames(Data)==occasion)]
score <- Data[,(colnames(Data)==score)]
if(!is.null(group)) group <- Data[,(colnames(Data)==group)]
}


if(!is.null(subset.ids))
{

if(!is.null(pct.rand)) stop("Since \'subset.ids\' was specified, do not also specify \'pct.random\'.")

Rows.to.Select <- matrix(NA, length(id), length(subset.ids))
for(i in 1:length(subset.ids))
{
Rows.to.Select[,i] <- id==subset.ids[i]
}
Rows.to.Select <- apply(Rows.to.Select, MARGIN=1, sum)==1

id <- id[Rows.to.Select]  
occasion <- occasion[Rows.to.Select] 
score <- score[Rows.to.Select]
if(!is.null(group)) group <- group[Rows.to.Select]
}

if(!is.null(group)) group <- unique(cbind(id, group))[,2]

if(!is.null(pct.rand))
{
if(pct.rand <=1 & pct.rand>0) {rand.samp <- sample(unique(id), ceiling(length(unique(id))*pct.rand), replace=FALSE)}
if(pct.rand >1 & pct.rand<=100) {rand.samp <- sample(unique(id), ceiling(length(unique(id))*pct.rand/100), replace=FALSE)}

Data.From.Rand <- matrix(NA, length(id), length(rand.samp))
for(i in 1:length(rand.samp))
{
Data.From.Rand[,i] <- id==rand.samp[i]
}

Rows.to.Select <- apply(Data.From.Rand, MARGIN=1, sum)==1
id <- id[Rows.to.Select]  
occasion <- occasion[Rows.to.Select] 
score <- score[Rows.to.Select]
if(!is.null(group)) group <- group[Rows.to.Select]
}

ID <- unique(id)
N <- length(ID)

# Up until here, the function has only organized that data that is to be used. 
####################################################################################

####################################
if(plot.points==TRUE)
{
if(!is.null(group)) if(!((length(pch)==1) | (length(pch)==N) | (length(pch)==length(unique(group))))) stop("The vector 'pch' must be of length 1, length N, or of a length equal to the number of groups.")
if(is.null(group)) if(!((length(pch)==1) | (length(pch)==N))) stop("The vector 'pch' must be of length 1 or length N (or of a length equal to the number of groups).")

if(!is.null(pch))
{
if(length(pch) == 1)
{
pch <- rep(pch, N)
}

if(!is.null(group) & length(pch)==1)
{
pch <- rep(pch, N)
}

if(!is.null(group) & length(pch)>1)
{
tmp.MBESS <- rep(NA, N)
which.pch <- rep(NA, N)
for(i in 1:N)
{
tmp.MBESS[i] <- which(group[i]==unique(group))
which.pch[i] <- pch[tmp.MBESS[i]]
}
pch <- which.pch
}

}
}

####################################
if(plot.points==TRUE)
{

if(!is.null(group)) if(!((length(points.col)==1) | (length(points.col)==N) | (length(points.col)==length(unique(group))))) stop("The vector 'points.col' must be of length 1, length N, or of a length equal to the number of groups.")
if(is.null(group)) if(!((length(points.col)==1) | (length(points.col)==N))) stop("The vector 'points.col' must be of length 1 or length N (or of a length equal to the number of groups).")

if(!is.null(points.col))
{
if(length(points.col) == 1)
{
points.col <- rep(points.col, N)
}


if(!is.null(group))
{
tmp.MBESS <- rep(NA, N)
which.points.col <- rep(NA, N)
for(i in 1:N)
{
tmp.MBESS[i] <- which(group[i]==unique(group))
which.points.col[i] <- points.col[tmp.MBESS[i]]
}
points.col <- which.points.col
}

}
}

#####################################

if(!is.null(group)) if(!((length(lines.col)==1) | (length(lines.col)==N) | (length(lines.col)==length(unique(group))))) stop("The vector 'lines.col' must be of length 1, length N, or of a length equal to the number of groups.")
if(is.null(group)) if(!((length(lines.col)==1) | (length(lines.col)==N))) stop("The vector 'lines.col' must be of length 1 or length N (or of a length equal to the number of groups).")

if(is.null(group))
{
if(length(lines.col)==1)
lines.col <- rep(lines.col, N)
if(length(lines.col)!=N) stop("There is a problem with the vector \'lines.col\' (does not match N or the number of groups).")
}

if(!is.null(group) & length(lines.col)==1)
{
lines.col <- rep(lines.col, N)
}

if(!is.null(group) & length(lines.col) > 1)
{
tmp.MBESS <- rep(NA, N)
which.line.col <- rep(NA, N)
for(i in 1:N)
{
tmp.MBESS[i] <- which(group[i]==unique(group))
which.line.col[i] <- lines.col[tmp.MBESS[i]]
}
lines.col <- which.line.col
}

###########################

if(!is.null(group)) if(!((length(lty)==1) | (length(lty)==N) | (length(lty)==length(unique(group))))) stop("The vector 'lty' must be of length 1, length N, or of a length equal to the number of groups.")
if(is.null(group)) if(!((length(lty)==1) | (length(lty)==N))) stop("The vector 'lty' must be of length 1 or length N (or of a length equal to the number of groups).")

if(is.null(group))
{
if(length(lty)==1)
lty <- rep(lty, N)
if(length(lty)!=N) stop("There is a problem with the vector \'lty\' (does not match N or the number of groups).")
}

if(!is.null(group) & length(lty)==1)
{
lty <- rep(lty, N)
}

if(!is.null(group) & length(lty) > 1)
{
tmp.MBESS <- rep(NA, N)
which.lty <- rep(NA, N)
for(i in 1:N)
{
tmp.MBESS[i] <- which(group[i]==unique(group))
which.lty[i] <- lty[tmp.MBESS[i]]
}
lty <- which.lty
}

# Optionally sets up the limits of the axes automatically.
if(is.null(ylim)) ylim <- c(min(score, na.rm=TRUE), max(score, na.rm=TRUE))
if(is.null(xlim)) xlim <- c(min(occasion, na.rm=TRUE), max(occasion, na.rm=TRUE))
####################################################################################
if(!is.null(predicted.data)) 
{
predicted.id <- predicted.data[,(names(predicted.data)==predicted.id)]
predicted.occasion <- predicted.data[,(names(predicted.data)==predicted.occasion)]
predicted.score <- predicted.data[,(names(predicted.data)==predicted.score)]
}

# Draws the plotting region.
if(All.in.One==FALSE)
{
layout(per.page.layout)
par(pty=pty, mar=mar, mgp=mgp, ...)
####################################################################################
# The following optionally plots the values of the observed data.
print("In order to display multiple graphic windows (e.g., if there are more than 12 participants), click \'Recording\' in the \'History\' meno of an open plot window. Using \'Page Up\' and \'Page Down\' will scroll through the graphic windows.")
for(i in 1:length(ID))
{
if(is.null(ylab)) plot(occasion[id==ID[i]], score[id==ID[i]], ylim=ylim, xlim=xlim, ylab=paste("ID=", ID[i], sep=""), xlab=xlab, font.lab=3, type="n", main=main, ...)
if(!is.null(ylab)) plot(occasion[id==ID[i]], score[id==ID[i]], ylim=ylim, xlim=xlim, ylab=ylab, xlab=xlab, font.lab=3, type="n", main=main, ...)
if(plot.points==TRUE) points(occasion[id==ID[i]], score[id==ID[i]], pch=pch[i], cex=points.cex, col=points.col[i], ...)
if(draw.which.lines=="Observed" | draw.which.lines=="observed") lines(occasion[id==ID[i]], score[id==ID[i]], pch=pch[i], lty=lty[i], lwd=lwd, col=lines.col[i], ...)
if(draw.which.lines=="Predicted" | draw.which.lines=="predicted") lines(predicted.occasion[predicted.id==ID[i]], predicted.score[predicted.id==ID[i]], pch=pch[i], lty=lty[i], lwd=lwd, col=lines.col[i], ...)
####################################################################################
}
}

if(All.in.One==TRUE)
{
plot(0, ylim=ylim, xlim=xlim, ylab=ylab, xlab=xlab, font.lab=3, type="n", main=main, ...)
for(i in 1:N)
{
if(plot.points==TRUE) points(occasion[id==ID[i]], score[id==ID[i]], pch=pch[i], cex=points.cex, col=points.col[i], ...)
if(draw.which.lines=="Observed" | draw.which.lines=="observed") lines(occasion[id==ID[i]], score[id==ID[i]], pch=pch[i], lty=lty[i], lwd=lwd, col=lines.col[i], ...)
if(draw.which.lines=="Predicted" | draw.which.lines=="predicted") lines(predicted.occasion[predicted.id==ID[i]], pch=pch[i], predicted.score[predicted.id==ID[i]], lty=lty[i], lwd=lwd, col=lines.col[i], ...)
####################################################################################
}
}
}
