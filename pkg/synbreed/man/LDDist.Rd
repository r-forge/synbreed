\name{LDDist}
\Rdversion{1.1}
\alias{LDDist}

\title{
LD versus distance Plot
}
\description{
Computation of pairwise LD as \eqn{R^2} for marker data and plot LD versus distance for each linkage group (chromosome)
}
\usage{
LDDist(marker,linkageGroup,pos,file=NULL,type="p",breaks=NULL,chr=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{marker}{
A \code{data.frame} of \code{matrix} containing the marker data. Could be the result of a call of \code{codeGeno}.
}
  \item{linkageGroup}{
A numeric vector of length \code{nrow(marker)} with the linkage group of each marker.
}
  \item{pos}{
A numeric vector of length \code{nrow(marker)} with physical position of each marker on the chromosome.
}
  \item{file}{
Optionally a path were the plot is saved as \code{pdf}.
}
  \item{type}{
Character string to specify the type of plot. Use \code{"p"} for a scatterplot, \code{"bars"} for stacked bars and \code{"nls"} for scatterplot together with nonlinear regression curve estimated by data. 
}
  \item{breaks}{
List containing breaks for stacked bars. Components are \code{dist} with breaks for distance on x-axis and \code{r2} for breaks on for r2 on y-axis. By default, 5 equal spaced categories for dist and r2 are used.  
}
  \item{chr}{
Character or numeric. If specified, only one plot for markers in linkage group \code{chr} is created.  
}
  \item{\dots}{
Further plot arguments
}
}
\details{
LD is computed as coefficient of determination \eqn{R^2}.  Matrix of \eqn{R^2} of markers is returned. By default, one plot for each linkage group is generated. If markers across all chromosomes should be combined
in one plot, use \code{chr="all"}. If markers of, e.g. chromosomes 5 to 10 should be combined in one plot, change \code{linkageGroup} indicator manually.
}

\seealso{
\code{\link{codeGeno}}, \code{\link{LDMap}}
}
\examples{
# simulated data without LD with 1000 equal spaced markers on 4 chromosomes
M<- matrix(sample(c(0,2),size=50000,replace=TRUE,prob=c(0.65,0.35)),nrow=50,ncol=1000)
colnames(M) <- paste("M",1:1000,sep="")
lgd <- rep(1:4,each=250)
pos <- sequence(c(250,250,250,250))
par(mfrow=c(2,2))
LDDist(M,lgd,pos)
par(mfrow=c(1,1))

# maize data example
data(maize)
marker <- codeGeno(maize.geno[,-1])

# scatterplot for chromosome 1
LDDist(marker,maize.marker.pos$chr,maize.marker.pos$pos,type="p",chr=1)

# stacked bars for chromosome 1
LDDist(marker,maize.marker.pos$chr,maize.marker.pos$pos,type="bars",chr=1)



}

\keyword{LD}
\keyword{genotypic data}
