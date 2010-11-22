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
LDDist(gpData,chr=NULL,type="p",breaks=NULL,file=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData}{
Object of class \code{gpData} where \code{geno} and \code{map} are available.
}
  \item{chr}{
Character or numeric. If specified, only one plot for markers in linkage group \code{chr} is created. If \code{chr="all"}, all pairwise LD between markes across chromosomes is used.  
}
  \item{type}{
Character string to specify the type of plot. Use \code{"p"} for a scatterplot, \code{"bars"} for stacked bars or \code{"nls"} for scatterplot together with nonlinear regression curve estimated by data. 
}
  \item{breaks}{
List containing breaks for stacked bars. Components are \code{dist} with breaks for distance on x-axis and \code{r2} for breaks on for r2 on y-axis. By default, 5 equal spaced categories for dist and r2 are used.  
}
  \item{file}{
Optionally a path to a file were plot is saved as \code{pdf}.
}
  \item{\dots}{
Further arguments for \code{plot}
}
}
\details{
LD is computed as coefficient of determination \eqn{R^2}.  Matrix of \eqn{R^2} of markers is returned. Missing values in genotypic data are allowed and pairwise complete observations of two markers are used to compute \eqn{R^2}.  By default, one plot for each linkage group is generated. If markers across all chromosomes should be combined
in one plot, use \code{chr="all"}. If markers of, e.g. chromosomes 5 to 10 should be combined in one plot, change \code{linkageGroup} indicator manually.
}
\references{
For nonlinear regression curve:
Hill WG, Weir BS (1988) Variances and covariances of squared linkage disequilibria in finite populations. Theor Popul Biol 33:54-78.
}
\seealso{
\code{\link{codeGeno}}, \code{\link{LDMap}}
}
\examples{
# maize data example
data(maize)
maize <- codeGeno(maize)

# scatterplot for chromosome 1
LDDist(maize,type="p",chr=1,pch=19,col=hsv(alpha=0.1,v=0))

# stacked bars for chromosome 1 with default categories
LDDist(maize,type="bars",chr=1)

# stacked bars for chromosome 1 with specified categories
LDDist(maize,type="bars",chr=1,breaks=list(dist=c(0,10,20,40,60,180),r2=c(1,0.6,0.4,0.3,0.1,0)))

}
\keyword{hplot}

