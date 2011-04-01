\name{LDDist}
\Rdversion{1.1}
\alias{LDDist}

\title{
LD versus distance Plot
}
\description{
Computation of pairwise LD measured as \eqn{r^2} for markers pairs within an object of class \code{gpData}. Function creates separate plots for every chromosome or whole genome. 
}
\usage{
LDDist(gpData,chr=NULL,type="p",breaks=NULL,file=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData}{
object of class \code{gpData} where \code{geno} and \code{map} are available.
}
  \item{chr}{
\code{character} or \code{numeric}. If specified, only one plot for markers in linkage group \code{chr} is created. If \code{chr="all"}, all pairwise LD between markers across chromosomes is used.  
}
  \item{type}{
Character string to specify the type of plot. Use \code{"p"} for a scatterplot, \code{"bars"} for stacked bars or \code{"nls"} for scatterplot together with nonlinear regression curve according to Hill and Weir (1988). 
}
  \item{breaks}{
\code{list} containing breaks for stacked bars (only for \code{type="bars"}). Components are \code{dist} with breaks for distance on x-axis and \code{r2} for breaks on for r2 on y-axis. By default, 5 equal spaced categories for dist and r2 are used.  
}
  \item{file}{
\code{character}. path to a file were plot is saved as \code{pdf} (optional).
}
  \item{\dots}{
Further arguments for \code{plot}
}
}
\details{
LD is computed as coefficient of determination \eqn{r^2}.  Matrix of \eqn{r^2} of markers is returned. Missing values in genotypic data are allowed and pairwise complete observations of two markers are used to compute \eqn{r^2}.  By default, one plot for each linkage group is generated. If markers across all chromosomes should be combined
in one plot, use \code{chr="all"}.
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

# stacked bars for chromosome 1 with user-defined categories
LDDist(maize,type="bars",chr=1,breaks=list(dist=c(0,10,20,40,60,180),
r2=c(1,0.6,0.4,0.3,0.1,0)))

}
\keyword{hplot}

