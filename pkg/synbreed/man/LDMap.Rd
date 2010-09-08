\name{LDMap}
\Rdversion{1.1}
\alias{LDMap}

\title{
LD Heatmap
}
\description{
Computation of pairwise LD as \eqn{R^2} for marker data and plot of LD heatmap for each linkage group (chromosome)
}
\usage{
LDMap(marker, linkageGroup, pos, file = NULL, ...)
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
  \item{\dots}{
Further arguments that could be passed to function \code{LDheatmap}.
}
}
\details{
LD is computed as coefficient of determination \eqn{R^2}. The plot is simply a call of function \code{LDheatmap} in package  \code{LDheatmap}.  
}
\value{
A \code{matrix} containing all pairwise values of \eqn{R^2} for all markers across linkage groups. 
}

\seealso{
\code{\link{codeGeno}}, \code{\link{LDheatmap}}, \code{\link{LDDist}}
}
\examples{
data(maize)
LDMap(maize.geno[,-1],maize.marker.pos$chr,maize.marker.pos$pos)
}

\keyword{LD}
\keyword{genotypic data}% __ONLY ONE__ keyword per line
