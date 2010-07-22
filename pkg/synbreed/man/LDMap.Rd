\name{LDMap}
\Rdversion{1.1}
\alias{LDMap}

\title{
LD Heatmap
}
\description{
Computation of pairwise LD for marker data and plotting LD heatmaps for each linkage group (chromosome)
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
A numeric vector of length \code{nrow(marker)} with physical position of the marker on the chromosome.
}
  \item{file}{
Optional a file were the plot is saved as \code{pdf}.
}
  \item{\dots}{
Further arguments that could be passed to function \code{LDheatmap}.
}
}
\details{
LD is computed as coefficient of determination \eqn{R^2}. The plot is simply a call of function \code{LDheatmap} in package  \code{LDheatmap}.
}

\seealso{
\code{\link{codeGeno}}, \code{\link{LDheatmap}}
}
\examples{
data(maize)
LDMap(maize.geno[,-1],maize.marker.pos$chr,maize.marker.pos$pos)
}

\keyword{plotting}
\keyword{genotypic data}% __ONLY ONE__ keyword per line
