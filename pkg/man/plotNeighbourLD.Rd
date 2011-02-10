\name{plotNeighbourLD}
\alias{plotNeighbourLD}
\title{
Plot neighbour linkage disequilibrium
}
\description{
A function to visualize LD between adjacent markers.
}
\usage{
plotNeighbourLD(LD, map, nMarker = TRUE, dense=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LD}{
object of class \code{list} containing pairwise LD matrices for each chromosome, e.g. output of function \code{LDMap}.
}
  \item{map}{
object of class \code{gpData} with object \code{map} or a \code{data.frame} with columns 'chr' (specifying the chromosome of the marker) and 'pos' (position of the marker within chromosome measured with genetic or physical distances).
}
  \item{nMarker}{
\code{logical}. Add number of markers?
}
  \item{dense}{                                       
\code{logical}. Should visualization for high-density visualization be used?
}

  \item{\dots}{
further graphical arguments for function \code{plot}
}
}
\details{
The graph is similar to \code{plotGenMap} with the option \code{dense=TRUE}, but here the LD between adjacent markers is plotted along the chromosomes. 
}
\value{
Plot of smoothed neighbour LD along each chromosome. One chromosome is displayed from the first to the last marker.
}

\author{
Theresa Albrecht
}
\note{
If the output of \code{LDMap} is used, the function \code{LDMap} has to be allocated to an object.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plotGenMap}}, \code{\link{LDMap}}
}
\examples{
data(maize)
maize2 <-codeGeno(maize)
LD <- LDMap(maize2,chr=1:10)
plotNeighbourLD(LD,maize2,nMarker=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

