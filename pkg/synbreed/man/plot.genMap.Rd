\name{plotGenMap}
\alias{plotGenMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot marker map
}
\description{
A function to visualize high-density marker maps.
}
\usage{
plotGenMap(map, dense = FALSE, nMarker = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{map}{
object of class \code{gpData} with object \code{map} or a \code{data.frame} with columns 'chr' (specifying the chromosome of the marker) and 'pos' (position of the marker within chromosome measured with genetic or physical distances)
}                                                                                  
  \item{dense}{
\code{logical}. Should visualization for high-density genetic maps be used?
}
  \item{nMarker}{
\code{logical}. Add number of markers?
}

  \item{\dots}{
further graphical arguments for function \code{plot}
}
}
\details{
In low density plot, single position of markers are drawn. In high-density plot, the distribution of the markers is visualized as a heatmap of density estimation. In this case, 0 is lowest observed density for single chromosome and 1 highest. A colorkey is added. 
}
\value{
Plot of the marker positions within  each chromosome. One chromosome is displayed from the first to the last marker.
}

\author{
Valentin Wimmer
}


\seealso{
\code{\link{create.gpData}}
}
\examples{
# low density plot
data(maize)
plotGenMap(maize)

# high density plot  
data(mice)
plotGenMap(mice,dense=TRUE,nMarker=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
