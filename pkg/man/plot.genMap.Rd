\name{plotGenMap}
\alias{plotGenMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot marker map
}
\description{
A function to visualize low and high-density marker maps.
}
\usage{
plotGenMap(map, dense = FALSE, nMarker = TRUE, bw=1, centr=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{map}{
object of class \code{gpData} with object \code{map} or a \code{data.frame} with columns 'chr' (specifying the chromosome of the marker) and 'pos' (position of the marker within chromosome measured with genetic or physical distances)
}                                                                                  
  \item{dense}{
\code{logical}. Should density visualization for high-density genetic maps be used?
}
  \item{nMarker}{
\code{logical}. Print number of markers for each chromosome?
}
  \item{bw}{
\code{numeric}. Bandwith to use for  \code{dense=TRUE} to control the resolution.
}
  \item{centr}{
\code{numeric} vector. Position for the centromeres in the same order as chromosomes in \code{map}. If \code{"maize"}, centromere positions of maize are used.
}
  \item{\dots}{
further graphical arguments for function \code{plot}
}
}
\details{
In the low density plot, the unique position of markers are plotted as horizontal lines. In the  high-density plot, the distribution of the markers is visualized as a heatmap of density estimation together with a color key. In this case, the number of markers within a interval of equal bandwith \code{bw} is counted. The high density plot is typically useful if number of makers exceeds 200 on average per chromosome.
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
