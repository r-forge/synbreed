\name{plot.genMap}
\alias{plot.genMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot genetic map
}
\description{
A function to visualize high-density genetic maps.
}
\usage{
plot.genMap(map, dense = FALSE, nMarker = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{map}{
object of class \code{gpData} with object \code{map} or a \code{data.frame} with columns 'chr' and 'pos'
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
Plots the genetic map for each chromosome. Chromosome is displayed from first to last marker.
}

\author{
Valentin Wimmer
}


\seealso{
\code{\link{create.gpData}}
}
\examples{
# low density plot
plot.genMap(maize)

# high density plot
plot.genMap(mice,TRUE,FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
