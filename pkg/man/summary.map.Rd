\name{summaryGenMap}
\alias{summaryGenMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing marker map information
}
\description{
This function could be used to summarize information from a marker map and returns a \code{data.frame} with one row for each chromosome.
}
\usage{
summaryGenMap(map)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{map}{
\code{data.frame} with columns \code{chr} and \code{pos} or a \code{gpData} object with nonempty slot \code{map}
}
}
\details{
Summary statistics of differnces are based on euclidian distances between markers with available position, i.e. \code{pos!=NA}.
}
\value{
A \code{data.frame} with one row for each chromosome and columns
  \item{noM}{number of markers}
  \item{range}{range of positions, i.e. difference between first and last marker}
  \item{avDist}{avarage distance of markers}
  \item{maxDist}{maximum distance of markers}
  \item{minDist}{minimum distance of markers}
}
\author{
Valentin Wimmer
}
\seealso{
\code{\link{create.gpData}}
}
\examples{
data(maize)
summaryGenMap(maize)
}