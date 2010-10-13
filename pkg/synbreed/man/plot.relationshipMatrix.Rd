\name{plot.relationshipMatrix}
\alias{plot.relationshipMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Heatmap for relationship Matrix
}
\description{
Visualistion for objects of class \code{relationshipMatrix}.  
}
\usage{
\method{plot}{relationshipMatrix}(x, \ldots) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class \code{relationshipMatrix}
}
  \item{\dots}{
further graphical arguments
}
}

\author{
Valentin Wimmer
}

\examples{
ped <- simul.pedigree(gener=4,7)
A <- kinshipPlants(ped)
plot(A)
}

\keyword{hplot}


