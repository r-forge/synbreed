\name{plot.relationshipMatrix}
\alias{plot.relationshipMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Heatmap for relationship Matrix
}
\description{
Visualization for objects of class \code{relationshipMatrix}.  
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
further graphical arguments passed to function \code{levelplot} in package \code{lattice}
}
}

\author{
Valentin Wimmer
}

\examples{
ped <- simul.pedigree(gener=4,7)
A <- kin(ped,ret="add")
plot(A)
}

\keyword{hplot}


