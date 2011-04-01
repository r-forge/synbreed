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
further graphical arguments passed to function \code{levelplot} in package \code{lattice}. To create equal colorkeys for two heatmaps, use \code{at=seq(from,to,length=9)}. 
}
}

\author{
Valentin Wimmer
}

\examples{
# small pedigree
ped <- simul.pedigree(gener=4,7)
gp <- create.gpData(pedigree=ped)
A <- kin(gp,ret="add")
plot(A)

# big pedigree
data(maize)
A <- kin(maize,ret="add")
U <- kin(codeGeno(maize),ret="realized")/4
# equal colorkeys
plot(A,at=seq(0,2,length=9))
plot(U,at=seq(0,2,length=9))
}

\keyword{hplot}


