\name{summary.relationshipMatrix}
\alias{summary.relationshipMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing relationship matrices
}
\description{
\code{summary} method for class ''relationshipMatrix''
}
\usage{
\method{summary}{relationshipMatrix}(relationshipMatrix) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{relationshipMatrix}{
object of class \code{relationshipMatrix}
}

}

\author{
Valentin Wimmer
}
\note{
}

\examples{
ped <- simul.pedigree(gener=4,7)
A <- kinshipPlants(ped)
summary(A)
}

\keyword{relationshipMatrix}

