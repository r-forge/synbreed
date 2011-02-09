\name{[.relationshipMatrix}
\alias{[.relationshipMatrix}
\docType{data}
\title{
Extract or replace part of relationship matrix
}
\description{
Extract or replace part of an object of class \code{relationshipMatrix}.
}
\usage{
\method{[}{relationshipMatrix}(x,\dots)
}
\arguments{
  \item{x}{
object of class ''relationshipMatrix''
}
  \item{...}{
indices
}
}
\examples{
M<- matrix(c(2,2,0,0,0,0,0,2,0,0,0,0,2,2,0,2,2,0,0,0),nrow=5)
M
vR <- vanRaden(M)
vR[1:2,1:2]
}
