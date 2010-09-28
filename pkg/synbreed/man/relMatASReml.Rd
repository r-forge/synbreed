\name{relMatASReml}
\alias{relMatASReml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inverse relationship matrix for ASReml
}
\description{
This function computes the inverse or generalized inverse relationship matrix and returns it in a format suitable for ASReml.
}                                                   
\usage{
relMatASReml(relationshipMatrix, ginv = FALSE, file = NULL, digits = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{relationshipMatrix}{
An object of class \code{relationshipMatrix} or any other relationship matrix 
}
  \item{ginv}{
logical. Should the generalized inverse be used?
}
  \item{file}{
path for writing the result. 
}
  \item{digits}{
numeric. The result is rounded to \code{digits}.
}
}
\details{
The results is ordered in the format for ASReml inverse numerator relationship matrizes (first ordered by row and then by column).
}
\value{
Result is writen to \code{file} or returned if \code{file=NULL}
}

\author{
Valentin Wimmer
}

\examples{
ped.raw <- matrix(c(1,  0,   0,   1, 
                 2,  0,   0,   1, 
                 3,  0,   0,   1,
                 4,  0,   0,   1,
                 5,  1,   2,   2,
                 6,  1,   3,   2, 
                 7,  1,   2,   2, 
                 8,  1,   2,   2,
                 9,  2,   2,   2,
                 10, 5,   6,   3,  
                 11, 5,   1,   3),byrow=TRUE,ncol=4)

ped <- create.pedigree(ped.raw[,1],ped.raw[,2],ped.raw[,3],ped.raw[,4])

A <- kinshipPlants(ped)
relMatASReml(A)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{relationshipMatrix}

