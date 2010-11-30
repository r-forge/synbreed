\name{kinshipPlants}
\alias{kinshipPlants}
\title{
Compute the numerator relationship matrix for plant pedigree}
\description{
This function could be used to compute the numerator relationship matrix from a kinship structure. This function could be used for plant pedigree because 
inbreeding is allowed in pedigree structure.
}
\usage{
kinshipPlants(ped, inbreed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped}{
An object of class \code{pedigree} .
}
  \item{inbreed}{
A logical. TRUE if genotype is fully homozygous and FALSE otherwise. If \code{NULL}, a genotype with Par1 == Par2 is seen as inbreed.
}

}
\details{
Note that this function allows inbreeds where the mother is identical to the father.
}

\value{
The numerator relationship matrix (A-matrix in the aninmal model) derived by the pedigree structure.
}

\author{
Hans-Juergen Auinger and Valentin Wimmer
}

\note{ 
This function is obsolete and only available due to downward compatibility. Better use the new function \code{kinship} with more possibilities.  
}

\seealso{
\code{\link{plot.relationshipMatrix}} ,\code{\link{kinship}}
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
}
