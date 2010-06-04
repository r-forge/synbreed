\name{kinshipPlants}
\alias{kinshipPlants}
\title{
Compute the numerator relationship matrix for plant pedigree}
\description{
This function could be used to compute the numerator relationship matrix from a kinship structure. This function could be used for plant pedigree because 
inbreeding is allowed in pedigree structure.
}
\usage{
kinshipPlants(id, inbreed, par1.id, par2.id)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{
A numeric vector containing the id's for each genotype
}
  \item{inbreed}{
A logical. TRUE if genotype is fully homozygous and FALSE otherwise
}
  \item{par1.id}{
A numeric vector with the id's for the first parent (mother)
}
  \item{par2.id}{
A numeric vector with the id's for the second parent (father)
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

\examples{
pedigree <- matrix(c(1,  0,   0,   1, 0,
                 2,  0,   0,   1, 0,
                 3,  0,   0,   1, 0,
                 4,  0,   0,   1, 0,
                 5,  1,   2,   2, 0,
                 6,  1,   3,   2, 0, 
                 7,  1,   2,   2, 0,
                 8,  1,   2,   2, 0,
                 9,  2,   2,   2, 1,
                 10, 5,   6,   3, 0, 
                 11, 5,   1,   3, 0),byrow=TRUE,ncol=5)
colnames(pedigree) <- c("ID", "P1", "P2", "gener", "line")

A <- kinshipPlants(pedigree[, 1], pedigree[, 5], pedigree[, 2], pedigree[,3])
                 
                  

}