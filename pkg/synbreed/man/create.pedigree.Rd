\name{create.pedigree}
\alias{create.pedigree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create pedigree object
}
\description{
This function can be used to create a \code{pedigree} object. This function allows inlcuding inbreeds in the pedigree.
}
\usage{
create.pedigree(ID, Par1, Par2, gener=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ID}{
vector of unique IDs identifying e.g. each genotype. Should be numeric.
}
  \item{Par1}{
vector of IDs identifying parent 1 (with animals: father) 
}
  \item{Par2}{
vector of IDs identifying parent 2 (with animals: mother)
}
  \item{gener}{
vector identifying the generation. If \code{NULL} gener will be 0 for unknown parents and max(gener(Par1),gener(Par2))+1 otherwise.
}
}
\details{
Missing values for pedigree should be coded with 0 for numeric ID or \code{NA} for character ID. 
}
\value{
An object of class \code{pedigree} with the same contents as above. Note that \code{gener} starts from 0 and pedigree is sorted by generation.
}                                                                             
\author{
Valentin Wimmer
}

\seealso{
\code{\link{plot.pedigree}}
}
\examples{
# example with given generation
ped.raw <- matrix(c(1,  0,   0,   1, 
                 2,  0,   0,   1,
                 3,  0,   0,   1,
                 4,  0,   0,   1,
                5,  1,   2,   2, 
                 6,  1,   3,   2,
                 7,  1,   2,   2,
                 8,  1,   2,   2) ,
                 byrow = TRUE, ncol = 4)
ped <- create.pedigree(ped.raw[,1],ped.raw[,2],ped.raw[,3],ped.raw[,4])
ped

# example without given generation
ped.raw2 <- matrix(c(1,  0,   0,   
                 2,  0,   0,  
                 3,  0,   0,  
                 4,  0,   0,  
                5,  1,   2,   
                 6,  1,   3,  
                 7,  1,   2,  
                 8,  1,   2) ,
                 byrow = TRUE, ncol = 3)
ped2 <- create.pedigree(ped.raw2[,1],ped.raw2[,2],ped.raw2[,3])
ped2
}
\keyword{pedigree}