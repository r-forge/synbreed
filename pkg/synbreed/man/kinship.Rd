\name{kinship}
\alias{kinship}
\title{
Kinship based on pedigree
}
\description{
This function computes the relationship matrix between animals or plants based on pedigree information. It is possible to account for doubled haploid (DH) lines as their inbreeding coefficient is set 1.

}
\usage{
kinship(ped, DH = NULL, ret = c("add", "kin", "dom", "gam"))
}
\arguments{
  \item{ped}{
object of class \code{pedigree} with \eqn{n} individuals or object of class \code{gpData}
}
  \item{DH}{
vector with logicals of length \eqn{n}. \code{TRUE} or 1 if individual is a DH line and \code{FALSE} or 0 otherwise. 
%If \code{ped} is of class \code{gpData}, \code{DH} is a character giving the name of a variable in \code{covar} that contains this information for every individual.
}
  \item{ret}{
character. The type of relationship matrix to be returned. See 'Details'.
}
}
\details{
In the first step, the gametic relationship is build. The gametic relationship is of order 2\eqn{n} as each individual A has two alleles (\eqn{A1} and \eqn{A2}). One of these alleles is passed to a progeny. The gametic relationship
is defined as the matrix of probabilities that two genes are identical by descent (IBD).  Note that the diagonal elements of the gametic relationship matrix are 1. The off-diagonals of individuals with unknown pedigree are 0. If \code{ret="gam"} is specified, the gametic relationship matrix constructed by pedigree is returned.

The gametic relationship matrix can be used to set up other types of relationship matrices. If \code{ret="add"}, the additive numerator relationship matrix is returned. The additive relationship of individuals A (alleles \eqn{A1,A2}) and B (alleles \eqn{B1,B2}) is given by  the entries of the gametic relationship matrix
\deqn{0.5\cdot \left[(A1,B1) + (A1,B2) + (A2,B1) + (A2,B2)\right],}{0.5*[(A1,B1) + (A1,B2) + (A2,B1) + (A2,B2)],}
where \eqn{(A1,B1)} denotes the element [A1,B1] in the gametic relationship matrix. If \code{ret="kin"}, the kinship matrix is returned which is half of the additive relationship matrix.

If \code{ret="dom"}, the dominance relationship matrix is returned. The dominance relationship matrix between individuals A (\eqn{A1,A2}) and B (\eqn{B1,B2}) in case of no inbreeding is given by 
\deqn{\left[(A1,B1) \cdot (A2,B2) + (A1,B2) \cdot (A2,B1)\right],}{[(A1,B1) * (A2,B2) + (A1,B2) * (A2,B1)],} 
where \eqn{(A1,C1)} denotes the element [A1,C1] in the gametic relationship matrix. 
}
\value{
An object of class "relationshipMatrix".
}
\author{
Valentin Wimmer
}
\note{
This is a new version of function \code{kinshipPlants}.
}

\seealso{
\code{\link{plot.relationshipMatrix}} 
}
\examples{

# Animal Example
ped1 <- create.pedigree(c("A","B","C","D","E"),c(0,0,"A","A","D"),c(0,0,"B","C","B"))
kinship(ped1,DH=rep(0,5),"add")

# Plants Example with DH lines
id <- c("A","T","DHa","DHb","TCa","TCb")
par1 <- c(0,0,"A","A","DHa","DHb")
par2 <- c(0,0,"A","A","T","T")
ped2 <- create.pedigree(id,par1,par2)
ped2
# additive relationship
kinship(ped2,DH=c(0,1,1,1,0,0),"add")
# dominance relationship
kinship(ped2,DH=c(0,1,1,1,0,0),"dom")

}
