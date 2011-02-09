\name{write.relationshipMatrix}
\Rdversion{1.1}
\alias{write.relationshipMatrix}
\title{
Writing relationshipMatrix in table format
}
\description{
This function could be used to write an object of class "relationshipMatrix" in table format to file in a way that it could be used by other software, i.e. WOMBAT or ASReml.  The table has three columns, the row, the column and the entry of the (inverse) relationshipMatrix.
}
\usage{
write.relationshipMatrix(relationshipMatrix, file = NULL, 
sorting = "WOMBAT", type=c("ginv", "inv", "none"), digits = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{relationshipMatrix}{
Object of class "relationshipMatrix" 
}
  \item{file}{
Path were the output should be writen . If \code{NULL} the result is returned in R.
}
  \item{sorting}{
Type of sorting. Use "WOMBAT" for 'row-wise' sorting of the table and "ASReml" for 'column-wise' sorting.
}
  \item{type}{
A character string indicating which form of \code{relationshipMatrix} should be returned. One of "ginv" (Moore-Penrose generalized inverse), "inv" (inverse), or "none" (no inverse). 
}
  \item{digits}{
Numeric.  The result is rounded to \code{digits}.
}
}
\details{
Note that "WOMBAT" only uses the generalized inverse relationship matrix and expects a file with the name "ranef.gin", where 'ranef' is the name of the random effect with option 'GIN' in the 'MODEL' part of the parameter file.
For ASREML, either the relationship could be saved as "*.grm" or its generalized inverse as "*.giv".
  
}
\references{
Meyer, K. (2006) WOMBAT - A tool for mixed model analyses in quantitative genetics by REML, J. Zhejinag Uni SCIENCE B 8: 815-821.

Gilmour, A., Cullis B., Welham S., and Thompson R. (2000) ASREML. program user manual. NSW
Agriculture, Orange Agricultural Institute, Forest Road, Orange, Australia .
}
\author{
Valentin Wimmer
}

\examples{
ped.raw <- matrix(c(1,  0,   0,   
                 2,  0,   0,   
                 3,  0,   0,  
                 4,  0,   0,  
                 5,  1,   2,   
                 6,  1,   3,   
                 7,  1,   2,   
                 8,  2,   2 ),byrow=TRUE,ncol=3)

ped <- create.pedigree(ped.raw[,1],ped.raw[,2],ped.raw[,3])

A <- kinship(ped,ret="add")
write.relationshipMatrix(A,type="ginv")
}
\keyword{IO}
