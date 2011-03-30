\name{add.individuals}
\alias{add.individuals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adding new individuals to objects of class gpData
}
\description{
This function adds new data for individuals to an object of class \code{gpData} by adding new phenotypes, genotypes and pedigree.
}
\usage{
add.individuals(gpData, pheno = NULL, geno = NULL,
 pedigree = NULL, covar = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData}{
object of class \code{gpData}
}
  \item{pheno}{
\code{data.frame} with new rows for phenotypes
}
  \item{geno}{
\code{matrix} with new rows for genotypic data
}
  \item{pedigree}{
\code{data.frame} with new rows for pedigree data
}      
  \item{covar}{
\code{data.frame} with new rows for \code{covar} information
}
}
\details{
Make sure that \code{colnames} in \code{geno}, \code{pheno} and \code{pedigree} match with existing names in \code{gpData} object. 
}
\value{
object of class \code{gpData} with new individuals
}
\author{
Valentin Wimmer
}

\seealso{
\code{\link{add.markers}}, \code{\link{discard.individuals}}
}
\examples{
# adding one new DH line to maize data
data(maize)                                          
newDHpheno <- data.frame(Trait=1000,row.names="newDH")
# simulating genotypic data
newDHgeno <- matrix(sample(c(0,1),ncol(maize$geno),replace=TRUE),nrow=1)
rownames(newDHgeno) <- "newDH"
# new pedigree
newDHpedigree <- data.frame(ID="newDH",Par1=0,Par2=0,gener=0)
# new covar information
newDHcovar <- data.frame(family=NA,DH=1,tbv=1000,row.names="newDH")

maize2 <- add.individuals(maize,newDHpheno,newDHgeno,newDHpedigree,newDHcovar)  
summary(maize2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
