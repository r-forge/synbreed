\name{discard.markers}
\alias{discard.markers}
\alias{discard.individuals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subsetting class 'gpData'
}
\description{
Both functions could be used to get subsets of class \code{gpData}. To discard columns in \code{geno}, use function \code{discard.markers} to discard markers from \code{gpData}. Note that markers will also be removed from \code{map}, if available. To discard rows of \code{geno}, use function \code{discard.individuals} to discard individuals (genotypes) from \code{gpData}. Note that individuals will also be removed from \code{covar}, \code{pheno} and \code{pedigree}. 
}                      
\usage{
discard.markers(gpData, which)
discard.individuals(gpData, which)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData}{
object of class \code{gpData}
}
  \item{which}{
character vector either identifying the \code{colnames} of markers in \code{geno} to discard (in \code{discard.markers}) or the \code{rownames} of individuals to discard (in \code{discard.individuals}).
}
}
\value{
Object of class \code{gpData} 
}
\author{
Valentin Wimmer
}

\seealso{
\code{\link{create.gpData}}
}
\examples{
# example data
set.seed(311)
pheno <- data.frame(Yield = rnorm(10,100,5),gdc=rbeta(10,1,1))
rownames(pheno) <- letters[1:10]
geno <- matrix(sample(c("A","A/B","B",NA),size=120,replace=TRUE,
prob=c(0.6,0.2,0.1,0.1)),nrow=10)
rownames(geno) <-  letters[1:10]
colnames(geno) <- paste("M",1:12,sep="")
# one snp is not mapped (M5)
map <- data.frame(chr=rep(1:3,each=4),pos=rep(1:12))
map  <- map[-5,]
rownames(map) <- paste("M",c(1:4,6:12),sep="")
gp <- create.gpData(pheno=pheno,geno=geno,map=map)
summary(gp)

# remove unmapped SNP from gpData
gp2 <- discard.markers(gp,"M5")
summary(gp2)

# discard genotypes with missing values
gp3 <- discard.individuals(gp,names(which(rowSums(is.na(gp$geno))>0)))
summary(gp3)

 
}
\keyword{manip}
