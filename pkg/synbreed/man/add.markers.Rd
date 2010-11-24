\name{add.markers}
\alias{add.markers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adding new markers to 'gpData' object
}
\description{
This function adds new markers to genotypic data and updates genetic map.
}
\usage{
add.markers(gpData, geno, map)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData}{
object of class \code{gpData}
}
  \item{geno}{
\code{matrix} with new columns
}
  \item{map}{
\code{data.frame} with columns 'chr' and 'pos' and new rows for new markers 
}
}
\details{
Make sure that \code{rownames} in \code{geno} match with \code{rownames} genotypic data in \code{gpData} object. 

}
\value{
object of class \code{gpData} with new markers
}

\author{
Valentin Wimmer
}

\seealso{
\code{\link{add.individuals}}, \code{\link{create.gpData}}
}
\examples{
# creating gpData object
# phenotypic data
pheno <- data.frame(Yield = rnorm(10,100,5), leaf.damage=rbeta(10,1,1))
rownames(pheno) <- 1:10
# genotypic data
geno <- matrix(sample(c(1,0,2,NA),size=120,replace=TRUE,
prob=c(0.6,0.2,0.1,0.1)),nrow=10)
rownames(geno) <-  1:10
colnames(geno) <- paste("M",1:12,sep="")
# genetic map
map <- data.frame(chr=rep(1:3,each=4),pos=rep(1:12))
# as gpData object
gp <- create.gpData(pheno,geno,map)


# new data
geno2 <- matrix(c(0,0,1,1,1,2,2,1,1,2,1,2,0,2,1,1,1,2,2,2),ncol=2)
rownames(geno) <- 1:10
map2 <- data.frame(pos=c(0.3,5),chr=c(1,2))
rownames(map2) <- colnames(geno2) <- c("M13","M14")

# adding new markers
gp2 <- add.markers(gp,geno2,map2)
summary(gp2)
                                                                                       
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
