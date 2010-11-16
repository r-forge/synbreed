\name{gpData2data.frame}
\alias{gpData2data.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Merge of phenotypic and genotypic data
}
\description{
Create a \code{data.frame} out of phenotypic and genotypic data in object of class \code{gpData} by merging datasets.
}
\usage{
gpData2data.frame(gpData, phenoNo = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData}{
object of class \code{gpData}
}
  \item{phenoNo}{
numeric or character. Which phenotypes should be in the \code{data.frame}
}
}
\details{
Resulting \code{data.frame} only includes individuals with phenotypes and genotypes.
}
\value{
A \code{data.frame} with phenotypes in the first column(s) and marker matrix in subsequent columns.
}
\author{
Valentin Wimmer
}
\note{
No implementation for replicated individuals atm
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{create.gpData}}
}
\examples{
# example data
set.seed(311)
pheno <- data.frame(Yield = rnorm(12,100,5),gdc=rbeta(12,1,1))
rownames(pheno) <- letters[4:15]
geno <- matrix(sample(c("A","A/B","B",NA),size=120,replace=TRUE,
prob=c(0.6,0.2,0.1,0.1)),nrow=10)
rownames(geno) <-  letters[1:10]
colnames(geno) <- paste("M",1:12,sep="")
# one snp is not mapped (M5)
map <- data.frame(chr=rep(1:3,each=4),pos=rep(1:12))
map  <- map[-5,]
rownames(map) <- paste("M",c(1:4,6:12),sep="")

# create 'gpData' object
gp <- create.gpData(pheno=pheno,geno=geno,map=map)
summary(gp)

# as data.frame
gpData2data.frame(gp,1:2)
                            

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

