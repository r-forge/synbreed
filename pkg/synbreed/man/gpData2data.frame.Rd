\name{gpData2data.frame}
\alias{gpData2data.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Merge of phenotypic and genotypic data
}
\description{
Create a \code{data.frame} out of phenotypic and genotypic data in object of class \code{gpData} by merging datasets. Multiple observations for each individual in \code{pheno} could be reshape from "wide" into "long" format using a grouping variable. 
}
\usage{
gpData2data.frame(gpData,phenoNo=1,Rep=NULL,onlyPheno=!is.null(Rep),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData}{
object of class \code{gpData}
}
  \item{phenoNo}{
numeric or character. Which phenotypes should be in the \code{data.frame}? Only for \code{Rep=NULL}.
}
  \item{Rep}{
list. If multiple measures for each individual should be combined in the \code{data.frame} in "long" format, \code{Rep} is a named list with elements specifying variables in \code{pheno} that identify multiple observations, i.e. each element in \code{Rep} represents a level of the grouping variable. Name of list element is used as coloum name for the grouping variable. See \code{details} and \code{Examples}.
}
  \item{onlyPheno}{
scalar logical. Only return phenotypic data. Typically only useful, if \code{Rep} is used.
}
  \item{...}{
further arguments to be used in function \code{reshape}. The argument \code{times} could be useful to rename levels of grouping variable.
}
}
\details{
Resulting \code{data.frame} only includes individuals with phenotypes and genotypes. For multiple measures, phenotypic data in object \code{gpData} is arranged with replicates for each individual in seperate columns, i.e. i n so called "wide" format. With \code{gpData2data.frame} this coulb be reshaped to "long" format with multiple observations in one column. In this case, one column for the phenotype and 2 addtional columns for the \code{id} and the levels of the grouping variable are added.  
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
 \code{\link{create.gpData}}, \code{\link{reshape}}
}
\examples{
# example data with unrepeated observations
set.seed(311)

# simulating genotypic and phenotypic data 
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



# example with repeated observations
set.seed(311)

# simulating genotypic and phenotypic data 
pheno <- data.frame(loc1=rnorm(10,1,2),loc2=rnorm(10,2,0.2),loc3=rbeta(10,2,4))
geno <- matrix(rep(c(1,0,2),10),nrow=10)
colnames(geno) <- c("M1","M2","M3")

# creatin 'gpData' obbject
gp <- create.gpData(pheno=pheno,geno=geno)

# reshape of phenotypic data and merge of genotypic data,
# levels of grouping variable loc are named "a", "b" and "c"
gpData2data.frame(gp,onlyPheno=FALSE,
Rep=list(loc=c("loc1","loc2","loc3")),times=letters[1:3])
                            

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

