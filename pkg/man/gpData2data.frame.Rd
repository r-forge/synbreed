\name{gpData2data.frame}
\alias{gpData2data.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Merge of phenotypic and genotypic data
}
\description{
Create a \code{data.frame} out of phenotypic and genotypic data in object of class \code{gpData} by merging datasets using the common id. The shared data set could either include individuals with phenotypes and genotypes (default) or additional unphenotyped or ungenotyped individuals. In the latter case, the missing observations are filled by \ode{NA}. Multiple observations for each individual in \code{pheno} could be reshape from "wide" into "long" format using a grouping variable. 
}
\usage{
gpData2data.frame(gpData,phenoNo=1,Rep=NULL,onlyPheno=!is.null(Rep),all.pheno=FALSE,all.geno=FALSE,...)
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
  \item{all.pheno}{
scalar logical. Include all individuals with genotypes in the \code{data.frame} and fill the genotypic data with \code{NA}.
}
  \item{all.geno}{
scalar logical. Include all individuals with phenotypes in the \code{data.frame} and fill the phenotypic data with \code{NA}.
}
  \item{...}{
further arguments to be used in function \code{reshape}. The argument \code{times} could be useful to rename levels of grouping variable.
}
}
\details{                                        
Argument \code{all.geno} can be used to predict the genetic value of individuals without phenotypic records using the \code{BLR} package. Here, the genetic value of individuals with \code{NA} as phenotype is predicted by the marker profile.


For multiple measures, phenotypic data in object \code{gpData} is arranged with replicates for each individual in seperate columns, i.e. i n so called "wide" format. With \code{gpData2data.frame} this coulb be reshaped to "long" format with multiple observations in one column. In this case, one column for the phenotype and 2 addtional columns for the \code{id} and the levels of the grouping variable are added.  
}
\value{
A \code{data.frame} with phenotypes in the first column(s) and marker matrix in subsequent columns.
}
\author{
Valentin Wimmer
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
# different subset of individuals in pheno and geno

# create 'gpData' object
gp <- create.gpData(pheno=pheno,geno=geno,map=map)
summary(gp)

# as data.frame with individuals with genotypes and phenotypes
gpData2data.frame(gp,1:2)
# as data.frame with all individuals with phenotypes
gpData2data.frame(gp,1:2,all.pheno=TRUE)
# as data.frame with all individuals with gotypes
gpData2data.frame(gp,1:2,all.geno=TRUE)



# example with repeated observations
set.seed(311)

# simulating genotypic and phenotypic data 
pheno <- data.frame(loc1=rnorm(10,1,2),loc2=rnorm(10,2,0.2),loc3=rbeta(10,2,4))
geno <- matrix(rep(c(1,0,2),10),nrow=10)
colnames(geno) <- c("M1","M2","M3")

# creating 'gpData' obbject
gp <- create.gpData(pheno=pheno,geno=geno)

# reshape of phenotypic data and merge of genotypic data,
# levels of grouping variable loc are named "a", "b" and "c"
gpData2data.frame(gp,onlyPheno=FALSE,
Rep=list(loc=c("loc1","loc2","loc3")),times=letters[1:3])
                            

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

