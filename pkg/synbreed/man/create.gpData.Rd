\name{create.gpData}
\alias{create.gpData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{                                                                           
Create genomic prediction data object
}                  
\description{                                                                                                            
This function takes a set of data and combines all in one object of class \code{gpData}. Main content of this object are the phenotypic data, the genotypic data (marker data), the genetic map and the pedigree. It is strongly recommended to use this class for further analysis with \code{synbreed} package. Note that not all types of data must be specified.
}                                                                                                              
\usage{
create.gpData(pheno = NULL, geno = NULL, map = NULL, 
pedigree = NULL, family = NULL, covar = NULL,
rm.unmapped = FALSE, map.unit = "cM")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pheno}{
\code{data.frame} with individuals in rows and one or more phenotypes in columns. \code{rownames} should identify genotypes and \code{colnames} different phenotypes. If no \code{rownames} are used, default rownames are used. In case of multiple measures for each individual an additional column for every replication should be added.
}
  \item{geno}{
\code{matrix} with individuals in rows and markers in columns. Genotypes could be code arbitrarily. Missing values should be coded as \code{NA}. \code{rownames} should identify genotypes and \code{colnames} markers. If no \code{rownames} are used, it is assumed that genotypes are in the same order as in \code{pheno}.  If no \code{colnames} are used, \code{rownames} of \code{map} are used if dimension matches. If no \code{map} is specified, default marker names are used. 
}                                                                                        
  \item{map}{                                                                                                
\code{data.frame} with one row for each marker and two columns. First columns should give the chromosome marker belongs to and second column the position (physical or genetic) on the chromosome. \code{rownames} should give the marker names which should match with marker names in \code{geno}. Note that order and number of markers must not be identical to order in \code{geno}. If this is the case, gaps in map are filled with \code{NA} and order is identical to \code{geno}.
}
  \item{pedigree}{
Object of class \code{pedigree}.
}
  \item{family}{                                        
\code{data.frame} assigning individuals to families with names of individuals in \code{rownames} This information could be used for replacing of missing values with function \code{codeGeno}.
}
  \item{covar}{
\code{data.frame} with further covariates for all individuals that either appear in \code{pheno}, \code{geno} or \code{pedigree$ID},  e.g. sex or age. \code{rownames} must be specified to identify individuals. Typically this element is not specified by the user.
}
  \item{rm.unmapped}{
Logical. Should unmapped markers be removed from genotypic data? (Not yet implemented)
}
  \item{map.unit}{
Character. Unit of position in \code{map}, i.e. 'cM' for genetic distance or 'bp' for physical distance.
}
}
\details{
The class \code{gpData} is designed to provide a unified framework for data related to genomic prediction analysis. Of course this class could be used for other purposes too. For a pedigree based prediction model only \code{pheno} and \code{pedigree} are necessary. In an analysis of LD only \code{geno} and \code{map} are needed.

In an object of class \code{gpData} different individuals in \code{pheno}, \code{geno} and \code{pedigree} are possible. In this case the \code{id} in \code{covar} comprises all individuals that either appear in \code{pheno}, \code{geno} and \code{pedigree}. Two addtional columns in \code{covar} named \code{phenotyped} and \code{genotyped} identify individuals that appear in the corresponding object.   
}
\value{
Object of class \code{gpData} which is a list with the following items
\item{covar}{\code{list} with information on individuals}
\item{pheno}{\code{data.frame} with phenotypic data ordered by \code{row.names(pheno)}}
\item{geno}{\code{matrix} marker matrix containing genotypic data. Columns (marker) are in the same order as in \code{map}. Rows  ordered by \code{row.names(geno)}}
\item{map}{\code{data.frame} with columns 'chr' and 'pos' and markers sorted by 'pos' within 'chr'}
\item{info}{\code{list} with addtional information on data (coding of data, unit in \code{map})}
}
\author{
Valentin Wimmer
}
\note{
In case of missing row names or column names in one item, information is taken from other items and a warning is given. 
}
                                
\seealso{
\code{\link{codeGeno}}, \code{\link{summary.gpData}}, \code{\link{gpData2data.frame}}
}
\examples{
set.seed(123)
# 9 plants with 2 phenotypes 
n <- 9  # only for n > 6
pheno <- data.frame(Yield = rnorm(n,100,5), leaf.damage=rbeta(n,1,1))
rownames(pheno) <- sample(1:n)

# marker matrix
geno <- matrix(sample(c("AA","AB","BB",NA),size=n*12,replace=TRUE,
prob=c(0.6,0.2,0.1,0.1)),nrow=n)
rownames(geno) <-  sample(1:n)
colnames(geno) <- paste("M",1:12,sep="")

# genetic map
# one snp is not mapped (M5)
map <- data.frame(chr=rep(1:3,each=4),pos=rep(1:12))
map  <- map[-5,]
rownames(map) <- paste("M",c(1:4,6:12),sep="") 

# simulate pedigree
ped <- simul.pedigree(3,c(3,3,n-6))

# combine in one object
gp <- create.gpData(pheno,geno,map,ped)
summary(gp)
}
\keyword{manip}

