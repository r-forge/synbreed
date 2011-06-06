\name{create.gpData}
\alias{create.gpData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{                                                                           
Create genomic prediction data object
}                  
\description{                                                                                                            
This function combines all raw data sources in a single, unified data object of class \code{gpData}. This is a \code{list} with elements for phenotypic, genotypic, pedigree and further covariate data. Moreover, the marker map is an element. Any element is optional
}                                                                                                              
\usage{                                                                                         
create.gpData(pheno = NULL, geno = NULL, map = NULL, 
pedigree = NULL, family = NULL, covar = NULL,
reorderMap = TRUE, map.unit = "cM")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pheno}{
\code{data.frame} with individuals organized in rows and traits organized in columns. Unique \code{rownames} identify individuals and unique \code{colnames} identify different traits. If no \code{rownames} are available, they are take from element \code{geno} (if available and if dimensions match). In case of multiple obsercations for each individual, use on column for every replication (such as locations or years).
}
  \item{geno}{
\code{matrix} with individuals organized in rows and markers organized in columns. Genotypes could be code arbitrarily. Missing values should be coded as \code{NA}. Unique \code{rownames} identify individuals and unique \code{colnames} markers. If no \code{rownames} are available, they are taken from element \code{pheno} (if available and if dimensions match).  If no \code{colnames} are used, the \code{rownames} of \code{map} are used if dimension matches.  
}                                                                                        
  \item{map}{                                                                                                
\code{data.frame} with one row for each marker and two columns. First columns gives the chromosome and second column the position  on the chromosome in centimorgan or the physical distance relative to the refernece sequence in basepairs. Unique \code{rownames} give the marker names which should match with marker names in \code{geno}. Note that order and number of markers must not be identical to order in \code{geno}. If this is the case, gaps in map are filled with \code{NA} and order is identical to \code{geno}.
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
  \item{reorderMap}{
\code{logical}. Should markers in \code{geno} and \code{map} be reordered by chromosome number and position within chromosome according to \code{map}? 
}
  \item{map.unit}{
Character. Unit of position in \code{map}, i.e. 'cM' for genetic distance or 'bp' for physical distance.
}
}
\details{
The class \code{gpData} is designed to provide a unified framework for data related to genomic prediction analysis. Every data source can be omitted. In this case, the corresponding argument must be \code{NULL}. By default (argument \code{reorderMap}), markers in \code{geno} are ordered by their position in \code{map}. Individuals are ordered in alphabetical order.

In an object of class \code{gpData} different individuals may occur in \code{pheno}, \code{geno} and \code{pedigree} are possible. In this case the \code{id} in \code{covar} comprises all individuals that either appear in \code{pheno}, \code{geno} and \code{pedigree}. Two additional columns in \code{covar} named \code{phenotyped} and \code{genotyped} identify individuals that appear in the corresponding object.   
}
\value{
Object of class \code{gpData} which is a list with the following items
\item{covar}{\code{list} with information on individuals}
\item{pheno}{\code{data.frame} with phenotypic data ordered by \code{rownames(pheno)}}
\item{geno}{\code{matrix} marker matrix containing genotypic data. Columns (marker) are in the same order as in \code{map}. Rows  ordered by \code{rownames(geno)}}
\item{pedigree}{object of class \code{pedigree}}
\item{map}{\code{data.frame} with columns 'chr' and 'pos' and markers sorted by 'pos' within 'chr'}
\item{info}{\code{list} with additional information on data (coding of data, unit in \code{map})}
}
\author{
Valentin Wimmer
}
\note{
In case of missing row names or column names in one item, information is substituted from other elements (assuming the same order of individuals/markers) and a warning is given. 
}
                                
\seealso{
\code{\link{codeGeno}}, \code{\link{summary.gpData}}, \code{\link{gpData2data.frame}}
}
\examples{
set.seed(123)
# 9 plants with 2 phenotypes 
n <- 9  # only for n > 6
pheno <- data.frame(Yield = rnorm(n,200,5), Height=rbeta(n,100,1))
rownames(pheno) <- sample(1:n)

# marker matrix
geno <- matrix(sample(c("AA","AB","BB",NA),size=n*12,replace=TRUE,
prob=c(0.6,0.2,0.1,0.1)),nrow=n)
rownames(geno) <-  sample(1:n)
colnames(geno) <- paste("M",1:12,sep="")

# genetic map
# one SNP is not mapped (M5)
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

