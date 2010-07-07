\name{codeGeno}
\alias{codeGeno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to code biallelic marker data  
}
\description{
This function could be used to code biallelic marker data in character format (e.g. alleles coded with "ade","gua",...) into numeric format. The minor allele is coded as 2, the major allele is coded as 0.
Imputing of missing values can be done with the population structure for the data (see details).  Selection of markers can be done according to the minor allele frequency and fraction of missing values.                                     
}
\usage{
codeGeno(data, impute = FALSE, popStruc = NULL, maf=NULL, nmiss=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A \code{data.frame} or a \code{matrix} containing the  marker data. Each column corresponds to a biallelic marker, each row corresponds to a genotype. Coding of alleles is arbitrary, missing values have to be coded as \code{NA}.
}
  \item{impute}{
A logical specifing whether imputing of missing values should be done.
}
  \item{popStruc}{
A numeric or character vector of length \code{nrow(data)} giving the population id for each genotype. 
}
  \item{maf}{
A numeric. Markers with minor allele frequency < \code{maf} are discarded.
}
  \item{nmiss}{
A numeric. Markers with more than \code{nmiss} missing values are discarded.
}
}
\details{
Imputing of missing values is done in the following way:
\describe{
\item{with population structure}{
Suppose an observation \eqn{i} is missing (NA) for a marker \eqn{j} in population \eqn{k}. If marker \eqn{j} is fixed in population \eqn{k}, the imputed value will be the fixed allele. If marker \eqn{j} is segregating for the population \eqn{k}, 
the value is 0 with probability 0.5 and 1 with probability 0.5.}
\item{without population structure}{The missing values for a marker \eqn{j} are sampeled from the allele distribution of marker \eqn{j}.} 
}
}
\value{
A \code{data.frame} of the same dimension as \code{data} containg the recoded marker data.
}

\author{
Valentin Wimmer
}
\note{
If there are no oberservations for a population no imputing is performed and a error message is printed but computation is continued.
}

\examples{
snp9 <- matrix(c(
  "gua",   "ade",   "gua",   "cyt",   "ade",   "thy",   "cyt",   "gua",  NA,
  "gua",   "ade",   "gua",   "cyt",   "ade",   "thy",   "gua",   "ade",  NA,
  "ade",   "ade",   "gua",   "cyt",   "ade",   "thy",   "gua",   "ade",  NA,
  "ade",   "ade",   "gua",   "cyt",   "ade",   "gua",   "cyt",   "gua",  NA,
  "ade",   "ade",   "gua",   "cyt",   "ade",   NA,      "cyt",   "ade",  NA,
  "gua",   "ade",   "gua",   "cyt",   "gua",   "thy",   "cyt",   "gua",  NA,
  "ade",   "ade",   NA,      "cyt",      NA,   "thy",   "gua",   "ade",  "ade",
  "ade",   NA,      NA,      "cyt",   "gua",   "thy",   "cyt",   "gua",  "ade",
  "ade",   NA,      "ade",   "cyt",   "gua",   "thy",   "gua",   "ade",  NA),ncol=9,byrow=TRUE)
  
colnames(snp9) <- paste("SNP",1:9,sep="")
snp9 <- data.frame(snp9)
pop <- c(rep("A",7),rep("B",2))

codeGeno(snp9,TRUE,pop)
# use only markers with less than 20% missing values
codeGeno(snp9,TRUE,pop,nmiss=0.2)
                               
# example with maize data
data(maize)
# set 200 random choosen values to NA
marker <- maize.geno[,-1]
# coding without imputation
marker <- codeGeno(marker)
set.seed(123)
ind1 <- sample(1:nrow(marker),200)
ind2 <- sample(1:ncol(marker),200)
original <- marker[cbind(ind1,ind2)]

marker[cbind(ind1,ind2)] <- NA
# 25 populations
pop <- rep(1:25,each=50)
marker1 <- codeGeno(marker,impute=TRUE,pop)
# compare in crosstable
(t1 <- table(original,marker1[cbind(ind1,ind2)]) )
# sum of correct replacements
sum(diag(t1))/sum(t1)

# compare with wrong familiy information
pop2 <- sample(1:25,size=1250,replace=TRUE)
marker2 <- codeGeno(marker,impute=TRUE,pop2)
(t2 <- table(original,marker2[cbind(ind1,ind2)]) )
# sum of correct replacements
sum(diag(t2))/sum(t2)
}
\keyword{genotypic data}
