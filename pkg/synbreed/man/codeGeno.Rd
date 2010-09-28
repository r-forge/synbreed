\name{codeGeno}
\alias{codeGeno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Code biallelic marker data and select markers  
}
\description{
This function could be used to code biallelic marker data in character format (e.g. alleles coded with "ade","gua", ... , or "A", "B", ...) into numeric format. Heterozygous genotypes could be code arbitrary, e.g. "ade/gua".
Imputing of missing values can be done with the population structure for the data (see details).  Selection of markers can be done according to the minor allele frequency and fraction of missing values.                                     
}
\usage{
codeGeno(data, impute = FALSE, popStruc = NULL, maf=NULL, nmiss=NULL, label.heter=NULL, replace.value=NULL)
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
A numeric or character vector of length \code{nrow(data)} giving the population id for each genotype with no missing values.
}
  \item{maf}{
A numeric. Markers with minor allele frequency < \code{maf} are discarded.
}
  \item{nmiss}{
A numeric. Markers with more than \code{nmiss} missing values are discarded.
}
  \item{label.heter}{
This is either a character to identify heterozygous genotypes or a function which evaluates if an element of the marker matrix is a heterozygous genotype.
}
  \item{replace.value}{
A numeric value assigned to missing values in data. Only specified if no impution is done 
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
A \code{data.frame} of the same dimension as \code{data} containing the recoded marker data.  The genotype which is homozygousfor the minor allele is coded as 2, the other homozygous is coded as 0 and heterozygous genotype is coded as 1.
}

\author{
Valentin Wimmer
}
\note{
Information about imputing is printed on the screen (approximative run time, number of imputations by family structure and number of random imputations).
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
  "ade",   NA,      NA,      "cyt",    "gua",   "thy",   "cyt",   "gua",  "ade",
  "ade",   NA,      "ade",   "cyt",   "gua",   "thy",   "gua",   "ade",  NA),ncol=9,byrow=TRUE)
  
colnames(snp9) <- paste("SNP",1:9,sep="")
snp9 <- data.frame(snp9)
pop <- c(rep("A",7),rep("B",2))

codeGeno(snp9,TRUE,pop)
# use only markers with less than 20% missing values
codeGeno(snp9,TRUE,pop,nmiss=0.2)

# example with heterozygous values
snp6 <- matrix(c(
  "gua",   "ade",   "gua",   "cyt",   "ade",   "thy",   
  "gua",   "ade",   "gua",   "cyt",   "ade",   "thy",  
  "ade/gua",   "ade",   "gua",   "cyt",   "ade",   "thy/gua",   
  "ade",   "ade",   "gua",   "cyt",   "ade",   "gua",  
  "ade/gua",   "ade",   "gua",   "cyt",   "ade",   NA,     
  "gua",   "ade",   "gua",   "cyt",   "gua",   "thy",
  "ade",   "ade",   NA,      "cyt",      NA,   "thy",
  "ade",   NA,      NA,      "cyt",   "gua",   "thy",
  "ade",   NA,      "ade",   "cyt",   "gua",   "thy"),ncol=6,byrow=TRUE)
colnames(snp6) <- paste("SNP",1:6,sep="")
snp6 <- data.frame(snp6)  
  
# two equivalent solutions
codeGeno(snp6,FALSE,label.heter="/")     
codeGeno(snp6,FALSE,label.heter=function(x) nchar(x)>3)   
                               
# example with maize data
data(maize)
# set 200 random choosen values to NA
marker <- maize.geno[,-1]
# coding without imputation
marker <- codeGeno(marker)
marker <- as.matrix(marker)
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
