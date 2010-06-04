\name{codeGeno}
\alias{codeGeno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to code biallelic marker data  
}
\description{
This function could be used to code biallelic marker data in character format (e.g. alleles coded with "ade","gua",...) into numeric format. The minor allele is coded as 2, the major allele is coded as 0.
Imputing of missing values can be done with the population structure for the data (see details).                                       
}
\usage{
codeGeno(data, impute = FALSE, popStruc = NULL)
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
A numeric or character vector of length \code{nrow(data)} giving the population id for each genotype. \code{popStruc} is needed if \code{impute=TRUE}. 
}
}
\details{
Imputing of missing values is done in the following way:
Suppose an observation \eqn{i} is missing (NA) for a marker \eqn{j} in population \eqn{k}. If marker \eqn{j} is fixed in population \eqn{k}, the imputed value will be the fixed allele. If marker \eqn{j} is segregating for the population \eqn{k}, 
the value is sampled from the allele distibution of population \eqn{k}. If there are no oberservations for a population no imputing is performed and a message is printed.
}
\value{
A \code{data.frame} of the same dimension as \code{data} containg the recoded marker data.
}

\author{
Valentin Wimmer
}
\note{

}

\seealso{
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
}
