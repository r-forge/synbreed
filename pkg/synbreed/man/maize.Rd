\name{maize}
\alias{maize}
\docType{data}
\title{
Simulated maize data
}
\description{
This is a simulated dataset of a maize breeding sheme. Data comprises 1250 doubled haploid (DH) lines that were genotyped with 1117 polymorphic SNP markers and phenotyped in a testcross with a single tester for one quantitative trait. Markers are distributed along all 10 chromosomes of maize. Pedigree information starts with basis population and is available up to 15 generations. The 1250 lines belong to 25 full sib families with 50 individuals in each family. In the simulation of true breeding values (TBV), 1000 biallelic quantitative trait loci (QTL) with equal and additive (no dominance or epistasis) effects were generated. True breeding values for individuals were calculated according to 
\deqn{tbv=\sum_{k=1}^{1000} QTL_k}{TBV=\sum QTL(k)}
where \eqn{QTL_k}{QTL(k)} is the effect of the \eqn{k}-th QTL. Phenotypic values were simulated according to
\deqn{y_i=tbv_i + \epsilon_i}{trait=tbv+e}
where \eqn{\epsilon \sim N(0,I \sigma^2)}{e = N(0,sigma2)}. The value for \eqn{\sigma^2}{sigma2} was choosen in a way that a given plot heritability of \eqn{h^2=0.197}{h2=0.197} is realized. Note that true breeding values for 1250 phenotyped lines are stored as \code{tbv} in \code{covar} of \code{gpData} object. Reported phenotypic values of lines are adjusted values testcross means evaluated in 3 replications.
}                                                                                                                                                                                                                 
\usage{data(maize)}                                                                                
\format{
Object of class \code{gpData}
}
                                
\examples{
data(maize)
summary(maize)                                       
}                                  
\keyword{datasets}
                                  