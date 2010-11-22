\name{maize}
\alias{maize}
\docType{data}
\title{
Simulated maize data
}
\description{
This is a simulated dataset from breeding sheme for maize. 1250 doubled haploid (DH) lines were genotyped with 1117 ploymorphic SNP markers and phenotyped in a testcross with a single tester for one quantitative trait. Markers are spread along all 10 chromosomes of maize. Pedigree information is available up to 15 generations. The 1250 lines belong to 25 full sib families with 50 individuals in each family. For simulation 1000 QTL with equal and additive (no dominance or epistasis) effects were used. True breeding values for individuals were simulated according to 
\deqn{TBV=\sum_{i=1}^{1000} QTL_i + \epsilon_i,}{TBV=\sum QTL + e,}
where \eqn{\epsilon \sim N(0,I \sigma^2)}{e = N(0,sigma2)}. The value for \eqn{\sigma^2}{sigma2} was choosen in a way that a given plot heritability of \eqn{h^2=0.197}{h2=0.197} is realized. Note that true breeding values for 1250 phenotypes lines are stored as \code{tbv} in \code{covar} of \code{gpData} object. Reported phenotypic values are testcross means of lines from 3 replications.
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
                                  