\name{vanRaden}
\alias{vanRaden}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genomic relationship matrix according to vanRaden
}
\description{
This function computes the genomic relationship matrix with marker information from genotypic data with \eqn{n} genotypes and \eqn{p} markers. 
}                                
\usage{
vanRaden(marker)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{marker}{
 A \code{matrix} or \code{data.frame} with genotypes in rows and markers in columns. Alleles should be coded 0/2 for homozygous and 0/1/2 for heterozygous loci, i.e. \code{marker} could be the result  of a call of \code{codeGeno}.
}                                                                                                                                                                                                 
}
\details{
The computation is the same for homozygous and heterozygous loci. The formula for the genomic relationship matrix according to vanRaden is
\deqn{G = \frac{ZZ'}{2\sum p_i(1-p_i)}}{ZZ'/(2\sum pi(1-pi))} where \eqn{Z=M-P} and \eqn{M} is the marker matrix and \eqn{P} contains the allele frequencies multiplied by 2. \eqn{p_i} is the allele frequency of marker \eqn{i}. 
}
\value{
An object of class \code{relationshipMatrix} which is a symmetric \eqn{nxn} matrix containing the relationship coefficients.
}
\references{
%% ~put references to the literature/web site here ~
}                                                      
\author{
Valentin Wimmer and Theresa Albrecht
}


\seealso{
\code{\link{codeGeno}}
}
\examples{
M<- matrix(c(2,2,0,0,0,0,0,2,0,0,0,0,2,2,0,2,2,0,0,0),nrow=5)
M
vanRaden(M)
}
\keyword{relationshipMatrix}

