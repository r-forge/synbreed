\name{rogers}
\alias{rogers}
\title{
Genomic relatedness based on marker similarity
}
\description{
This function computes the genomic relationship matrix with marker information from genotypic data with \eqn{n} genotypes and \eqn{m} biallelic markers. Kinship is reciprocal to Rogers's distance
}
\usage{
rogers(marker, correction = c("none","Hayes", "Melchinger"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{marker}{
A \code{matrix} or \code{data.frame} with genotypes in rows and markers in columns. Alleles should be coded 0/2 for homozygous loci and 1 for heterozygous, i.e. \code{marker} could be the result  of a call of \code{codeGeno}.
}
  \item{correction}{
Character. Use \code{"Hayes"} for transformation of rogers distance into kinship coefficient by formula of Hayes and Goddard (2008) or use \code{"Melchinger"} for transformation according to Melchinger (1991) or any abbreviations of these names. If \code{"none"} (default) no correction is performed.
}
}
\details{
Relationship coefficient is computed according to Roger's distance 
\deqn{f=1-d}
where \eqn{d} is Roger's distance which is computed as
\deqn{d=\frac{1}{m}\sum_{i=1}^m \sqrt{1/2 \sum_{j=1}^{n_i}(p_{ij}-q_{ij})^2}}{d=1/m sum(sqrt(\sum(pij-qij)2))}
where \eqn{m} is the number of markers and \eqn{n_i}{ni} the number of alleles for marker \eqn{i}. \eqn{p_{ij}}{pij} and \eqn{q_{ij}}{qij} are the allele frequencies of allele \eqn{j} for marker \eqn{i}.
Kinship based on Rogers distance is the simple matching coefficient which is \eqn{f=1-d} (correction = \code{"none"}). 

Transformation to kinship coefficient according to Hayes and Goddard (2008)
\deqn{f = \frac{s-s_{min}}{1-s_{min}}}{f=s-smin/(1-smin)}
where \eqn{s=1-d}.

Transformation according to Melchinger (1991) is
\deqn{f=1-\frac{d}{\bar{d}}}{f=1-d/dquer}

In both cases the relationship \eqn{2f} which is twice the kinship returned by the function.
}
\value{
An object of class \code{relationshipMatrix} which is a symmetric \eqn{n \times n}{n x n} matrix containing the relationship coefficients (two times the kinship coefficient).
}
\references{
Rogers, J., 1972 Measures of genetic similarity and genetic distance. In Studies in genetics
VII, volume 7213. Univ. of Texas, Austin, 

Melchinger, A.E., M.M. Messmer, M. Lee, W.L. Woodman, and K.R. Lamkey: Diversity and relationships among U.S. maize inbreds revealed by restriction fragment length polymorphisms. Crop Sci. 31 

Hayes, B. J., and M. E. Goddard. 2008. Technical note: Prediction of breeding values using marker derived relationship matrices.
J. Anim. Sci. 86
}
\author{
Valentin Wimmer
}                    
\seealso{
\code{\link{vanRaden}}, \code{\link{codeGeno}}
}
\examples{
M<- matrix(c(2,2,0,0,0,0,0,2,0,0,0,0,2,2,0,2,2,0,0,0),nrow=5)
colnames(M) <- paste("M",1:4,sep="")
rownames(M) <- 1:5
M
rogers(M,correction="Hayes")
rogers(M,correction="Melchinger")
}
