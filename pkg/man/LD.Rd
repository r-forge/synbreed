\name{pairwiseLD}
\alias{pairwiseLD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pairwise LD between markers
}
\description{
Compute pairwise LD between markers in  an object of class \code{gpData}. Return values is a \code{data.frame} with one row per marker pair or a \code{matrix} with all marker pairs. Additional, the euclidian distance between markers is taken from the marker map. 
}
\usage{
pairwiseLD(gpData, chr = NULL, type = c("data.frame", "matrix", "both"), LD.measure = "r2", rm.unmapped = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData}{
object of class \code{gpData} with elements \code{geno} and \code{map}
}                                 
  \item{chr}{
\code{numeric} vector or \code{"all"}. Return value is a list with one element for each chromosome in \code{chr}. If \code{chr="all"}, LD between markers across the whole genome is computed
}                                      
  \item{type}{
\code{character}. Specifies the type of return value.
}
  \item{LD.measure}{
\code{character}. Specifies the type of LD measure. Only \code{"r2"} possible at the moment.
}
  \item{rm.unmapped}{
\code{logical}. Remove markers with unknown postion from the LD analysis?
}
}
\details{
LD according to \eqn{r^2}{r2} is given by
\deqn{D= p_{AB} - p_Ap_B }{D = p(AB) - p(A)*p(B)} and
\deqn{r^2=\frac{D^2}{p_Ap_B+p_ap_b}}{r2=D2/p(A)p(B)p(a)p(b)}
where \eqn{p_{AB}}{p(AB)} is defined as the observed probability of haplotype \eqn{AB}, \eqn{p_A=1-p_a} and \eqn{p_B=1-p_b} the observed probabilities of alleles \eqn{A} and \eqn{B}. 
}
\value{
For \code{type="data.frame"} a list with one element for each chromosome. Each element is a \code{data.frame} with columns \code{marker1}, \code{marker2}, \code{r2} and \code{distance} for all \eqn{p(p-1)/2} marker pairs.

For \code{type="matrix"} a list with one element for each chromosome. Each element is a list of 2: a \eqn{p \times p}{p x p} \code{matrix} with pairwise LD and a  \eqn{p \times p}{p x p} \code{matrix} with pairwise distances.

For \code{type="both"} a list with both types described above is returned.
}
\references{
Hill WG, Robertson A (1968). Linkage Disequilibrium in Finite Populations. Theoretical
and Applied Genetics, 6(38), 226 - 231.
}
\author{
Valentin Wimmer
}

\seealso{
\code{\link{LDDist}}, \code{\link{LDMap}}
}
\examples{
data(maize)
maizeC <- codeGeno(maize)
head(pairwiseLD(maizeC,chr=1))
}
