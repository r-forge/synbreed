\name{pairwiseLD}
\alias{pairwiseLD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pairwise LD between markers
}
\description{
Compute pairwise LD between markers in  an object of class \code{gpData}. For the general case, a gateway to the software PLINK (Purcell et al. 2007) is established to estimate the LD. A within-R solution is available for marker data with only 2 genotypes, i.e. homozgous inbred lines Return values is a \code{data.frame} with one row per marker pair or a \code{matrix} with all marker pairs. Additional, the euclidian distance between position of markers can be computed
}
\usage{
pairwiseLD(gpData, chr = NULL, type = c("data.frame", "matrix"),use.plink=FALSE, ld.threshold=0, rm.unmapped = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData}{
object of class \code{gpData} with elements \code{geno} and \code{map}
}                                 
  \item{chr}{
\code{numeric}. Return value is a list with pairwise LD of all markers for each chromosome in \code{chr}.
}                                      
  \item{type}{
\code{character}. Specifies the type of return value.
}
  \item{use.plink}{
\code{logical}. Should the software PLINK be used for the computation?
}
  \item{ld.threshold}{
\code{numeric}. Threshold for the LD. Only pairwise LD>\code{ld.threshold} is reported when PLINK is used. This argument is can only be used for \code{"type=data.frame"}.  
}
  \item{rm.unmapped}{
\code{logical}. Remove markers with unknown postion in \code{map} before using PLINK?
}
}
\details{
The function \code{write.plink} to prepare the input files and the script for PLINK. The executive PLINK file must be available (e.g. in the working directory). The function \code{pairwiseLD} starts PLINK and reads the results. The evaluation is performed separate for every chromosome. The measure for LD is \eqn{r^2}{r2}. This is defined as 
\deqn{D= p_{AB} - p_Ap_B }{D = p(AB) - p(A)*p(B)} and
\deqn{r^2=\frac{D^2}{p_Ap_Bp_ap_b}}{r2=D2/p(A)p(B)p(a)p(b)}
where \eqn{p_{AB}}{p(AB)} is defined as the observed probability of haplotype \eqn{AB}, \eqn{p_A=1-p_a} and \eqn{p_B=1-p_b} the observed probabilities of alleles \eqn{A} and \eqn{B}. 
If the number of markers is higt, a threshold for the LD can be used. In this case, only pairwise LD above the threshold is reported (argument \code{--ld-window-r2 in PLINK}).
}
\value{
For \code{type="data.frame"} a list with one element for each chromosome. Each element is a \code{data.frame} with columns \code{marker1}, \code{marker2}, \code{r2} and \code{distance} for all \eqn{p(p-1)/2} marker pairs.

For \code{type="matrix"} a list with one element for each chromosome. Each element is a list of 2: a \eqn{p \times p}{p x p} \code{matrix} with pairwise LD and a  \eqn{p \times p}{p x p} \code{matrix} with pairwise distances.

}
\references{
Hill WG, Robertson A (1968). Linkage Disequilibrium in Finite Populations. Theoretical
and Applied Genetics, 6(38), 226 - 231.

Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007)
PLINK: a toolset for whole-genome association and population-based linkage analysis. American Journal of Human Genetics, 81.
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
maizeLD <- pairwiseLD(maizeC,chr=1,type="data.frame")
}
