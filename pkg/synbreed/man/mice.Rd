\name{mice}
\alias{mice}
\docType{data}
\title{
Heterogenous stock mice population
}
\description{
Data set comprises public available data of 2527 (993 males and 947 females) heterogenous stock mice derived from eight inbred strains (A/J, AKR/J, BALBc/J, CBA/J, C3H/HeJ, C57BL/6J, DBA/2J and LP/J) followed by 50 generations of pseudorandom mating. Pedigree is available on parents of phenotyped individuals. All individuals are labeled with a unique ID, starting with  \code{A048005080}. For all individuals, family, sex (females=0, males=1), month of birth (1-12), birthyear, coat color, cage density and litter is available and stored in \code{covar}. 

 The measured phenotypes are descibed in Solberg et al. (2006). Here, the body weight at age of 6 weeks and growth slope are available. The heritabilities of these traits are 0.74 and 0.30, respectively (Valdar et al, 2006b). Data was taken from \url{http://mus.well.ox.ac.uk/GSCAN/HS_PHENOTYPES/Weight.txt}.
 
 Genotypic data consists of 12545 biallelic SNP markers and is available for 1940 individuals. Raw genotypic data from \url{http://mus.well.ox.ac.uk/GSCAN/HS_GENOTYPES/All} is given in the \code{Ped-File Format} with two columns for each marker. Both alleles were combined to a single genotype for each marker in \code{mice} data. The SNPs are mapped in a sex-averaged genetic map with distances given in centimorgan (Shifman et al. (2006)). SNPs are mapped across all 19 autosomes and X-chromosome where distances between adjacent markers variied form 0 to 3 cM.  
}                                                                                                                                                                                                                 
\usage{data(mice)}
\format{
Object of class 'gpData'
}
\source{Wellcome Trust Centre for Human Genetics, Oxford University, data available from \url{http://gscan.well.ox.ac.uk}}   
\references{
Shifman S, Bell JT, Copley RR, Taylor MS, Williams RW, et al. (2006) A High-Resolution Single Nucleotide Polymorphism Genetic Map of the Mouse Genome. PLoS Biol 4(12): e395. doi:10.1371/journal.pbio.0040395

Solberg L.C. et al, A protocol for high-throughput phenotyping, suitable for quantitative trait analysis in mice. Mamm. Genome 17, 129-146 (2006)

Valdar W, Solberg LC, Gauguier D, Burnett S, Klenerman P, Cookson WO, Taylor MS, Rawlins JN, Mott R, Flint J. (2006a) Genome-wide genetic association of complex traits in heterogeneous stock mice. Nat Genet. 2006 Aug;38(8):879-87.

Valdar W, Solberg LC, Gauguier D, Cookson WO, Rawlins NJ, Mott R, Flint J.(2006b) Genetic and environmental effects on complex traits in mice. Genetics. 2006 Aug 3;
}                             
\examples{
data(mice)
summary(mice)                                       
}                                  
\keyword{datasets}
