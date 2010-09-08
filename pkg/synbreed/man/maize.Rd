\name{maize}
\alias{maize}
\docType{data}
\title{
Simulated maize data
}
\description{
This datatset contains genotypic data, phenotypic data and pedigree data for a simulated dataset of 1250 doubled haploid (DH) lines of maize.  The DH lines derived on 25 crosses with 50 DHs for each cross. Pedigree information is available for up to three generations.
}                                                                                                                                                                                                                 
\usage{data(maize)}
\format{
  4 datasets
  \describe{
    \item{\code{maize.geno}}{A \code{data.frame} containing genotypic data of 696 biallelic SNP markes for the 1250 genotpyes identified with \code{ID}}
    \item{\code{maize.pheno}}{A \code{data.frame} containing the phenotypic data with two columns \code{ID} and \code{Trait}.}
    \item{\code{maize.ped}}{A \code{data.frame} containing the pedigree data for the 1250 genotypes and 51 individuals for pedigree \describe{
    \item{\code{ID}}{Identification number for the genotypes. Note that the first 51 rows are not genotyped and phenotype and the last 1250 rows are the genotyped and phenotype individuals.}
    \item{\code{Par1}}{Parent 1 in pedigree. Zero if unknown parent 1}\item{\code{Par2}}{Parent 2 in pedigree. Zero if unknown parent 2} \item{\code{DH}}{logical. \code{TRUE} if the genotype is a DH line}}}
    \item{\code{maize.marker.pos}}{A \code{data.frame} with two columns for the physical positions of the 696 marker \describe{ \item{\code{chr}}{A numeric giving the chromosome of the marker} \item{\code{pos}}{A numeric giving the position on the chromosome in cM}  }}
  }
}
\details{
}
\source{
}
\references{
}                                  
\examples{
data(maize)
# create.pedigree
ped <- create.pedigree(maize.ped$ID,maize.ped$Par1,maize.ped$Par2)
#  plot small part of the pedigree
plot(ped[c(7,10,27,55:58,8,29,190:195),])
# compute relationship matrix according to kinship
Akin <- kinshipPlants(ped,inbreed=maize.ped$DH)
plot(Akin)

# code genotypic data
# do not use first column (ID column)
marker <- codeGeno(maize.geno[,-1])
# histogramm of minor allele frequency
hist(colMeans(marker)/2,xlab="minor allele frequency",main="")
                                              

# compute the genomic relationship matrix
# 1) vanRaden
vR <- vanRaden(marker)                   
plot(vR)     
# 2) rogers
rog <- rogers(marker,correction="Hayes")           
plot(rog)
                                         
}                                  
\keyword{datasets}
