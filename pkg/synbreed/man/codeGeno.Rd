\name{codeGeno}
\alias{codeGeno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{                                            
Recoding of genotypic data, imputing of missing values and preselection of markers 
}
\description{
Function to code genotypic data in an object of class \code{gpData}. Raw genotypic data could be in any format (e.g. alleles coded by alleles "A/T","G/C", ... , or by genotypes "AA", "BB", "AB"). Function is suitable for biallelic markers with a maximum of 3 genotypes per locus. Raw data is recoded as the number of counts of the minor allele, i.e. 0, 1 and 2. 
Imputing of missing values can be done with a fixed value, random sampling from allele distribution or family information (see details). Additional preselection of markers can be done according to the minor allele frequency and/or fraction of missing values.                                     
}
\usage{
codeGeno(gpData,impute=FALSE,impute.type=c("fix","random","family"),
replace.value=NULL,maf=NULL,nmiss=NULL,label.heter="AB",
keep.identical=TRUE,verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData}{
object of class \code{gpData} with marker matrix with arbitrary coding. Missing values have to be coded as \code{NA}.
}           
  \item{impute}{
logical scalar. Should missing value be replaced by imputing?
}
  \item{impute.type}{
character with one out of "fix", "random", or "family". See details.
}
  \item{replace.value}{
numeric value to replace missing values in case \code{impute.type="fix"}. 
}
  \item{maf}{
numeric scalar. Markers with minor allele frequency < \code{maf} are discarded, thus  \code{maf} in [0,0.5]. If \code{map} in \code{gpData} is available, markers are also removed from \code{map}.
}
  \item{nmiss}{
numeric scalar. Markers with more than \code{nmiss} fraction of missing values are discarded, thus  \code{nmiss} in [0,1]. If \code{map} in \code{gpData} is available, markers are also removed from \code{map}.
}
  \item{label.heter}{
This is either a scalar or vector of characters to identify heterozygous genotypes or a function which evaluates if an element of the marker matrix is a heterozygous genotype. Defining a function is useful, if number of unique heterozygous genotypes is high, i.e. if genotypes are coded by alleles. Note that heterozygus values must be identified unambiguously by \code{label.heter}, e.g. "/" is not valid if genotypes are coded "A/A", "A/B" and "B/B". Use \code{label.heter=NULL} if there are only homozygous genotypes, i.e. in DH lines, to speed up computation.
}
  \item{keep.identical}{
logical scalar. Should duplicated markers be kept, i.e. markers which have a LD of 1? 
}
 \item{verbose}{
	logical scalar. If \code{TRUE} verbose output is generated during the steps of the algorithm. This is usefull if selection of markers is performed as in this case information on discarded markers is given.
 }
}                                         
\details{
Coding of genotypic data is done in the following order (depending on choice of arguments not all steps are performed):
1. Discarding markers due with fraction > \code{nmiss} of missing values

2. Recoding alleles from character/factor/numeric as the number of the minor alleles, i.e. 0, 1 and 2

3. Replace of missing values by \code{replace.value}  or impute missing values according to one of the following methods:

Imputing is done according to \code{impute.type}
\describe{
\item{"fix"}{All missing values are imputed by \code{replace.value}}. Note that in most cases only value 1 should be choosen.
\item{"random"}{The missing values for a marker \eqn{j} are sampeled from the marginal allele distribution of marker \eqn{j}. With 2 possible genotypes, i.e. 0 and 2, values are sampled from distribution with probabilites \eqn{P(x=0)=1-p} and \eqn{P(x=2)=p}, where \eqn{p} is the minor allele frequency of marker \eqn{j}. To use this distribution for the sampling of missing values, specify \code{label.heter=NULL}. In case of 3 genotypes, i.e. with heterozygous genotypes, values are sampled from distribution \eqn{P(x=0)=(1-p)^2}, \eqn{P(x=1)=p(1-p)} and \eqn{P(x=2)=p^2}}. 
\item{"family"}{
Suppose an observation \eqn{i} is missing (NA) for a marker \eqn{j} in population \eqn{k}. If marker \eqn{j} is fixed in population \eqn{k}, the imputed value will be the fixed allele. If marker \eqn{j} is segregating for the population \eqn{k},                                          
the value is 0 with probability 0.5 and 1 with probability 0.5. To impute with family information, a column named "family" in element \code{covar} of \code{gpData} is necessary. This column should contain a \code{charater} or \code{numeric} to identify family of all genotyped individuals. Family information for nongenotyped individuals is not used an could be set to any value, i.e. NA}                          
}

4. Recoding of alleles after imputation, if necessary due to changes after imputation

5. Discarding markers with a minor allele frequency of <= \code{maf}

6. Discarding duplicated markers if \code{keep.identical=FALSE}

7. Restoring original data format (\code{gpData}, \code{matrix} or \code{data.frame}) 
}
\value{
An object of class \code{gpData} containing the recoded marker matrix. If \code{maf} or \code{nmiss} were specified or \code{keep.identical=FALSE}, dimension od \code{geno} and \code{map} may be reduced due to selection of markers  The genotype which is homozygous for the minor allele is coded as 2, the other homozygous is coded as 0 and heterozygous genotype is coded as 1.
}

\author{
Valentin Wimmer
}
\note{
Information about imputing is printed on the screen (approximative run time, number of imputations by family structure and number of random imputations).
}

\examples{
# create marker data for 9 SNPs and 10 homozygous individuals 
snp9 <- matrix(c(
  "AA",   "AA",   "AA",   "BB",   "AA",   "AA",   "AA",   "AA",  NA,
  "AA",   "AA",   "BB",   "BB",   "AA",   "AA",   "BB",   "AA",  NA,
  "AA",   "AA",   "BB",   "BB",   "AA",   "AA",   "AA",   "BB",  NA,
  "AA",   "AA",   "BB",   "BB",   "AA",   "AA",   "AA",   "AA",  NA,
  "AA",   "AA",   "BB",   "BB",   "AA",   "BB",   "BB",   "BB",  NA,
  "AA",   "AA",   "BB",   "BB",   "AA",   NA,     "BB",   "AA",  NA,
  "BB",   "AA",   "BB",   "BB",   "BB",   "AA",   "BB",   "BB",  NA,
  "AA",   "AA",    NA,    "BB",    NA,    "AA",   "AA",   "AA",  "AA",
  "AA",    NA,     NA,    "BB",   "BB",   "BB",   "BB",   "BB",  "AA",
  "AA",    NA,    "AA",   "BB",   "BB",   "BB",   "AA",   "AA",  NA),
  ncol=9,byrow=TRUE)

# set names for markers and individuals
colnames(snp9) <- paste("SNP",1:9,sep="")  
rownames(snp9) <- paste("ID",1:10+100,sep="")

# create object of class 'gpData'
# two families A and B
fam <- data.frame(family=c(rep("A",7),rep("B",3)))
rownames(fam) <- paste("ID",1:10+100,sep="")

gp <- create.gpData(geno=snp9,family=fam)  

# code genotypic data  
gp.coded <- codeGeno(gp)

# impute missing values by family information
gp.imputed <- codeGeno(gp,impute=TRUE,impute.type="family")


# example with heterogenous stock mice
data(mice)
summary(mice)
# heterozygous values must be labeled  (takes some minutes)
# mice.coded <- codeGeno(mice,label.heter=function(x) substr(x,1,1)!=substr(x,3,3)) 

                               
# example with maize data and imputing by family
data(maize)
# first only recode alleles
maize.coded <- codeGeno(maize)

# set 200 random choosen values to NA
set.seed(123)
ind1 <- sample(1:nrow(maize.coded $geno),200)
ind2 <- sample(1:ncol(maize.coded $geno),200)
original <- maize.coded$geno[cbind(ind1,ind2)]

maize.coded$geno[cbind(ind1,ind2)] <- NA
# imputing of missing values by family structure
maize.imputed <- codeGeno( maize.coded,impute=TRUE,impute.type="family")


# compare in a crosstable
imputed <- maize.imputed$geno[cbind(ind1,ind2)]
(t1 <- table(original,imputed) )
# sum of correct replacements
sum(diag(t1))/sum(t1)

# compare with random imputation
maize.random <- codeGeno( maize.coded,impute=TRUE,impute.type="random")
imputed2 <- maize.random$geno[cbind(ind1,ind2)]
(t2 <- table(original,imputed2) )
# sum of correct replacements
sum(diag(t2))/sum(t2)
}
\keyword{manip}
