\name{gpData2cross}
\alias{gpData2cross}
\alias{cross2gpData}
\title{
Conversion between objects of class 'cross' and 'gpData'
}
\description{
Conversion of object class \code{gpData} to object class \code{cross} (F2 intercross) of package \code{qtl} and vice versa. Function \code{codeGeno} is applied in \code{cross2gpData}, if not done before. 
}
\usage{
gpData2cross(gpData,...)
cross2gpData(cross)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData}{
object of class \code{gpData} with non-empty elements for \code{pheno}, \code{geno} and \code{map}. 
}
  \item{cross}{
object of class \code{cross}.
}
  \item{...}{
further arguments for function \code{codeGeno}. Only used in \code{gpData2cross}.
}
}
\details{
In \code{cross}, genotypic data is splitted into chromosomes while in \code{gpData} genotypic data comprises all chromosomes. 
Note that coding of genotypic data differs between classes. In \code{gpData}, genotypic data is coded as the number of copies of the minor allele, i.e. 0,1 and 2. Thus, function \code{codeGeno} should be applied to \code{gpData} before using \code{gpData2cross} to ensure correct coding. In \code{cross}, coding for F2 intercross is: AA = 1, AB = 2, BB = 3. When using  \code{gpData2cross} or \code{cross2gpData}, resulting genotypic data has correct format.
}
\value{
Object of class \code{cross} or \code{gpData} for function  \code{gpData2cross} or \code{cross2gpData}, respectively.
}
\references{
Broman, K. W. and Churchill, S. S. (2003). R/qtl: Qtl mapping in experimental crosses. Bioinformatics, (19):889-890.
}
\author{
Valentin Wimmer and Hans-Jürgen Auinger
}

\seealso{
\code{\link{create.gpData}},  \code{\link[qtl]{read.cross}} ,  \code{\link{codeGeno}}
}
\examples{
# from gpData to cross
data(maize)
maizeC <- codeGeno(maize)
maize.cross <- gpData2cross(maizeC)
# descriptive statistics
summary(maize.cross)
plot(maize.cross)

# search for QTL on chr 1
maize.cross <- calc.genoprob(maize.cross, step=2.5)
result <- scanone(maize.cross, pheno.col=1,method="em")
# display of LOD curve
plot(result) 


# from cross to gpData
data(fake.f2)
fake.f2.gpData <- cross2gpData(fake.f2)
summary(fake.f2.gpData)




}
