\name{gpMod}
\alias{gpMod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genomic predictions models for objects of class gpData
}
\description{
This function can be used to fit genomic prediction models to an object of class \code{gpData}. The possible models are Best Linear Unbiased Prediction (BLUP) using a pedigree-based or marker-based genetic relationship matrix and Bayesian Lasso (BL) or Bayesian Ridge regression (BRR).
BLUP models are fitted using the implementation of the \code{regress} package (Clifford and McCullagh, 2009). The Bayesian regression models are fitted using the implementaion of the \code{BLR} package (de los Campos and Perez, 2010). The covariance structure in the BLUP model is defined by an opject of class \code{relationshipMatrix}. The training set for the model fit consists of all individuals with phenotypes and genotypes. All data is restricted to individuals from the training set.


}
\usage{
gpMod(gpData, model = c("BLUP", "BL", "BRR"), kin = NULL, trait = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData}{
object of class \code{gpData}
}
  \item{model}{
\code{character}. The model to use.
}
  \item{kin}{
object of class \code{relationshipMatrix}. Kinship structure to be used in the BLUP model or as addtional polygenic effect \eqn{u} in the Bayesian regression models.
}
  \item{trait}{
\code{numeric} or \code{character}. The name or number of the trait to fit the model
}
  \item{\dots}{
further arguments to be used by the genomic prediction models, i.e. prior values and MCMC options for the \code{BLR} function (see \code{?BLR}).
}
}
\details{
Only a subset of the individuals - the training set - is used for the model fit. This contains all individuals with phenotypes and genotypes. If \code{kin} does not match the dimension of the training set (if, e.g. ancestors are included), the respective rows and columns from the trainings set are choosen.
}
\value{
Object of class \code{gpMod} which is a list of
\item{fit}{The model fit returned by the genomic prediction method}
\item{model}{The model type, see 'Arguments'}
\item{y}{The phenotypic records for the individuals in the training set}
\item{g}{The predicted genetic values for the individuals in the training set}
\item{m}{Predicted SNP effects (for BL and BRR)}
\item{kin}{Matrix \code{kin}}
}
\not{
The verbose output of the \code{BLR} function is written to a file \code{BLRout.txt} to prevent the screen output from overload
}
\references{
Clifford D, McCullagh P (2009). regress: Gaussian Linear Models with Linear Covariance
Structure. R package version 1.1-2, URL http://www.csiro.au.

Gustavo de los Campos and Paulino Perez Rodriguez, (2010). BLR:
Bayesian Linear Regression. R package version 1.2.
http://CRAN.R-project.org/package=BLR

}
\author{
Valentin Wimmer
}
\seealso{
\code{\link{kin}}, \code{\link{crossVal}}
}
\examples{
data(maize)
maizeC <- codeGeno(maize)

# pedigree-based (expected) kinship matrix
A <- kin(maizeC,ret="kin",DH=maize$covar$DH)

# marker-based (realized) relationship matrix
U <- kin(maizeC,ret="realized")/4

# BLUP models
mod1 <- gpMod(maizeC,model="BLUP",kin=A)
mod2 <- gpMod(maizeC,model="BLUP",kin=U)

# Bayesian Lasso
prior <- list(varE=list(df=3,S=35),lambda = list(shape=0.52,rate=1e-4,value=20,type='random'))
mod3 <- gpMod(maizeC,model="BL",prior=prior,nIter=6000,burnIn=1000,thin=5)

summary(mod1)
summary(mod2)
summary(mod3)

}
