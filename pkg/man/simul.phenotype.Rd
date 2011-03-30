\name{simul.phenotype}
\Rdversion{1.1}
\alias{simul.phenotype}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of a field trial with single trait
}
\description{
Simulation of a records of a field trial for given pedigree or relationship structure.
The field trial consists of multiple locations and randomized complete block design within locations. 
A single quantitative trait is simulated according to the specified contribution of relationship among individual, location and block effects to the phenotypic values.
}
\usage{
simul.phenotype(pedigree = NULL, A = NULL, mu = 100, 
vc = NULL, Nloc = 1, Nrepl = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pedigree}{
Object of class "pedigree" 
}
  \item{A}{
Object of class "relationshipMatrix"
}
  \item{mu}{
Numeric; Overall mean of the trait.
}
  \item{vc}{
List containing the variance components. \code{vc} consists of elements \code{sigma2e}, \code{sigma2a}, \code{sigma2l}, \code{sigma2b} with the variance components of the residual, the additive genetic effect,
the location effect and the block effect.
}
  \item{Nloc}{
Numeric; Number of locations in the field trial.
}
  \item{Nrepl}{
Numeric; Number of complete blocks within location.
}
}
\details{
Either \code{pedigree} or \code{A} must be specified. If \code{pedigree} is given, pedigree information is used to set up numerator relationship matrix with function \code{kinship}. If unrelated individuals should be used for simulation,
use identity matrix for \code{A}. True breeding values for $N$ individuals is  simulated according to following distribution
\deqn{tbv \sim N(0,\bf A \sigma_a^2)}{tbv = chol(A)*sigma2a*rnorm(N,0,1)}
If now location or block effects should appear, set  \code{sigma2l} or \code{sigma2b} to zero.
}
\value{
A \code{data.frame} with containing the simulated values for trait and the following variables
\item{ID}{Factor identifying the individuals. Names are extracted from \code{pedigree} or rownames of \code{A}}
\item{Loc}{Factor for Location}
\item{Block}{Factor for Block within Location}
\item{Trait}{Numeric values for trait}
\item{TBV}{Simulated values for true breeding values of individuals}
Result is sorted for locations within individuals.
}
\author{
Valentin Wimmer
}
\note{
In future multiple trait model will be implemented
}
\seealso{
\code{\link{simul.pedigree}}
}
\examples{
ped <- simul.pedigree(gener=5)
varcom <- list(sigma2e=25,sigma2a=36,sigma2l=9,sigma2b=4)
# field trial with 3 locations and 2 blocks within locations
data.simul <- simul.phenotype(ped,mu=10,vc=varcom,Nloc=3,Nrepl=2)
head(data.simul)
# analysis of variance
anova(lm(Trait~ID+Loc+Loc:Block,data=data.simul))
}

