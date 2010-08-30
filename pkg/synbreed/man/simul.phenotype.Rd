\name{simul.phenotype}                                                                                                                   
\alias{simul.phenotype}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of field trial
}
\description{                                                                                                                             
This function could be used to simulate for a given pedigree results of a quantitative trait in a field trial with several replications of each genotype and several locations.
}
\usage{
simul.phenotype(pedigree, h2, loc = 1, repl = 1, seed = 123, genotypeEff = NULL, locEff = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pedigree}{
An object of class \code{pedigree} given the pedigree structure
}
  \item{h2}{    
Numeric: heretability of the trait  
}
  \item{loc}{
Numeric: number of locations 
}
  \item{repl}{
Numeric: number of replications in each location 
}
  \item{seed}{
Numeric: starting value for use in \link{set.seed}
}
  \item{genotypeEff}{
Numeric vector of length N containing the base level in simulation for each genotype in the order given by \code{pedigree}. If \code{NULL}, the same value 100 is used for every genotype in pedigree.
}
  \item{locEff}{
Numeric vector of length \code{loc} containg the effect in simulation for each location. If \code{NULL}, the effects will be drawn from a normal distribution.
}
}
\details{
Simulation will be done for the whole pedigree, not only the last generation.
}
\value{
A \code{data.frame} with one column for the \code{ID} of pedigree, one with simulated values of the trait and one column with the number of the location
}
\author{
Valentin Wimmer
}
\note{
in developement 
}


\seealso{
\code{\link{simul.pedigree}}
}
\examples{
# simulation of pedigree
p <- simul.pedigree()

# same effect for each genotype
dat1 <- simul.phenotype(p,h2=0.5,loc=2,repl=2)
head(dat1)

# different effect for each genotype
dat2 <- simul.phenotype(p,h2=0.5,loc=2,repl=2,genotypeEff=1:nrow(p))
head(dat2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{simulation}
