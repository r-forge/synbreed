\name{plot.pedigree}
\alias{plot.pedigree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualization of pedigree
}
\description{
A function to visualize pedigree structure by a graph. Each genotype is represented as vertex and direct offsprings are linked by an edge. 
}
\usage{
\method{plot}{pedigree}(x, effect = NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
object of class \code{pedigree} or object of class \code{gpData} with element \code{pedigree}
}
  \item{effect}{
vector of length \code{nrow(pedigree)} with effect to for the x axis
}
  \item{...}{
Other arguments for function \code{igraph.plotting}
}
}                                                                    
\details{
 The pedigree is structured top-bottom. In the first line the first generation is printed. Links over more than one generation are possible as well as genotypes with only one (known) parent. Usually, no structure in one generation is plotted. If an \code{effect} is given, the genotypes are ordered by this effect and an labeled axis is drawn.
}                                
\value{
A named graph visualizing the pedigree structure
}

\author{
Valentin Wimmer
}
\note{
This function uses the plotting method for graphs in the libraray \code{igraph}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{create.pedigree}}, \code{\link{simul.pedigree}}
}
\examples{
# example with 9 individuals        
id <- 1:9
par1 <- c(0,0,0,0,1,1,1,4,7)
par2 <- c(0,0,0,0,2,3,2,5,8)
gener <- c(0,0,0,0,1,1,1,2,3)

# create pedigree object 
ped <- create.pedigree(id,par1,par2,gener)
plot(ped)               
}
\keyword{hplot}

