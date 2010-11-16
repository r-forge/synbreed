\name{plot.pedigree}
\alias{plot.pedigree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualisaton of pedigree
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
object of class \code{pedigree} or object of class \code{gpData} with filled slot \code{pedigree}
}
  \item{effect}{
vector of length \code{nrow(pedigree)} with effect to be drawn
}
  \item{...}{
Other arguments for function \code{igraph.plotting}
}
}                                                                    
\details{
 The pedigree is structured top-bottom. In the first line the first generation is printed. Links over more than one generation are possible as well as genotypes with only one (known) parent. Usually, no structure in one generation is plotted. If an \code{effect} is given, the genotypes are ordered by this effect and an labeled axis is drawn.
}                                
\value{
A named graph visualising the pedigree structure
}

\author{
Valentin Wimmer
}
\note{
This function uses the plotting method for graphs in the libaray \code{igraph}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{create.pedigree}}, \code{\link{simul.pedigree}}
}
\examples{
ped.raw <- matrix(c(1,  0,   0,   1, 
                 2,  0,   0,   1,
                 3,  0,   0,   1,
                 4,  0,   0,   1,
                5,  1,   2,   2, 
                 6,  1,   3,   2,
                 7,  1,   2,   2,
                 8,  1,   2,   2) ,
                 byrow = TRUE, ncol = 4)
ped <- create.pedigree(ped.raw[,1],ped.raw[,2],ped.raw[,3],ped.raw[,4])
plot(ped)                 
}
\keyword{hplot}

