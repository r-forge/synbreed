\name{LDDist}
\Rdversion{1.1}
\alias{LDDist}

\title{
LD versus distance Plot
}
\description{
Computation of pairwise LD as \eqn{R^2} for marker data and plot LD versus distance for each linkage group (chromosome)
}
\usage{
LDDist(marker, linkageGroup, pos, file = NULL, fit=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{marker}{
A \code{data.frame} of \code{matrix} containing the marker data. Could be the result of a call of \code{codeGeno}.
}
  \item{linkageGroup}{
A numeric vector of length \code{nrow(marker)} with the linkage group of each marker.
}
  \item{pos}{
A numeric vector of length \code{nrow(marker)} with physical position of each marker on the chromosome.
}
  \item{file}{
Optionally a path were the plot is saved as \code{pdf}.
}
  \item{fit}{
Logical. Should a nonlinear curve be fitted for the data.
}
  \item{\dots}{
Further plot arguments.
}
}
\details{
LD is computed as coefficient of determination \eqn{R^2}.  Matrix of \eqn{R^2} of markers is returned. If only fitted curves should be plotted, add argument \code{type="n"}.
}

\seealso{
\code{\link{codeGeno}}, \code{\link{LDMap}}
}
\examples{
# simulated data without LD with 1000 equal spaced markers on 4 chromosomes
M<- matrix(sample(c(0,2),size=50000,replace=TRUE,prob=c(0.65,0.35)),nrow=50,ncol=1000)
colnames(M) <- paste("M",1:1000,sep="")
lgd <- rep(1:4,each=250)
pos <- sequence(c(250,250,250,250))
par(mfrow=c(2,2))
LDDist(M,lgd,pos)
par(mfrow=c(1,1))

# maize data example
data(maize)
marker <- codeGeno(maize.geno[,-1])
LDDist(marker,maize.marker.pos$chr,maize.marker.pos$pos,fit=TRUE)
}

\keyword{LD}
\keyword{genotypic data}
