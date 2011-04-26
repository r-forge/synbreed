\name{manhattanPlot}
\alias{manhattanPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Manhattan plot for SNP effects
}
\description{
Plot of SNP effect along the chromosome
}
\usage{
manhattanPlot(b, gpData = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{
numeric vector of effects to plot
}
  \item{gpData}{
object of class \code{gpData} with map position
}
  \item{\dots}{
further arguments for function \code{plot}
}
}

\author{
Valentin Wimmer
}

\examples{
data(mice)
b <- rexp(ncol(mice$geno),3)
manhattanPlot(b,mice)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
