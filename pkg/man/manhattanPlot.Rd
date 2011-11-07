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
manhattanPlot(b, gpData = NULL, colored = FALSE, add = FALSE, pch = 19, ylab = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{
numeric vector of effects to plot
}
  \item{gpData}{
object of class \code{gpData} with map position
}
  \item{colored}{
  Is the plot colored or in grey. The \code{default} is "\code{grey}"
}
  \item{add}{
  If \code{TRUE} the plot is added to an existing plot. The default is \code{FALSE}.
}
  \item{pch}{
  a vector of plotting characters or symbols: see \code{\link[graphics]{points}}. The default is an open circle.
}
  \item{ylab}{
  a title for the y axis: see \code{\link[graphics]{title}}.
}
  \item{\dots}{
further arguments for function \code{plot}
}
}

\author{
Hans-Jürgen Auinger
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
