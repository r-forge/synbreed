\name{LDMap}
\Rdversion{1.1}
\alias{LDMap}

\title{
LD Heatmap
}
\description{
Computation of pairwise LD as \eqn{R^2} for marker data and plot of LD heatmap for each linkage group (chromosome)
}
\usage{
LDMap(LDmat,gpData,file=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LDmat}{
Object of class \code{gpData} where \code{geno} and \code{map} are available.
}
  \item{gpData}{
Object of class \code{gpData} that was used to infer \code{LDmat}.
}
  \item{file}{
Optionally a path to a file were plot is saved as \code{pdf}.
}
  \item{\dots}{
Further arguments that could be passed to function \code{LDheatmap}.
}
}

\seealso{
\code{\link{pairwiseLD}}, \code{\link[LDheatmap]{LDheatmap}}, \code{\link{LDDist}}
}
\examples{
data(maize)
maizeC <- codeGeno(maize)
\dontrun{
# LD for chr 1
maizeLD <- pairwiseLD(maizeC,chr=1,type="matrix")
LDMap(maizeLD,maizeC)
}
}
\keyword{hplot}
