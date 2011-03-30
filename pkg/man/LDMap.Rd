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
LDMap(gpData,chr=NULL,file=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData}{
Object of class \code{gpData} where \code{geno} and \code{map} are available.
}
  \item{chr}{
Character or numeric. If specified, only one plot for markers in linkage group \code{chr} is created. If \code{chr="all"}, all pairwise LD between markers across chromosomes is used.  
}
  \item{file}{
Optionally a path to a file were plot is saved as \code{pdf}.
}
  \item{\dots}{
Further arguments that could be passed to function \code{LDheatmap}.
}
}
\details{
LD is computed as coefficient of determination \eqn{R^2}. The plot is simply a call of function \code{LDheatmap} in package  \code{LDheatmap}.  
}
\value{
A list with elements 
\item{LD}{list with a \code{matrix} containing all pairwise values of \eqn{R^2} or markers for each element of \code{chr}}
\item{distance}{list with a \code{matrix} containing all pairwise euclidian distances of markers for each element of \code{chr}}
}
\seealso{
\code{\link{codeGeno}}, \code{\link[LDheatmap]{LDheatmap}}, \code{\link{LDDist}}
}
\examples{
data(maize)
maize <- codeGeno(maize)
LDMap(maize,chr=1)
}
\keyword{hplot}
