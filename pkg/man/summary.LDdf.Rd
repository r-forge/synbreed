\name{summary.LDdf}
\Rdversion{1.1}
\alias{summary.LDdf}
\alias{summary.LDmat}

\title{
Summary for LD objects
}
\description{
\code{Summary} method for class ''LDdf''  and  ''mat''
}
\usage{
\method{summary}{LDdf}(LDdf) 
\method{summary}{LDmat}(LDmat) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LDdf}{
object of class \code{LDdf} which is the output of function \code{pairwiseLD} and argument \code{type="data.frame"}
}
  \item{LDmat}{
object of class \code{LDmat} which is the output of function \code{pairwiseLD} and argument \code{type="matrix"}
}
}
\details{
Returns for each chromosome: Number of markers; mean, minimum and maximum LD measured as \eqn{r^2}{r2};  fraction of markers with \eqn{r^2>0.2}{r2>0.2};  maximum distance of markers
}
\author{
Valentin Wimmer
}

\seealso{
\code{\link{pairwiseLD}}, ~~~
}
\examples{
  \dontrun{
  data(maize)
  maizeC <- codeGeno(maize)
  maizeLD <- pairwiseLD(maizeC,chr=1:10,type="data.frame")
  maizeLDm <- pairwiseLD(maizeC,chr=1:10,type="matrix")
  summary(maizeLD)
  summary(maizeLDm)
  }
}

\keyword{methods}

