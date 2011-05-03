\name{summary.gpMod}
\alias{summary.gpMod}
\alias{print.summary.gpMod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary for class gpMod
}
\description{
S3 \code{summary} method for objects of  class \code{gpMod} }
\usage{
summary.gpMod(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
object of class \code{gpMod}
}
  \item{\dots}{
not used
}
}
\seealso{
\code{\link{gpMod}}
}
\examples{
data(maize)
maizeC <- codeGeno(maize)
# marker-based (realized) relationship matrix
U <- kin(maizeC,ret="realized")/4

# BLUP model
mod <- gpMod(maizeC,model="BLUP",kin=U)
summary(mod)
}

\keyword{methods}
