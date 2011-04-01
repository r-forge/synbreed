\name{summary.relationshipMatrix}
\alias{summary.relationshipMatrix}
\alias{print.summary.relationshipMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing relationship matrices
}
\description{
\code{summary} method for class ''relationshipMatrix''
}
\usage{
\method{summary}{relationshipMatrix}(object,...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
object of class ''relationshipMatrix''
}
  \item{...}{
not used
}
}

\author{
Valentin Wimmer
}

\examples{
data(maize)
U <- kin(codeGeno(maize),ret="realized")
summary(U)
}
\keyword{methods}



