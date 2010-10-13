\name{summary.pedigree}
\alias{summary.pedigree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing pedigree information
}
\description{
\code{summary} method for class ''pedigree''
}
\usage{
\method{summary}{pedigree}(object,...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
object of class ''pedigree''
}
\item{...}{
not used
}
}

\author{
Valentin Wimmer
}

\examples{
# plant pedigree
ped <- simul.pedigree(gener=4,7)
summary(ped)

# animal pedigree
ped <- simul.pedigree(gener=4,7,animals=TRUE)
summary(ped)
}
\keyword{methods}

