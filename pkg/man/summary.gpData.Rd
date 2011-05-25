\name{summary.gpData}
\alias{summary.gpData}
\alias{print.summary.gpData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary for class gpData
}
\description{
S3 \code{summary} method for objects of class \code{gpData}
}
\usage{
\method{summary}{gpData}(object,...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
object of class \code{gpMod}
}
\item{...}{
not used
}
}

\author{
Valentin Wimmer
}

\examples{
# Example from Henderson (1977)
dat <- data.frame(y=c(132,147,156,172),time=c(1,2,1,2),animal=c(1,2,3,4))
ped <- create.pedigree(ID=c(6,5,1,2,3,4),Par1=c(0,0,5,5,1,6),Par2=c(0,0,0,0,6,2))
gp <- create.gpData(pheno=dat,pedigree=ped)
A <- kin(gp,ret="add")

# h2 =0.5, vare = 1 -> varg = 1
mod1 <- list(fit=list(sigma=c(1,1)),kin=A,model="BLUP",trainingSet=paste(1:4),y=c(132,147,156,172))
class(mod1) <- "gpMod"
predict(mod1,c("5","6"))

# 'by hand'
y <- c(132,147,156,172)
X <- matrix(1,ncol=1,nrow=4)
Z <- diag(6)[-c(1,2),]
AI <- solve(A)
RI <- diag(4)
 
res <- MME(X,Z,AI,RI,y)
res$b + res$u[1:2]

}
\keyword{methods}

