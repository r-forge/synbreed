summary.relationshipMatrix <- function(object,...){
     relMat <- object
     ans <- list(dim=c(nrow=nrow(relMat),ncol=ncol(relMat)),rank=qr(relMat)$rank,range.off.diagonal=c(min=min(relMat[upper.tri(relMat,diag=FALSE)]),max=max(relMat[upper.tri(relMat,diag=FALSE)])),nUnique=length(unique(relMat[upper.tri(relMat,diag=TRUE)])),inbr.coef=summary(diag(relMat)-1))
     class(ans) <- "summary.relationshipMatrix"
     ans
}

print.summary.relationshipMatrix <- function(x,...){
    cat("\t dimension                   ",x$dim[1],"x",x$dim[2],"\n")
    cat("\t rank                        ",x$rank,"\n")
    cat("\t range of off-diagonal values",x$range.off.diagonal[1],"--",x$range.off.diagonal[2],"\n")
    cat("\t number of unique values     ",x$nUnique,"\n")
    cat("\t range of inbreeding coef.   ",x$inbr.coef[1],"--",x$inbr.coef[6],"\n")
}