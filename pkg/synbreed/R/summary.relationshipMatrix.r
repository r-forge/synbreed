summary.relationshipMatrix <- function(relationshipMatrix){
      relMat <- relationshipMatrix
     cat("Dimension         :",nrow(relMat),"x",ncol(relMat),"\n")
     cat("Rank              :",sum(diag(relMat%*%ginv(relMat))),"\n")
     cat("Range             :",min(relMat),"--",max(relMat),"\n")
     cat("# of unique values:",length(unique(relMat[upper.tri(relMat,diag=TRUE)])),"\n")
}