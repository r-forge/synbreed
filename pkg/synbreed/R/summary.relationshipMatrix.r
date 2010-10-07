summary.relationshipMatrix <- function(relationshipMatrix,...){
      relMat <- relationshipMatrix
     cat("Dimension         :",nrow(relMat),"x",ncol(relMat),"\n")
     cat("Rank              :",qr(relMat)$rank,"\n")
     cat("Range             :",min(relMat),"--",max(relMat),"\n")
     cat("# of unique values:",length(unique(relMat[upper.tri(relMat,diag=TRUE)])),"\n")
}