# computing genomic relationship matrix according to vanRaden

vanRaden <- function(marker){
    if(!is.matrix(marker)) marker <- as.matrix(marker)
    # M coded with 0 / 2
    M <- marker
    n <- nrow(M)
    p <- ncol(M)
    maf <- colMeans(M,na.rm=TRUE)
    P <- matrix(rep(maf,each=n),ncol=p)
    Z <- M - P
    Zq <- tcrossprod(Z)
    G <- Zq/(2*sum(maf/2*(1-maf/2)))
    class(G) <- "relationshipMatrix"
    return(G)
}
