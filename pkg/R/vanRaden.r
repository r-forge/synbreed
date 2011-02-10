# computing genomic relationship matrix according to vanRaden (2008)

vanRaden <- function(marker){
    
    # extract information from arguments
          if(any(class(marker)=="gpData")){
             if(!marker$info$codeGeno) stop("use function 'codeGeno' before using 'vanRaden'") 
             marker <- marker$geno
          }
          else marker <- marker

    if(!is.matrix(marker)) marker <- as.matrix(marker)
    
    # M supposed to be coded with 0,1,2
    M <- marker
    n <- nrow(M)
    p <- ncol(M)
    
    # 2* minor allele frequency as expectation
    maf <- colMeans(M,na.rm=TRUE)
    P <- matrix(rep(maf,each=n),ncol=p)
    
    # compute realized relationship matrix G
    Z <- M - P
    Zq <- tcrossprod(Z)
    G <- Zq/(2*sum(maf/2*(1-maf/2)))
    
    # return
    class(G) <- "relationshipMatrix"
    return(G)
}
