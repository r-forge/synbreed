kinshipPlants <- function (ped,inbreed=NULL) {

    # read data from pedigree
    id <- ped$ID
    par1.id <- ped$Par1
    par2.id <- ped$Par2

    if (is.null(inbreed)) inbreed <- (par1.id == par2.id) &  (par1.id != 0) &   (par2.id != 0)
    else inbreed <- as.logical(inbreed) #TRUE or 1 if genotype ist fully homozygous, otherwise FALSE or 0
    library(kinship)
    n <- length(id)
    if (any(duplicated(id)))
        stop("All id values must be unique")
    kmat <- diag(n + 1)
    kmat[n + 1, n + 1] <- 0
    pdepth <- kindepth(id, par1.id, par2.id)
    drow <- match(par1.id, id, nomatch = n + 1)
    mrow <- match(par2.id, id, nomatch = n + 1)
    for (depth in 0:max(pdepth)) {
        indx <- (1:n)[pdepth == depth]
        for (i in indx) {
            par1 <- drow[i]
            par2 <- mrow[i]
            kmat[i, ] <- kmat[, i] <- (kmat[par1, ] + kmat[par2,
                ])/2
                if (!inbreed[i])
                    kmat[i, i] <- 1 + kmat[par1, par2]/2 else
                    kmat[i, i] <- 2
        }
    }
    kmat <- kmat[1:n, 1:n]
    dimnames(kmat) <- list(id, id)
    class(kmat) <- "relationshipMatrix"
    return(kmat)
}
