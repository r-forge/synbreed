LDDist <- function(marker,linkageGroup,pos,file=NULL,fit=FALSE,...){

    if (length(linkageGroup)!=ncol(marker)) stop("'linkageGroup' must be of length ncol(marker)")
    if (length(pos)!=ncol(marker)) stop("'pos' must be of length ncol(marker)")
    
    lg <- unique(linkageGroup)
    
    # function for fit according to ...
    smooth.fit <- function(overallDist,overallr2,n){
    
      nonlinearoverall <- nls(overallr2 ~ ((10 + p*overallDist)) / ((2+p*overallDist) * (11 + p*overallDist) ) *
      ( 1 + ( (3+ p*overallDist) * (12 + 12 * p + p^2*overallDist^2)) / ( n*(2+p*overallDist) * (11 + p*overallDist))),
      start=list(p=1))

      p <- coef(nonlinearoverall)

      fitcurve <- function(x,p,n) {
        ((10 + p*x)) / ((2+p*x) * (11 + p*x) ) *
        ( 1 + ( (3+ p*x) * (12 + 12 * p + p^2*x^2)) / ( n*(2+p*x) * (11 + p*x)))
      }

      curve(fitcurve(x,p=p,n=n), from=min(overallDist), to = max(overallDist), add=TRUE,col=2,lwd=2)
    }
    

    # initialize return data list
    ret <- list()

    if(!is.null(file)) pdf(file)
    # compute distances within each linkage group
    for (i in 1:length(lg)){

       markeri <- marker[,linkageGroup==lg[i]]
       p <- ncol(markeri)
       mn <- colnames(markeri)
       posi <- pos[linkageGroup==lg[i]]
       ld.r2i <- cor(markeri,method="spearman")^2
       ld.r2i <- ld.r2i[lower.tri(ld.r2i)] # column-wise
       rowi <- rep(1:p,times=(p:1)-1)
       coli <- p+1 - sequence(1:(p-1))
       coli <- coli[length(coli):1]
       disti <- abs(posi[rowi] - pos[coli])

       ret[[i]] <- data.frame(marker1=mn[rowi],marker2=mn[coli],r2=ld.r2i,dist=disti)

       plot(r2~dist,data=ret[[i]],main=paste("Linkage Group",i),...)
       if(fit) smooth.fit(ret[[i]][,4],ret[[i]][,3],n=nrow(markeri)) 

    }
     if(!is.null(file)) dev.off()
    
    invisible(ret)

}



