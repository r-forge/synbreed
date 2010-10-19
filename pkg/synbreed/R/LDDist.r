LDDist <- function(marker,linkageGroup,pos,file=NULL,type="p",breaks=NULL,chr=NULL,...){

    # catch errors
    if (length(linkageGroup)!=ncol(marker)) stop("'linkageGroup' must be of length ncol(marker)")
    if (length(pos)!=ncol(marker)) stop("'pos' must be of length ncol(marker)")
    
    # set marker names if not given
    if(is.null(colnames(marker)))  colnames(marker) <- paste("M",1:ncol(marker),sep="") 

    # select chromosomes if 'chr' is specified
    lg <- unique(linkageGroup)
    if(!is.null(chr)){ 
        lg <- chr
        if(chr=="all") linkageGroup <- rep("all",length(linkageGroup))
    }
    
    # function for fit according to ...
    smooth.fit <- function(overallDist,overallr2,n){
    
      nonlinearoverall <- nls(overallr2 ~ ((10 + p*overallDist)) / ((2+p*overallDist) * (11 + p*overallDist) ) *
      ( 1 + ( (3+ p*overallDist) * (12 + 12 * p + p^2*overallDist^2)) / ( n*(2+p*overallDist) * (11 + p*overallDist))),
      start=list(p=1))

      p <- coef(nonlinearoverall)
      x <- NA
      fitcurve <- function(x,p,n) {
        ((10 + p*x)) / ((2+p*x) * (11 + p*x) ) *
        ( 1 + ( (3+ p*x) * (12 + 12 * p + p^2*x^2)) / ( n*(2+p*x) * (11 + p*x)))
      }

      curve(fitcurve(x,p=p,n=n), from=min(overallDist), to = max(overallDist), add=TRUE,...)
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
       # compute LD as R2
       ld.r2i <- cor(markeri,method="spearman")^2
       ld.r2i <- ld.r2i[lower.tri(ld.r2i)] # column-wise
       # index vectors for LD matrix
       rowi <- rep(1:p,times=(p:1)-1)
       coli <- p+1 - sequence(1:(p-1))
       coli <- coli[length(coli):1]

       # distance between markers
       disti <- abs(posi[rowi] - pos[coli])

       # create dataset with information from above
       ret[[i]] <- data.frame(marker1=mn[rowi],marker2=mn[coli],r2=ld.r2i,dist=disti)

       # create plots

       # scatterplot
       if(type=="p") plot(r2~dist,data=ret[[i]],main=paste("Linkage Group",lg[i]),...)

       # scatterplot with nls curve
       if(type=="nls"){
               plot(r2~dist,data=ret[[i]],main=paste("Linkage Group",lg[i]),...) 
               smooth.fit(ret[[i]][,4],ret[[i]][,3],n=nrow(markeri))
       }

       # stacked histogramm
       if(type=="bars"){
          # use default breaks
          if(is.null(breaks)){
             breaks.dist <- seq(from=min(ret[[i]]$dist),to=max(ret[[i]]$dist),length=6)
             breaks.r2 <- seq(from=1,to=0,by=-0.2) 
          }
          # use user- specified breaks
          else{
             breaks.dist <- breaks$dist
             breaks.r2 <- breaks$r2
          }
          cut.dist <- cut(ret[[i]]$dist,breaks=breaks.dist)
          cut.r2 <- cut(ret[[i]]$r2,breaks=breaks.r2)
          
          # create matrix with relative frequencies
          tab.abs <- table(cut.r2,cut.dist)
          colSum <- matrix(rep(colSums(tab.abs),nrow(tab.abs)),nrow=nrow(tab.abs),byrow=TRUE)

          # baplot out of frequency matrix
          barplot((tab.abs/colSum)[nrow(tab.abs):1,],col=grey(1:nrow(tab.abs)/nrow(tab.abs)),space=c(.2),main=paste("Linkage Group",lg[i]),xlim=c(0,ncol(tab.abs)+2.8),...)
          legend(ncol(tab.abs)+1.2,0.95,fill=grey(1:nrow(tab.abs)/nrow(tab.abs))[nrow(tab.abs):1],legend=levels(cut.r2),title="LD (r2)",cex=1)


}

    }
     if(!is.null(file)) dev.off()
    
    invisible(ret)
}



