pairwiseLD <- function(gpData,chr=NULL,type=c("data.frame","matrix","both"),LD.measure="r2",rm.unmapped=TRUE){

    # catch errors

    if(is.null(gpData$geno)) stop("no genotypic data available")
    if(!gpData$info$codeGeno) stop("use function 'codeGeno' before")
    if(is.null(gpData$map))  stop("no map information available")
    type <- match.arg(type)

    # extract information from gpData  (only use mapped markers)
    if(rm.unmapped){
      mapped <- !(is.na(gpData$map$chr) & is.na(gpData$map$pos)) 
    }
    else mapped  <- rep(TRUE,ncol(gpData$geno))
    marker <- gpData$geno[,mapped]
    linkageGroup <- gpData$map$chr[mapped]
    pos <- gpData$map$pos[mapped]
    
    # select chromosomes if 'chr' is specified
    lg <- unique(linkageGroup)
    if(!is.null(chr)){ 
        lg <- chr
        if(any(chr=="all")) linkageGroup <- rep("all",length(linkageGroup))
    }
    
    # initialize return LD value data.frame list
    retList <- list()
    # initialize return LD value matrix list
    retMat <- list()

  
    # loop over all chromosomes (linkage groups)
    for (i in 1:length(lg)){
       
       # read information from data
       markeri <- marker[,linkageGroup==lg[i]]
       p <- ncol(markeri)
       mn <- colnames(markeri)
       posi <- pos[linkageGroup==lg[i]]
       
       # compute LD as R2 (missing values are allowed)
       ld.r2 <- cor(markeri,method="spearman",use="pairwise.complete.obs")^2
       colnames(ld.r2) <- rownames(ld.r2) <- colnames(marker)[linkageGroup == lg[i]]
       ld.r2i <- ld.r2[lower.tri(ld.r2)] # column-wise
       
       # index vectors for LD matrix
       rowi <- rep(1:p,times=(p:1)-1)
       coli <- p+1 - sequence(1:(p-1))
       coli <- coli[length(coli):1]

       # distance between markers
       disti <- abs(posi[rowi] - posi[coli])
       # matrix of distances 
       distance <- as.matrix(dist(pos[linkageGroup == lg[i]],diag=FALSE,upper=FALSE))
       colnames(distance) <- rownames(distance) <- colnames(marker)[linkageGroup == lg[i]]

       # create dataset with information from above in a data.frame
       retList[[lg[i]]] <- data.frame(marker1=mn[rowi],marker2=mn[coli],r2=ld.r2i,dist=disti)
       # and as a matrix
       retMat$LD[[lg[i]]] <- ld.r2           # omit lower/upper triangle?
       retMat$distance[[lg[i]]] <- distance  
    }   


      # return values 
      if(type=="both") return(list(dataFrame=retList,matrix=retMat))
      if(type=="data.frame") return(retList)
      if(type=="matrix") return(retMat)
}