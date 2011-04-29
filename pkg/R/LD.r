pairwiseLD <- function(gpData,chr=NULL,type=c("data.frame","matrix"),distances=TRUE,ld.threshold=0,rm.unmapped=TRUE){

    # catch errors
    if(is.null(gpData$geno)) stop("no genotypic data available")
    if(!gpData$info$codeGeno) stop("use function 'codeGeno' before")
    if(is.null(gpData$map))  stop("no map information available")
    
    type <- match.arg(type)
    if(ld.threshold>0 & type=="matrix") warning("'ld.threshold not used for type='matrix'")
    if(any(is.na(gpData$geno))) stop("no missing values allowed, try to impute using 'codeGeno'")

    # extract information from gpData  (only use mapped markers)
    if(rm.unmapped){
      mapped <- !(is.na(gpData$map$chr) & is.na(gpData$map$pos)) 
    }
    else mapped  <- rep(TRUE,ncol(gpData$geno))
    
    # restrict data to mapped markers
    linkageGroup <- gpData$map$chr[mapped]
    pos <- gpData$map$pos[mapped]
    names(pos) <- rownames(gpData$map)[mapped]
   
    
    # select chromosomes if 'chr' is specified
    lg <- unique(linkageGroup)
    if(!is.null(chr)){ 
        lg <- chr
        if(any(chr=="all")) stop("option chr='all' not yet possible") #linkageGroup <- rep("all",length(linkageGroup))
    }                                                
  
    # initialize return LD value data.frame list
    retList <- list()
    # initialize return LD value matrix list
    retMat <- list()

  
    # loop over all chromosomes (linkage groups)
    for (i in 1:length(lg)){
       
       
       # call PLINK to compute the LD as r2
       sel <- rownames(gpData$map[gpData$map$chr!= lg[i],])
       gpTEMP <- discard.markers(gpData,which=sel)
       pre <- paste("chr",lg[i],sep="")
       write.plink(gpTEMP,type=type,ld.threshold=ld.threshold,prefix=pre) 
       system(paste("plink --script ",pre,"plinkScript.txt",sep="")) 
       
       # read data from PLINK
       if(type=="matrix") ld.r2 <- as.matrix(read.table(paste(pre,".ld",sep="")))
       if(type=="data.frame") ld.r2.df <- read.table(paste(pre,".ld",sep=""),header=TRUE,stringsAsFactors=FALSE)
       

       # distances between markers
       if (distances) {
         if(type=="matrix"){
          distance <- as.matrix(dist(pos[linkageGroup == lg[i]],diag=FALSE,upper=FALSE))
          colnames(distance) <- rownames(distance) <- colnames(ld.r2 ) <- rownames(ld.r2 ) <-names(pos)[linkageGroup == lg[i]]
         }
        if(type=="data.frame"){
          distance <- abs(pos[ld.r2.df$SNP_A]-pos[ld.r2.df$SNP_B]) 
        }
       }
       else distance <- rep(NA)
       
       # create dataset with information from above in a data.frame
       if(type=="data.frame") retList[[lg[i]]] <- with(ld.r2.df,data.frame(marker1=SNP_A,marker2=SNP_B,r2=R2,dist=distance,stringsAsFactors=FALSE))
       # and as a matrix
       if(type=="matrix")retMat$LD[[lg[i]]] <- ld.r2           # omit lower/upper triangle?
       if(type=="matrix")retMat$distance[[lg[i]]] <- distance  
    }   


      # return values 
      #if(type=="both") return(list(dataFrame=retList,matrix=retMat))
      if(type=="data.frame") {
        class(retList) <- "LDdf"
        return(retList)
      }  
      if(type=="matrix"){
         class(retMat) <- "LDmat"
         return(retMat)
      }   
}
