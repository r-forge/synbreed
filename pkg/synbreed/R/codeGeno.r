# coding alles in 0 (major) and 2 (minor)

codeGeno <- function(data,impute=FALSE,popStruc=NULL){

  if(class(data)!= "data.frame" & class(data) != "matrix") stop("wrong data format")

     M <- ncol(data)
     n <- nrow(data)

  if(!is.logical(impute)) stop("impute has to be logical")
  if(impute){
    if(is.null(popStruc)) stop("population structure need to impute missing values")
    if(length(popStruc)!=n) stop("population structure must have equal length as obsersvations in data")
    if(any(is.na(popStruc))) stop("no missing values allowd in popStruc")
   }


  
  codeNumeric <- function(x){
    # names of alleles ordered by allele frequency
    alleles <-  names(table(x)[order(table(x),decreasing=TRUE)])
    if (length(alleles)>2) stop("only biallelic marker allowed")
    return((as.numeric(factor(x,levels=alleles))-1)*2)
   }

  res <- apply(data,2,codeNumeric)
  res <- matrix(res,nrow=n,ncol=M)
  res <- data.frame(res)
  colnames(res) <- colnames(data)
  rownames(res) <- rownames(data)
  # coding of SNPs finished

  # impute missing values
  if(impute){
    for (j in 1:M){
        if(sum(!is.na(res[,j]))>0){
        if(j==1) ptm <- proc.time()[3]
        try({
        # compute population structure  as counts
        poptab <- table(popStruc,res[,j])
        rS <- rowSums(poptab)
        #if(any(rS==0)) warning(paste("Note: No data for at least one population in column",j,"\n",sep=" "))

        # and as frequencies
        poptabF <- poptab/rS       

        # continue only if there are missing values
        if(sum(is.na(res[,j]))>0 ){
          # compute otherstatisticx
          major.allele <- unlist(attr(poptab,"dimnames")[[2]][apply(poptab,1,which.max)])

          # look if SNP is segregating  for this population
          polymorph <- apply(poptab,1,length) >1 & ( apply(poptab,1,min) != 0)
         
        
          # count missing vlalues
          nmissfam <- tapply(is.na(res[,j]),popStruc,sum)
          
          # must be a named list
          names(major.allele) <- names(polymorph)
          
          
          # loop over all families
          
         
          
          for ( i in row.names(poptab)[nmissfam>0] ){
            # impute values
            res[is.na(res[,j]) & popStruc == i ,j] <- as.numeric(ifelse(polymorph[i],sample(c(0,2),size=nmissfam[i],prob=poptabF[i,],replace=TRUE),rep(major.allele[i],nmissfam[i])))
          }
        
        }
        
        if(j==1) cat("Approximative run time ",(proc.time()[3] - ptm)*M," seconds \n",sep=" ")
     })
    }
    }
  }  
  # code again if allele frequeny chanched to to imputing
  if(any(colMeans(res,na.rm=TRUE)>1)) res[,which(colMeans(res,na.rm=TRUE)>1)] <- 2 - res[,which(colMeans(res,na.rm=TRUE)>1)] 
  res <- as.matrix(res,nrow=nrow(x))
  return(res)
}
