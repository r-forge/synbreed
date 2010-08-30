# coding alles in 0 (major) and 2 (minor)

codeGeno <- function(data,impute=FALSE,popStruc=NULL,maf=NULL,nmiss=NULL,replace.value=NULL){

  if(class(data)!= "data.frame" & class(data) != "matrix") stop("wrong data format")
  if(impute & !is.null(replace.value)) stop("No replacing of values in case of imputing")

    # number of genotypes
    n <- nrow(data)
   
  if(!is.logical(impute)) stop("impute has to be logical")
  if(impute & !is.null(popStruc)){
    if(length(popStruc)!=n) stop("population structure must have equal length as obsersvations in data")
    if(any(is.na(popStruc))) stop("no missing values allowd in popStruc")
   }

   # remove makrers with more than nmiss fraction of missing values
   if(!is.null(nmiss)) dataRaw <- data[,apply(is.na(data),2,mean,na.rm=TRUE)<nmiss]
   else dataRaw <- data
   
     M <- ncol(dataRaw)
     

  
  codeNumeric <- function(x){
    # names of alleles ordered by allele frequency
    alleles <-  names(table(x)[order(table(x),decreasing=TRUE)])
    if (length(alleles)>2) stop("only biallelic marker allowed")
    return((as.numeric(factor(x,levels=alleles))-1)*2)
   }

  res <- apply(dataRaw,2,codeNumeric)
  res <- matrix(res,nrow=n,ncol=M)
  res <- data.frame(res)
  colnames(res) <- colnames(dataRaw)
  rownames(res) <- rownames(dataRaw)
  # coding of SNPs finished
  
  # number of missin values
  nmv <- sum(is.na(res))
   
   # initialize counter  
   cnt1 <- 0    # for nr. of imputations with family structure
   cnt2 <- 0    # for nr. of random imputations

  # if no  imputation should be performed, replace missing values
  # accorind to specified value
  if(!impute & !is.null(replace.value)){
     res[is.na(res)] <- replace.value
  } 

  # impute missing values according to population structure
  if(impute & !is.null(popStruc)){
  
    
  
    # loop over all markers
    for (j in 1:M){
        if(sum(!is.na(res[,j]))>0){
        if(j==1) ptm <- proc.time()[3]
        try({
        # compute population structure  as counts
        poptab <- table(popStruc,res[,j])
        rS <- rowSums(poptab)
        #if(any(rS==0)) warning(paste("Note: No data for at least one population in column",j,"\n",sep=" "))

        # and as frequencies
        #poptabF <- poptab/rS       

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
            res[is.na(res[,j]) & popStruc == i ,j] <- as.numeric(ifelse(polymorph[i],sample(c(0,2),size=nmissfam[i],prob=c(0.5,0.5),replace=TRUE),rep(major.allele[i],nmissfam[i])))
            # update counter
            ifelse(polymorph[i],cnt2 <- cnt2 + nmissfam[i],cnt1 <- cnt1 + nmissfam[i])  
             
          }
        
        }
        
        if(j==1) cat("approximative run time ",(proc.time()[3] - ptm)*M," seconds \n ... \n",sep=" ")
     })
    }
    }
  }  
  
   # impute missing values with no population structure
    if(impute & is.null(popStruc)){
        for (j in 1:M){
             cnt2 <- cnt2 + sum(is.na(res[,j]))
             if(j==1) ptm <- proc.time()[3]
              if (length(table(res[,j]))==1)   res[is.na(res[,j]),j] <- 0
              else res[is.na(res[,j]),j] <- sample(c(0,2),size=sum(is.na(res[,j])),prob=table(res[,j])/n,replace=TRUE)
             if(j==1) cat("Approximative run time ",(proc.time()[3] - ptm)*M," seconds \n",sep=" ")
        }   
   }    
  
  
  # code again if allele frequeny changed to to imputing
  if(impute){
    if(any(colMeans(res,na.rm=TRUE)>1)){
      res[,which(colMeans(res,na.rm=TRUE)>1)] <- 2 - res[,which(colMeans(res,na.rm=TRUE)>1)] 
      res <- as.matrix(res,nrow=nrow(x))
    }
  }
  # remove markers with minor allele frequency < maf
  if(!is.null(maf)) res <- res[,apply(res,2,mean,na.rm=TRUE)>2*maf]

  
  # print summary of imputation
  if(impute){
    cat(paste("Total number of missing values                :",nmv,"\n"))
    cat(paste("Number of imputations by family structure     :",cnt1,"\n"))
    cat(paste("Number of random imputations                  :",cnt2,"\n"))
    cat(paste("Approximative fraction of correct imputations :",round((cnt1+0.5*cnt2)/(cnt1+cnt2),3),"\n"))
  }
  
  return(res)
}
