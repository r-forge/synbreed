# coding alles in 0 (major) and 2 (minor)

codeGeno <- function(data,impute=FALSE,popStruc=NULL,maf=NULL,nmiss=NULL,label.heter=NULL,replace.value=NULL){

  if(class(data)!= "data.frame" & class(data) != "matrix") stop("wrong data format")
  if(impute & !is.null(replace.value)) stop("No replacing of values in case of imputing")

    # number of genotypes
    n <- nrow(data)
    
    # keep names of data object
   if(is.data.frame(data)){
     cnames <- colnames(data)
     rnames <- row.names(data)
     dataRaw <- matrix(unlist(data),nrow=n)
   }
   else dataRaw <- data

  # catch errors 
  if(!is.logical(impute)) stop("impute has to be logical")
  if(impute & !is.null(popStruc)){
    if(length(popStruc)!=n) stop("population structure must have equal length as obsersvations in data")
    if(any(is.na(popStruc))) stop("no missing values allowd in popStruc")
   }

   # remove makrers with more than nmiss fraction of missing values
   if(!is.null(nmiss)){
    which.miss <- apply(is.na(dataRaw),2,mean,na.rm=TRUE)<nmiss 
    dataRaw <- dataRaw[,which.miss]

    if(is.data.frame(data)) cnames <- cnames[which.miss]

    }
   else dataRaw <- dataRaw

   
   # number of markers
   M <- ncol(dataRaw)    

   # identify heterozygous
   if(!is.null(label.heter)){
    if (is.character(label.heter)) is.heter <- function(x) return(length(grep(label.heter,as.vector(x)))>0)
    else{
      if (is.function(label.heter))  is.heter <- label.heter
      else stop("label.heter must be a character or a function")
      } 
   
   is.heter <- Vectorize(is.heter,USE.NAMES = FALSE)
   heter.ind <- which(matrix(is.heter(dataRaw),nrow=n,ncol=M),arr.ind=TRUE)
   
   # set values as NA (instead there would be more than two alleles at each locus)
   # this would cause promblems in recoding
   # values are restored after recoding
   
   dataRaw[heter.ind] <- NA

   }

   
   
   
   # function to redoce alleles within one locus  
  codeNumeric <- function(x){
    # names of alleles ordered by allele frequency
    alleles <-  names(table(x)[order(table(x),decreasing=TRUE)])
    if (length(alleles)>2) stop("only biallelic marker allowed (check for heterozygous genotypes)")
    return((as.numeric(factor(x,levels=alleles))-1)*2)
   }

  res <- apply(dataRaw,2,codeNumeric)
  res <- matrix(res,nrow=n,ncol=M)
  
  # coding of SNPs finished
          
  
  # start imputing of values
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
             if(j==1) cat("approximate run time ",(proc.time()[3] - ptm)*M," seconds \n",sep=" ")
        }   
   }    
  
  
  # code again if allele frequeny changed to to imputing
   if(impute){
    if(any(colMeans(res,na.rm=TRUE)>1)){
      res[,which(colMeans(res,na.rm=TRUE)>1)] <- 2 - res[,which(colMeans(res,na.rm=TRUE)>1)]     
    }
  }
  
  
  # restore heterozygous values
  if(!is.null(label.heter)) res[heter.ind] <- 1
  

 

  
  # remove markers with minor allele frequency < maf
  if(!is.null(maf)){
    which.maf <- apply(res,2,mean,na.rm=TRUE)>2*maf 

    res <- res[,which.maf]
    if(is.data.frame(data)){
      cnames <- cnames[which.maf]
    }
    
    
  }

    # keep structure for return object
   if(is.matrix(data)) res <- as.matrix(res,nrow=n)
   if(is.data.frame(data)){
      res  <- as.data.frame(res)
      colnames(res) <- cnames
      row.names(res) <- rnames
   }


  # print summary of imputation
  if(impute){
    cat(paste("total number of missing values                :",nmv,"\n"))
    cat(paste("number of imputations by family structure     :",cnt1,"\n"))
    cat(paste("number of random imputations                  :",cnt2,"\n"))
    cat(paste("approximate fraction of correct imputations :",round((cnt1+0.5*cnt2)/(cnt1+cnt2),3),"\n"))
  }
  
  return(res)
}
