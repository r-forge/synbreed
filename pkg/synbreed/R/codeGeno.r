
# coding genotypic data
codeGeno <- function(gpData,impute=FALSE,impute.type=c("fix","random","family"),replace.value=NULL,maf=NULL,nmiss=NULL,label.heter="AB",keep.identical=TRUE,verbose=FALSE){

  # read information from arguments
  
  # check for class 'gpData'
  if(class(gpData)=="gpData"){
    if(is.null(gpData$geno)) stop("no genotypic data available")
    else data <- gpData$geno
    # family information (population structure) for genotypic data
    popStruc <- gpData$covar$family[gpData$covar$genotyped] 
    if(gpData$info$codeGeno & is.null(label.heter)){
       warning("assuming heterozygous genotypes coded as 1. Use 'label.heter' to specify if that is not the case")
       label.heter <- "1"
    } 
  }
  # atm other formats are supported too
  else{
   data <- gpData
   popStruc <- NULL
   gpData$map <- NULL
  }                                
  #  catch errors
  if(class(data)!= "data.frame" & class(data) != "matrix") stop("wrong data format")
  

  # number of genotypes
  n <- nrow(data)

   
   # keep names of data object
   cnames <- colnames(data)
   rnames <- row.names(data)
   dataRaw <- matrix(unlist(data),nrow=n)


  # elements from control list
  #impute <- control$impute
  if(impute)  impute.type <- match.arg(impute.type)
  #maf <- control$maf
  #nmiss <- control$nmiss
  #label.heter <- control$label.heter
  #replace.value <- control$replace.value
  #keep.identical <- control$keep.identical
  #print(keep.identical)
  
  
  # catch errors 
  if (impute){
  if(!is.logical(impute)) stop("impute has to be logical")
  if(impute.type=="fix" & is.null(replace.value)) stop("'replace.value' must be given for impute.type='fix'")
  # imputing with family information
  if(impute.type=="family" & is.null(popStruc)) stop(paste("family information needed, but '",substitute(gpData),"$covar$family' is empty",sep=""))
  if(impute.type=="family" & !is.null(popStruc)){
    if(length(popStruc)!=n) stop("population structure must have equal length as obsersvations in genotypic data")
    if(any(is.na(popStruc))) stop("no missing values allowed for family information")
  }
  }
 
   # preparation of genotypic data
  
   # remove makrers with more than nmiss fraction of missing values
   if(!is.null(nmiss)){
    if(nmiss<0 | nmiss>1) stop("'nmiss' must be in [0,1]")
    which.miss <- apply(is.na(dataRaw),2,mean,na.rm=TRUE)<=nmiss 
    dataRaw <- dataRaw[,which.miss]
    cat("step 1 :",sum(!which.miss),"marker(s) removed with >",nmiss*100,"% missing values \n")
    cnames <- cnames[which.miss]
    # update map
    if(!is.null(gpData$map)) gpData$map <- gpData$map[which.miss,]
    }
    else{
      dataRaw <- dataRaw
      if (verbose) cat("step 1 : No markers removed due to fraction of missing values \n")
    }
    


   # identify heterozygous genotypes
   if(!is.null(label.heter)){
    if (is.character(label.heter)) is.heter <- function(x) return(length(grep(label.heter,as.vector(x)))>0)
    else{
      if (is.function(label.heter))  is.heter <- label.heter
      else stop("label.heter must be a character or a function")
      } 
   
   is.heter <- Vectorize(is.heter,USE.NAMES = FALSE)
   heter.ind <- which(matrix(is.heter(dataRaw),nrow=n),arr.ind=TRUE)
   
   # set values as NA (instead there would be more than two alleles at each locus)
   # this would cause problems in recoding
   # values are restored after recoding
   dataRaw[heter.ind] <- NA
   sumNA.heter <- sum(is.na(dataRaw[heter.ind] ))
   }
   else sumNA.heter <- 0

   
   # function to recode alleles within one locus  
  codeNumeric <- function(x){
    # names of alleles ordered by allele frequency
    alleles <-  names(table(x)[order(table(x),decreasing=TRUE)])
    if (length(alleles)>2) stop("only biallelic marker allowed (check for heterozygous genotypes)")
    return((as.numeric(factor(x,levels=alleles))-1)*2)
   }

  if (verbose) cat("step 2 : Recoding alleles \n")
  res <- apply(dataRaw,2,codeNumeric)
  res <- matrix(res,nrow=n)
  
  # coding of SNPs finished

  # start of imputing
  if(impute){
  
  # number of markers
   M <- ncol(res)          
  
  # start imputing of values
  
   # number of missing values - number of heterozygous values
   nmv <- sum(is.na(dataRaw)) - sumNA.heter
 
  # initialize counter  
   cnt1 <- 0   # for nr. of imputations with family structure
   cnt2 <- 0    # for nr. of random imputations
   

  # if impute.type="fix", replace missing values according to specified value
  if(impute.type=="fix"){  
     res[is.na(res)] <- replace.value
     if (verbose) cat("step 3 : Replace missing values by",replace.value," \n")
  } 

  # impute missing values according to population structure
  if(impute.type=="family"){
   if (verbose) cat("step 3 : Imputing of missing values by population structure \n")
   # initialize counter (- number of heterozygous values) 
   cnt1 <- - sumNA.heter    # for nr. of imputations with family structure
   cnt2 <- - sumNA.heter    # for nr. of random imputations
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
        
        if(j==ceiling(M/100)) cat("... approximative run time ",(proc.time()[3] - ptm)*99," seconds ... \n",sep=" ")
     })
    }
    }
  }  

   # impute missing values with no population structure
    if(impute.type=="random"){
       # initialize counter (- number of heterozygous values) 
       cnt2 <- - sumNA.heter    # for nr. of random imputations
        for (j in 1:M){
             cnt2 <- cnt2 + sum(is.na(res[,j]))
             if(j==1) ptm <- proc.time()[3]
              if (length(table(res[,j]))==1)   res[is.na(res[,j]),j] <- 0
              else res[is.na(res[,j]),j] <- as.numeric(sample(names(table(res[,j])),size=sum(is.na(res[,j])),prob=table(res[,j])/n,replace=TRUE))
             if(j==1) cat("approximate run time ",(proc.time()[3] - ptm)*M," seconds \n",sep=" ")
        }   
   }    
  
  
  # code again if allele frequeny changed to to imputing
    if(any(colMeans(res,na.rm=TRUE)>1)){
      if (verbose) cat("step 4 : Recode alleles due to imputation \n")
      res[,which(colMeans(res,na.rm=TRUE)>1)] <- 2 - res[,which(colMeans(res,na.rm=TRUE)>1)]     
    }
    else{
     if (verbose) cat("step 4 : No recoding of alleles necessary after imputation \n") 
    }
  
  # end of imputing
  }
  
  # restore heterozygous values
  if(!is.null(label.heter)) res[heter.ind] <- 1
  
  
  # remove markers with minor allele frequency < maf
  if(!is.null(maf)){
    if(maf<0 | maf>1) stop("'maf' must be in [0,1]")
    which.maf <- apply(res,2,mean,na.rm=TRUE)>=2*maf 
    if (verbose) cat("step 5 :",sum(!which.maf),"marker(s) removed with maf <",maf,"\n")
    res <- res[,which.maf]
    cnames <- cnames[which.maf] 
    # update map
    if(!is.null(gpData$map)) gpData$map <- gpData$map[which.maf,]
  }
  else {
   if (verbose) cat("step 5 : No markers discarded due to minor allele frequency \n")
  }
  # discard duplicated markers
  if(!keep.identical){
       which.duplicated <- duplicated(res,MARGIN=2)
       res <- res[,!which.duplicated]
       cat("step 6 :",sum(which.duplicated),"duplicated marker(s) discarded \n")
       cnames <- cnames[!which.duplicated]
       # update map 
       if(!is.null(gpData$map)) gpData$map <- gpData$map[!which.duplicated,]
  }
  else{
  if (verbose) cat("step 6 : No duplicated markers discarded \n")
     }
     
    # keep structure for return object
   if (verbose) cat("step 7 : Restoring original data format \n")
   if(is.matrix(data)){
     res <- matrix(res,nrow=n)
     rownames(res) <- rnames
     colnames(res) <- cnames
   }
   if(is.data.frame(data)){
      res  <- as.data.frame(res)
      colnames(res) <- cnames
      row.names(res) <- rnames
   }


  # print summary of imputation
  if(impute){
    cat("\n")
    cat("Summary of imputation \n")
    cat(paste("  total number of missing values                :",nmv,"\n"))
    cat(paste("  number of imputations by family structure     :",cnt1,"\n"))
    cat(paste("  number of random imputations                  :",cnt2,"\n"))
    cat(paste("  approximate fraction of correct imputations   :",round((cnt1+0.5*cnt2)/(cnt1+cnt2),3),"\n"))
  }
  
  # overwrite original genotypic data
  if(class(gpData)=="gpData") {
    gpData$geno <- res
    gpData$info$codeGeno <- TRUE
  }
  else gpData <- res

  return(gpData)
}