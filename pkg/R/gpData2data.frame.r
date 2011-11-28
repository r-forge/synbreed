# conversion of object class 'gpData' to data.frame

gpData2data.frame <- function(gpData,trait=1,onlyPheno=FALSE,all.pheno=FALSE,all.geno=FALSE,repl=NULL,...){
     
      # check for class
      if(class(gpData)!="gpData") stop("object '",substitute(gpData),"' not of class 'gpData'") 
      pheno <- abind(gpData$pheno, matrix(1:dim(gpData$pheno)[1]+10**ceiling(log10(dim(gpData$pheno)[1])), ncol=dim(gpData$pheno)[3], nrow=dim(gpData$pheno)[1]), along=2)
      IDs <- dimnames(pheno)[[1]]
      reps <- dimnames(pheno)[[3]]
      dimnames(pheno)[[2]][dim(pheno)[2]] <- "ID"
      # choose Traits
      if(all(is.numeric(trait))) trait <- (dimnames(pheno)[[2]])[trait]   
      pheno <- array(pheno[, c("ID", trait), ], dim=c(dim(pheno)[1], length(trait)+1, dim(pheno)[3]))
      dimnames(pheno) <- list(IDs, c("ID", trait), reps)
      # look for covariables
      if(!is.null(gpData$phenoCovars)){
        pheno <- abind(gpData$phenoCovars, pheno, along=2)
      }
      # append column for each replication
      if(dim(pheno)[3]>1){
        pheno <- abind(matrix(rep(reps, each=dim(gpData$phenoCovars)[1]), ncol=dim(gpData$pheno)[3], nrow=dim(gpData$pheno)[1]), pheno, along=2)
        dimnames(pheno)[[2]][1] <- "repl"
      }
      if(!is.null(repl)){ 
        if(all(repl %in% 1:dim(pheno)[3])) repl <- dimnames(pheno)[[3]][repl] 
        else if(!all(repl %in% dimnames(pheno)[[3]])) stop("wrong replication names used")
      } else repl <- dimnames(pheno)[[3]]
      pheno <- as.data.frame(apply(pheno[, ,repl], 2, cbind))
      for(i in names(gpData$info$attrPhenoCovars)){
        if(gpData$info$attrPhenoCovars[i] == "numeric") 
          pheno[, i] <- as(as.character(pheno[, i]), gpData$info$attrPhenoCovars[i])
      }
      if(!is.null(repl))
        for(i in colnames(pheno))
          if(class(pheno[, i]) == "factor")
            pheno[, i] <- as.factor(as.character(pheno[, i]))
      pheno$ID <- IDs[as.numeric(as.factor(as.character(pheno$ID)))]
      if(!onlyPheno){
        geno <- gpData$geno
        # merge genotypic and phenotypic data
        mergeData <- merge(pheno,geno,by.x="ID", by.y="row.names",all.x=all.pheno,all.y=all.geno)
        # omit row.names column
       } else mergeData <- pheno[, c("ID", colnames(pheno)[!colnames(pheno) %in% "ID"])]
     # sort by ID
     for(i in trait)
       mergeData[, i] <- as.numeric(mergeData[, i])
     if(all(mergeData$ID %in% 1:nrow(mergeData))) mergeData$ID <- as.numeric(mergeData$ID)
     if(is.null(mergeData$repl))
       mergeData <- orderBy(~ID,data=mergeData)  
     else{
       if(all(mergeData$repl %in% 1:dim(gpData$pheno)[3])){
         mergeData$repl <- as.numeric(mergeData$repl)
       }
       mergeData <- orderBy(~ID+repl,data=mergeData) 
     } 
     mergeData$ID <- as.factor(mergeData$ID)
     rownames(mergeData) <- NULL
     return(mergeData)
}          
