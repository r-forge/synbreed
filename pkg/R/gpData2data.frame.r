# conversion of object class 'gpData' to data.frame


gpData2data.frame <- function(gpData,trait=1,onlyPheno=FALSE,all.pheno=FALSE,all.geno=FALSE,...){
     
     # check for class
     if(class(gpData)!="gpData") stop("object '",substitute(gpData),"' not of class 'gpData'") 

      pheno <- abind(gpData$pheno, matrix(1:dim(gpData$pheno)[1], ncol=dim(gpData$pheno)[3], nrow=dim(gpData$pheno)[1]), along=2)
      IDs <- dimnames(pheno)[[1]]
      reps <- dimnames(pheno)[[3]]
      dimnames(pheno)[[2]][dim(pheno)[2]] <- "ID"
      # choose Traits
      if(all(is.numeric(trait))) trait <- (dimnames(pheno)[[2]])[trait]   
      pheno <- array(pheno[, c("ID", trait), ], dim=c(dim(pheno)[1], length(trait)+1, dim(pheno)[3]))
      dimnames(pheno) <- list(IDs, c("ID", trait), reps)
      # look for covariables
      if(!is.null(gpData$phenoCovars)) pheno <- abind(gpData$phenoCovars, pheno, along=2)
      # append column for each replication
      if(dim(pheno)[3]>1){
        pheno <- abind(matrix(rep(1:dim(gpData$pheno)[3], each=dim(gpData$phenoCovars)[1]), ncol=dim(gpData$pheno)[3], nrow=dim(gpData$pheno)[1]), pheno, along=2)
        dimnames(pheno)[[2]][1] <- "repl"
      }
      pheno <- as.data.frame(apply(pheno, 2, cbind))
      pheno$ID <- IDs[pheno$ID]
      if(!onlyPheno){
        geno <- gpData$geno
        # merge genotypic and phenotypic data
        mergeData <- merge(pheno,geno,by.x="ID", by.y="row.names",all.x=all.pheno,all.y=all.geno)
        # omit row.names column
       } else mergeData <- pheno[, c("ID", colnames(pheno)[!colnames(pheno) %in% "ID"])]
     # sort by ID
     for(i in trait)
       mergeData[, i] <- as.numeric(mergeData[, i])
     if(is.null(mergeData$repl))
       mergeData <- orderBy(~ID,data=mergeData)  
     else{
       mergeData$repl <- as.numeric(mergeData$repl)
       mergeData <- orderBy(~ID+repl,data=mergeData) 
     } 
     rownames(mergeData) <- NULL
     return(mergeData)
}          
