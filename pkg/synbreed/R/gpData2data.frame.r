# conversion of object class 'gpData' to data.frame


gpData2data.frame <- function(gpData,phenoNo=1,Rep=NULL,onlyPheno=!is.null(Rep),...){
     
     # check for class
     if(class(gpData)!="gpData") stop("object '",substitute(gpData),"' not of class 'gpData'") 

     # for single meausers
     if(is.null(Rep)){
      pheno <- gpData$pheno[phenoNo]
      if(!onlyPheno){
        geno <- gpData$geno
        # merge genotypic and phenotypic data
        mergeData <- merge(pheno,geno,by="row.names")
        rownames(mergeData) <- mergeData$Row.names
        # omit row.names column
        mergeData <- mergeData[,-1]
       }
      else mergeData <- pheno
     }
     
     # for repeated measures
     if(!is.null(Rep)){
       # only 1 variable for multiple measures
       if(length(Rep)>1) warning("only first grouping variable used")
      
      # only use specified variables
      pheno <- gpData$pheno[Rep[[1]]]
         
      # reshape from long to wide data format
      pheno <- reshape(pheno,varying=Rep[[1]],timevar=names(Rep)[1],v.names="trait",direction="long",...)
      
      # merge with genotypic data 
       if(!onlyPheno){
          geno <- gpData$geno
          mergeData <- merge(pheno,geno,by.x="id",by.y="row.names")
          rownames(mergeData) <- rownames(pheno) 
       }
       else  mergeData <- pheno
       
     
     
     }   
     
     return(mergeData)
}           

