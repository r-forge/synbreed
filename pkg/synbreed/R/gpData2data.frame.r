# conversion of object class 'gpData' to data.frame
# note : needed modification for replicated trials

gpData2data.frame <- function(gpData,phenoNo=1){
     
     # check for class
     if(class(gpData)!="gpData") stop("object '",substitute(gpData),"' not of class 'gpData'") 

     # discard nongenotyped and non-phenotyped individuals

     geno <- gpData$geno
     pheno <- gpData$pheno[phenoNo]

     # merge genotypic and phenotypic data
     mergeData <- merge(pheno,geno,by="row.names")
     rownames(mergeData) <- mergeData$Row.names
     # omit Row.names column
     mergeData <- mergeData[,-1]
     
     return(mergeData)
}           

