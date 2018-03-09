add.gpData <- function(gpData1, gpData2){
 if(substr(gpData1$info$version, 47, 50)<0.12) stop(paste("Recode ", substitute(gpData1), "! You have used an old version to create/code ", substitute(gpData1), sep=""))
 if(substr(gpData2$info$version, 47, 50)<0.12) stop(paste("Recode ", substitute(gpData2), "! You have used an old version to create/code ", substitute(gpData2), sep=""))
 if(!is.null(gpData1$pheno))
   pheno1 <- gpData2data.frame(gpData1, onlyPheno=TRUE, trait=1:dim(gpData1$pheno)[2], stringsAsFactors=TRUE) else pheno1 <- NULL
 if(!is.null(gpData1$pheno))
   pheno2 <- gpData2data.frame(gpData2, onlyPheno=TRUE, trait=1:dim(gpData1$pheno)[2], stringsAsFactors=TRUE) else pheno2 <- NULL
 if(is.null(pheno1)){
   if(is.null(pheno2)) pheno <- pheno1 else{
     pheno1
   }
 } else {

 }
 return(0) 
}
