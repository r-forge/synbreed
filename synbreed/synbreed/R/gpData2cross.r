# conversion to class 'cross'
# return F2 intercross

gpData2cross <- function(gpData,...){
     # check for class
     if(class(gpData)!="gpData") stop("object '",substitute(gpData),"' not of class 'gpData'")  
     # check on geno and map
     if(is.null(gpData$geno) | is.null(gpData$map)) stop("'geno' and 'map' needed in",substitute(gpData))
     else{
       # use codeGeno if not yet done
       if(!gpData$info$codeGeno) gpData <- codeGeno(gpData,...)
     
       # only use individuals with genotypes and phenotypes
       genoPheno <- gpData$covar$id[gpData$covar$genotyped & gpData$covar$phenotyped]
       # read information from gpData
       geno <- data.frame(gpData$geno[rownames(gpData$geno) %in% genoPheno,])
       pheno <- data.frame(gpData$pheno[rownames(gpData$pheno) %in% genoPheno,])
       map  <- gpData$map
       n <- nrow(geno)
     }
     # split markers (+pos) and genotypes on chromosomes
     genoList <- split(cbind(rownames(map),map$pos,t(geno)),map$chr)
     # result is a list       
     # function to bring each list element in right format
     addData <- function(x){                       
         ret <- list()
         Nm <- length(x)/(n+2)
         
         # elements of x:
         #  1:Nm: marker names
         #  (Nm+1):(2*Nm): marker positions
         #  rest: genotypes as vector
         
         # add 1 to genotypes
         # coding for F2 intercross: AA=1, AB=2, BB=3
         ret[["data"]] <- matrix(as.numeric(x[-(1:(2*Nm))])+1,nrow=n,ncol=Nm,byrow=TRUE,dimnames=list(NULL,x[1:Nm]))
         ret[["map"]]  <- as.numeric(x[(Nm+1):(2*Nm)])
         names(ret[["map"]]) <- x[1:Nm]
         # this may have to be changed
         class(ret) <- "A"
         ret
     }
     #                  
     # apply function to each list element
     genoList <- lapply(genoList,addData)
     # create object 'cross'
     cross <- list(geno=genoList,pheno=pheno)
     class(cross) <- c("f2","cross")
     cross                   
}