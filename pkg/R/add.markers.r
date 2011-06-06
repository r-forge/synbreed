# adding markers to gpData object

add.markers <- function(gpData,geno,map=NULL){
      # check if markers are allready in data
      if (any(colnames(geno) %in% c(colnames(gpData$geno),rownames(gpData$map))) | any(rownames(map) %in% c(colnames(gpData$geno),rownames(gpData$map)))){
        stop("some of the markers of ", substitute(geno)," are allready in ", substitute(gpData))
      }
      # take names form map if available
      if(is.null(colnames(geno)) & !is.null(map)) colnames(geno) <- rownames(map)
      
      # merge genotypic data
      geno <- merge(gpData$geno,geno,by.x="row.names",by.y="row.names",sort=FALSE)
      # first column as rownames and delete first column
      rownames(geno) <- geno[,1]
      geno <- data.matrix(geno[,-1],TRUE)
      # merge map
      map <- rbind(gpData$map,map)
    
      # need id in rownames for create.gpData
      rownames(gpData$covar) <- gpData$covar$id
      
      # create new gpData object
      ret <- create.gpData(pheno=gpData$pheno,geno=geno,map=map,pedigree=gpData$pedigree,covar=gpData$covar,map.unit=gpData$info$map.unit)
      return(ret)
}


# adding new individuals to gpData object

add.individuals <- function(gpData,pheno=NULL,geno=NULL,pedigree=NULL,covar=NULL){
      # check if individuals are allready in data
      if (any(rownames(pheno) %in% gpData$covar$id) | any(rownames(geno) %in% gpData$covar$id) ){
        stop("some of the individuals of are allready in ", substitute(gpData))
      }
      
      # merge phenotypic data
      if(!is.null(pheno)) if(any(colnames(pheno)!=colnames(gpData$pheno))) stop("different phenotypes (colnames) in 'gpData$pheno' and 'pheno'")
      pheno <- rbind(gpData$pheno,pheno)
      
      # merge genotypic data
      if(!is.null(geno)) if(any(colnames(geno)!=colnames(gpData$geno))) stop("different markers (colnames) in 'gpData$geno' and 'geno'")
      geno <- rbind(gpData$geno,geno)
      
      # merge pedigree
      pedigree <- rbind(gpData$pedigree,pedigree)
      # reorder if necessary
      if(!is.null(pedigree)) pedigree <- orderBy(~gener,pedigree)
      
      # merge covar  (not with first columns)
      
      if(any(colnames(covar) %in% c("genotyped","phenotyped","id"))) stop("not specify columns 'genotyped','phenotyped' and 'id' in 'covar' ")
      if(!is.null(covar)){
        cc <- data.frame(gpData$covar[,!colnames(gpData$covar) %in% c("genotyped","phenotyped","id")])
        rownames(cc) <- rownames(gpData$covar)
        colnames(cc) <- colnames(gpData$covar)[!colnames(gpData$covar) %in% c("genotyped","phenotyped","id")]
        covarUpdate <- rbind(cc,covar) 
      }  
      else covarUpdate <- gpData$covar
      
      # need id in rownames for create.gpData
      rownames(covarUpdate) <- c(as.character(gpData$covar$id),rownames(covar))


      # create new gpData object
      ret <- create.gpData(pheno=pheno,geno=geno,map=gpData$map,pedigree=pedigree,covar=covarUpdate,map.unit=gpData$info$map.unit)
      return(ret)
}
