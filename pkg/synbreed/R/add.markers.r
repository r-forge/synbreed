# adding markers to gpData object

add.markers <- function(gpData,geno,map){
      # check if markers are allready in data
      if (any(colnames(geno) %in% c(colnames(gpData$geno),rownames(gpData$map))) | any(rownames(map) %in% c(colnames(gpData$geno),rownames(gpData$map)))){
        stop("some of the markers of ", substitute(geno)," are allready in ", substitute(gpData))
      }
      # merge genotypic data
      geno <- merge(gpData$geno,geno,by.x="row.names",by.y="row.names",sort=FALSE)
      # first column as rownames and delete first column
      rownames(geno) <- geno[,1]
      geno <- geno[,-1]
      # merge map
      map <- rbind(gpData$map,map)
    
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
      pheno <- rbind(gpData$pheno,pheno)
      # merge genotypic data
      geno <- rbind(gpData$geno,geno)
      # merge pedigree
      pedigree <- rbind(gpData$pedigree)
      # merge covar  (not with first columns)
      if(any(colnames($covar) %in% c("genotyped","phenotyped","id")])) stop("not specify columns 'genotyped','phenotyped' and 'id' in 'covar' ")
      covar <- rbind(gpData$covar[,!colnames(gpData$covar) %in% c("genotyped","phenotyped","id")],covar)
      # create new gpData object
      ret <- create.gpData(pheno=pheno,geno=geno,map=gpData$map,pedigree=pedigree,covar=covar,map.unit=gpData$info$map.unit)
      return(ret)
}
