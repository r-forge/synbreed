# read genomic prediction data
create.gpData <- function(pheno=NULL,geno=NULL,map=NULL,pedigree=NULL,family=NULL,covar=NULL,reorderMap=TRUE,map.unit="cM"){
  
  # start with some checks on data
  # geno as matrix but not data.frame (storage) 
  if(!is.null(geno)){
     if(class(geno) == "data.frame"){
        geno <- matrix(unlist(geno),nrow=nrow(geno),ncol=ncol(geno),dimnames=dimnames(geno))
     }
  }

  # match geno and map
  if(!is.null(geno) & !is.null(map)){
    if(ncol(geno) != nrow(map)){  # different subsets of markers in geno and map
       if(is.null(colnames(geno)) & is.null(rownames(map))) stop("no marker names found in 'map' or 'geno'")
       else{ 
        # fill in gaps in map
        warning("not all markers in 'geno' mapped in 'map', gaps filled with 'NA' \n")
        map <- data.frame(chr=map$chr[match(colnames(geno),rownames(map))],pos=map$pos[match(colnames(geno),rownames(map))],row.names=colnames(geno))
       }
    }
    else{  # same subsets of markers in geno and map
      # missing colnames in geno
      if(is.null(colnames(geno))){  # no marker names in geno
        if (is.null(map)) stop("missing rownames in 'geno'")
        if(!is.null(rownames(map))){
          colnames(geno) <- rownames(map)
          warning("missing colnames in 'geno': assuming to be identical as rownames in 'map' \n")
        }
        else{ # default marker names
          warning("no marker names provide in 'geno' or 'map', using default names M1, M2, ... \n")
          colnames(geno) <- rownames(map) <- paste("M",1:ncol(geno),sep="")
        }
      }
      else{  # marker names in geno
      # missing rownames in map
         if(is.null(rownames(map))){
           rownames(map) <- colnames(geno)
           warning("missing rownames in 'map': assuming to be identical as colnames in 'geno' \n")  
         }
        # if(is.null(rownames(map)) & is.null(colnames(geno))){
        #  warning("missing marker names, setting default names M1, M2, ... ")
        #  colnames(geno) <- rownames(map) <- paste("M",1:ncol(geno),sep="")
        #} 
      }
    }
  }

  # match geno and pheno
  if(!is.null(geno) & !is.null(pheno)){
    if(is.null(rownames(pheno)) | is.null(rownames(geno))){
      if(nrow(pheno) == nrow(geno)){
        warning("assuming identical order of genotypes in 'pheno' and 'geno' \n")  
        if(is.null(rownames(pheno))) rownames(pheno) <- rownames(geno)
        else rownames(geno) <- rownames(pheno)
        if(is.null(rownames(pheno)) & is.null(rownames(geno))) rownames(pheno) <- rownames(geno) <- paste("ID",1:nrow(geno),sep="")
    }
    # now geno and pheno have rownames
    else stop("missing rownames for 'pheno' and 'geno'")
    }
   } 
    # sort geno and pheno by rownames (alphabetical order)
    if(!is.null(geno)) geno <- geno[order(row.names(geno)),]
    if(!is.null(pheno)){
      phenonames <- colnames(pheno)
      phenodf <- pheno[order(row.names(pheno)),]
      pheno <- data.frame(phenodf,row.names= row.names(pheno)[order(row.names(pheno))])
      colnames(pheno) <- phenonames
    }

  # sort markers by chromosome and position within chromosome
  if(!is.null(map)){
    if(any(colnames(map) != c("chr","pos"))) stop("colnames of 'map' must be 'chr' and 'pos'")
    if (reorderMap){
     # first order in alphabetical oder (important for SNPs with the same position)
     map <- map[order(as.character(rownames(map))),]
     map <- orderBy(~chr+pos,data=map)
      # sortcolumns in geno, too
     geno <- geno[,rownames(map)]
    }
  }

  # return object
  # geno as matrix
  if(!is.null(geno)) geno <- data.matrix(geno,TRUE)
  obj <- list(covar=NULL,pheno=pheno,geno=geno,map=map,pedigree=pedigree)
  
  # add information to element covar
  # sort all available individuals
  ids <- sort(unique(c(row.names(obj$pheno),rownames(obj$geno),obj$pedigree$ID))) 

  if(is.null(covar)) obj$covar <- data.frame(id=ids,stringsAsFactors=FALSE)
  else obj$covar$id <- ids 

  obj$covar$phenotyped <- obj$covar$id %in% rownames(obj$pheno)
  obj$covar$genotyped <- obj$covar$id %in% rownames(obj$geno)
  
  # family information for genotyped indviduals  
  if(!is.null(family)){
    colnames(family)[1] <- "family"
    obj$covar <- merge(obj$covar,family,by.x=1,by.y=0,all=TRUE)
  }
  
  # add covar from arguments, if available 
  if(!is.null(covar)){
    if(is.null(rownames(covar))) stop("missing rownames in covar")    
    # do not use any existing columns named 'genotyped', 'phenotyped' or 'id'
    covar <- covar[!colnames(covar) %in% c("genotyped","phenotyped","id")]
    # merge with existing data
    if(!is.null(covar)) obj$covar <- merge(obj$covar,covar,by.x=1,by.y=0,all=TRUE)
    else  obj$covar <- obj$covar
  }
  
  # further information
  obj$info$map.unit <- map.unit
  obj$info$codeGeno <- FALSE

  # set class of sub-object pedigree
  if(!is.null(obj$pedigree)) class(obj$pedigree) <- c("pedigree","data.frame")

  # return object of class 'gpData'
  class(obj) <- "gpData"
  return(obj)
}
