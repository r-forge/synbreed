summary.pedigree <- function(pedigree){
     n <- nrow(pedigree) 
     cat("# individuals    :",n,"\n")
     cat("# Par1 (sire)    :",length(unique(pedigree$Par1)),"\n")
     cat("# Par2 (dam)     :",length(unique(pedigree$Par2)),"\n")
     cat("# generations    :",length(unique(pedigree$gener)),"\n")
     cat("# unknow parents :",sum(pedigree$Par1==0)+sum(pedigree$Par2==0),"\n")
     cat("# inbred         :",sum((pedigree$Par1==pedigree$Par2) & pedigree$Par1 != 0 & pedigree$Par2 !=0 ),"\n")
     if(!is.null(pedigree$sex)){
      cat("# male/female    :",n*mean(pedigree$sex),"/",n*(1-mean(pedigree$sex)),"\n")
     }
}    
    