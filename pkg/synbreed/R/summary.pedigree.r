summary.pedigree <- function(object,...){
     ped <- object 
     n <- nrow(ped) 
     ans <- list(nID=n,nPar1=length(unique(ped$Par1)),nPar2=length(unique(ped$Par1)),nGener=length(unique(ped$gener)),nUnknownParents=sum(ped$Par1==0)+sum(ped$Par2==0))
     if(!is.null(ped$sex)) ans$sex <- c(males=sum(ped$sex),females=sum(1-ped$sex))
     class(ans) <- "summary.pedigree"
     ans
}

print.summary.pedigree <- function(x,...){
    cat("Number of \n")
    cat("\t individuals ",x$nID,"\n")
   if(!is.null(x$sex)){
      cat("\t males : ",x$sex[1],", females : ",x$sex[2],"\n")
      cat("\t Par 1 (sire) ",x$nPar1,"\n")
      cat("\t Par 2 (dam)  ",x$nPar2,"\n")
      }
   else{ 
    cat("\t Par 1       ",x$nPar1,"\n")
    cat("\t Par 2       ",x$nPar2,"\n")
    }
    cat("\t generations ",x$nGener,"\n")
} 



