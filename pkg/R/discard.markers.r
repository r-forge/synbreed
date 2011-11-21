# discard markers from class 'gpData'

discard.markers <- function(gpData,which){
  # update geno
  gpData$geno <- gpData$geno[,!colnames(gpData$geno) %in% which]
  # update map                 
  gpData$map <- gpData$map[!rownames(gpData$map) %in% which,]
  return(gpData)
}
                     
# discard individuals from class 'gpData'

discard.individuals <- function(gpData,which){
   # updata pheno
   gpData$pheno <- gpData$pheno[!rownames(gpData$pheno) %in% which, , ]
   # update geno
   gpData$geno <- subset(gpData$geno,!rownames(gpData$geno) %in% which) 
   # update pedigree 
   gpData$pedigree <- subset(gpData$pedigree,!gpData$pedigree$ID %in% which) 
   # update covar
   gpData$covar <- subset(gpData$covar,!gpData$covar$id %in% which)
   return(gpData)
}