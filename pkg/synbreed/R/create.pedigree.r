create.pedigree <- function(ID,Par1,Par2,gener){

     # value zero for unkonwn pedigree
     Par1[is.na(Par1)] <- 0
     Par1[is.na(Par2)] <- 0
     if(min(gener) != 0) gener <- gener-min(gener)
     pedigree <- cbind(ID,Par1,Par2,gener)
     pedigree <- data.frame(pedigree)
     colnames(pedigree) <- c("ID","Par1","Par2","gener")
     class(pedigree) <- c("pedigree","data.frame")
     pedigree[is.na(pedigree)] <- 0
     return(pedigree)
}
