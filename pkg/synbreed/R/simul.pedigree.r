simul.pedigree <- function(generations=2,ids=4){

        # if only one value is given, set samenumber for all generations
        if(length(ids)==1) ids <- rep(ids,times=generations)

        # initialisation
        gener <- rep(1:generations,times=ids)
        ID <- 1:sum(ids)
        Par1 <- Par2 <- rep(0,length(ID))
            
        # random mating any inbreeding lines, parents of F2, will only be in genotopes of F1
        for (i in 2:generations){
          Par1[gener==i] <- ID[sample(ID[gener==i-1],size=ids[i],replace=TRUE)]
          Par2[gener==i] <- ID[sample(ID[gener==i-1],size=ids[i],replace=TRUE)]
         }
        
        # return results
        ped <- data.frame(ID=ID,Par1=Par1,Par2=Par2,gener=gener) 
        class(ped) <- c("pedigree","data.frame") 
        return(ped)
}
