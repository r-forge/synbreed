simul.pedigree <- function(generations=2,ids=4,animals=FALSE){

        # if only one value is given, set samenumber for all generations
        if(length(ids)==1) ids <- rep(ids,times=generations)

        # initialisation
        gener <- rep(1:generations,times=ids)
        ID <- 1:sum(ids)
        Par1 <- Par2 <- rep(0,length(ID))
            
        # random mating for plants (inbreeds are likely)
        if(!animals){
          for (i in 2:generations){
            Par1[gener==i] <- ID[sample(ID[gener==i-1],size=ids[i],replace=TRUE)]
            Par2[gener==i] <- ID[sample(ID[gener==i-1],size=ids[i],replace=TRUE)]
          }
          ped <- data.frame(ID=ID,Par1=Par1,Par2=Par2,gener=gener-1) 
        }
        # define sire and dams for animals (no inbreeds)
        else{
           # define sex for 1st generation
           # 0 = female
           # 1 = male
           sex <- rep(0,length(ID))
           sex[gener==1] <- sample(rep(c(0,1),length=sum(gener==1)),sum(gener==1),replace=FALSE) 
           
           for (i in 2:generations){
            sex[gener==i] <- sample(rep(c(0,1),length=sum(gener==i)),sum(gener==i),replace=FALSE)

            Par1[gener==i] <- ID[sample(ID[(gener==i-1) & sex==1],size=ids[i],replace=TRUE)]
            Par2[gener==i] <- ID[sample(ID[(gener==i-1) & sex==0],size=ids[i],replace=TRUE)]
          }
           ped <- data.frame(ID=ID,Par1=Par1,Par2=Par2,gener=gener-1,sex=sex) 
        } 
        
        
        class(ped) <- c("pedigree","data.frame") 
        return(ped)
}
