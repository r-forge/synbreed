create.pedigree <- function(ID,Par1,Par2,gener=NULL){

     # only unique IDs
     n <- length(ID)
     if(length(unique(ID))!=n) step("ID is not unique")
 
     
     # add gener if NULL
     if(is.null(gener)){  
          gener <- rep(0,n)  
          for(i in 1:n){
            if(all(c(Par1[i],Par2[i]) == 0,na.rm=TRUE) | (is.na(Par1[i]) & is.na(Par2[i])) ) a <- 1# gener[i] <- 0 
             else {
             
             gener[i] <- max(gener[ID==Par1[i]],gener[ID==Par2[i]])+1
             }
          }

     }
     # value zero for unkonwn pedigree
     Par1[is.na(Par1)] <- 0
     Par1[is.na(Par2)] <- 0
     
     # gener starts from 0
     if(min(gener) != 0) gener <- gener-min(gener)
     #pedigree <- cbind(ID,Par1,Par2,gener)
     pedigree <- data.frame(ID=ID,Par1=Par1,Par2=Par2,gener=gener,stringsAsFactors=FALSE)
     colnames(pedigree) <- c("ID","Par1","Par2","gener")
     class(pedigree) <- c("pedigree","data.frame")
     
     # missing values as 0
     pedigree[is.na(pedigree)] <- 0
     
     
     
     return(pedigree)
}


ped <- read.table("M:\\Projekte\\SYNBREED\\Research\\pre_synbreed\\Goettingen\\reducedPedigree-2.txt",head=TRUE,stringsAsFactors=FALSE,sep="",strip.white=TRUE)
ped[ped=="   0"] <- NA
ped2 <- create.pedigree(ped$ID,ped$Par1,ped$Par2)