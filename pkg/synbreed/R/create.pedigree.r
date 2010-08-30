create.pedigree <- function(ID,Par1,Par2,gener=NULL){

     # only unique IDs
     n <- length(ID)
     if(length(unique(ID))!=n) stop("ID is not unique")
 
     # NA for unkonwn pedigree
     Par1[Par1=="0"] <- NA
     Par2[Par2=="0"] <- NA
     Par1[Par1==0] <- NA
     Par2[Par2==0] <- NA
     
     #Par1[is.na(Par1)] <- 0
     #Par1[is.na(Par2)] <- 0
     
     # add gener if NULL
     if(is.null(gener)){  
          gener <- rep(0,n)  
          for(i in 1:n){
            if(is.na(Par1[i]) & is.na(Par2[i])) gener[i] <- 0 
            else {
             
             gener[i] <- max(c(gener[ID==Par1[i]],gener[ID==Par2[i]]))+1
             }
          }

     }
 
     
     # gener starts from 0
     #if(min(gener) != 0) gener <- gener-min(gener)
     #pedigree <- cbind(ID,Par1,Par2,gener)
     pedigree <- data.frame(ID=ID,Par1=Par1,Par2=Par2,gener=gener,stringsAsFactors=FALSE)
     # sort by generation
     pedigree <- pedigree[order(gener,partial=ID),]
     colnames(pedigree) <- c("ID","Par1","Par2","gener")
     class(pedigree) <- c("pedigree","data.frame")
     
     # missing values as 0
     pedigree[is.na(pedigree)] <- 0
     
     
     
     return(pedigree)
}
