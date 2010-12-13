create.pedigree <- function(ID,Par1,Par2,gener=NULL){

     
     n <- length(ID)
     # only unique IDs
     if(length(unique(ID))!=n)  warning("ID is not unique, removing duplicated individuals")
     
     # further checks  
     if(length(Par1)!=n) stop("Par1 must have same length as ID") 
     if(length(Par2)!=n) stop("Par2 must have same length as ID") 
     if(!any(Par1 %in% ID) & is.null(gener)) stop("gener must be specified if pedigree is not complete")
       
     # NA for unkonwn pedigree
     Par1[Par1=="0"] <- NA
     Par2[Par2=="0"] <- NA
     Par1[Par1==0] <- NA
     Par2[Par2==0] <- NA
     
     # add gener if NULL
     if(is.null(gener)){  
          gener <- rep(0,n)  
          for(i in 1:n){
            if(is.na(Par1[i]) & is.na(Par2[i])) gener[i] <- 0 
            else {
               gener[i] <- max(c(gener[ID==Par1[i]],gener[ID==Par2[i]]),na.rm=TRUE)+1
             }
          }

     }
 
     
     # gener starts from 0
     pedigree <- data.frame(ID=ID,Par1=Par1,Par2=Par2,gener=gener,stringsAsFactors=FALSE)
     
     # removing duplicated entries
     pedigree <- pedigree[!duplicated(pedigree), ]
                                                 
     # sort by generation
     pedigree <- pedigree[order(gener,partial=ID),]
     colnames(pedigree) <- c("ID","Par1","Par2","gener")
     class(pedigree) <- c("pedigree","data.frame")
     
     # missing values as 0
     pedigree[is.na(pedigree)] <- 0
     
     return(pedigree)
}
