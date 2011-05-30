create.pedigree <- function(ID,Par1,Par2,gener=NULL,sex=NULL,add.ancestors=FALSE){
     n <- length(ID)
     # only unique IDs
     if(length(unique(ID))!=n)  warning("ID is not unique, removing duplicated individuals")
     
     # further checks  
     if(length(Par1)!=n) stop("Par1 must have same length as ID") 
     if(length(Par2)!=n) stop("Par2 must have same length as ID") 
    
     # NA for unkonwn pedigree
     Par1[Par1=="0"] <- NA
     Par2[Par2=="0"] <- NA
     Par1[Par1==0] <- NA
     Par2[Par2==0] <- NA  
     
     # add ancestors (if required)
     if(add.ancestors){
          ancestors <- c(Par1[!is.na(Par1) & !(Par1 %in% ID) ],Par2[!is.na(Par1) & !(Par1 %in% ID) ])
          ID <- c(ancestors,ID)
          Par1 <- c(rep(NA,length(ancestors)),Par1)
          Par2 <- c(rep(NA,length(ancestors)),Par2)
          if (!is.null(gener)) gener <- c(rep(NA,length(ancestors)),gener)
          if (!is.null(sex)) sex <- c(rep(NA,length(ancestors)),sex)
          # update n
          n <- length(ID)

     }
     
     # add gener if NULL
     if(is.null(gener)){  
          #if(any(! Par1 %in% c(0,ID) ) | any(! Par1 %in% c(0,ID) )) stop("parents without pedigree, try to use argument 'add.ancestors=TRUE'")
          gener <- rep(0,n)  
          for(i in 1:n){
            if(is.na(Par1[i]) & is.na(Par2[i])) gener[i] <- 0 
            else{
               gener[i] <- max(c(gener[ID==Par1[i]],gener[ID==Par2[i]]),na.rm=TRUE)+1
             }
          }
     }
     
  
     # gener starts from 0
     if(!is.null(sex)) pedigree <- data.frame(ID=ID,Par1=Par1,Par2=Par2,gener=gener,sex=sex,stringsAsFactors=FALSE)
     else pedigree <- data.frame(ID=ID,Par1=Par1,Par2=Par2,gener=gener,stringsAsFactors=FALSE)
     
     # removing duplicated entries
     pedigree <- pedigree[!duplicated(pedigree), ]
                                                 
     # sort by generation
     pedigree <- pedigree[order(pedigree$gener,partial=pedigree$ID),]
     class(pedigree) <- c("pedigree","data.frame")
     
     # missing values as 0
     pedigree[is.na(pedigree)] <- 0
     
     return(pedigree)
}
