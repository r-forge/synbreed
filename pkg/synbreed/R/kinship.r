kinship <- function(ped,DH=NULL,ret=c("add","kin","dom","gam")){

    # check for 'gpData'
    if(class(ped)=="gpData"){
      if (!is.null(DH)) DH <- ped$covar[ped$covar$id %in% ped$pedigree$ID ,DH]
      if (is.null(ped$pedigree)) stop("no pedigree found")
      else ped <- ped$pedigree
    }  
    
    
    # number of ids
    n <- nrow(ped)
    if(is.null(DH)) DH <- rep(0,n)
    if(!is.null(DH) & length(DH) != n) stop("DH must have same length as ped")


    # set up extended pedigree
    ID <- rep(seq_along(ped$ID),each=2)
    par1 <- pmatch(ped$Par1,ped$ID,nomatch = 0, duplicates.ok = TRUE)
    par2 <- pmatch(ped$Par2,ped$ID,nomatch = 0, duplicates.ok = TRUE)

    # set up gametic pedigree data.frame
    gamMat <- matrix(data=0,nrow=n*2,ncol=3,byrow=FALSE)
    gamMat[,1] <- ID

    # loop over ID
    for (i in 1:n){
        par1gam <- par1[i]
        par2gam <- par2[i]
        j <- (i-1)*2 + 1
        k <- j + 1
        #  parents of male genome contribution
        if(par1gam > 0){
           gamMat[j,2] <- (par1gam - 1)*2 + 1
           gamMat[j,3] <- (par1gam - 1)*2 + 2
        }
        #  parents of female genome contribution
        if(par2gam > 0){
           gamMat[k,2] <- (par2gam - 1)*2 + 1
           gamMat[k,3] <- (par2gam - 1)*2 + 2
        }
    }  # end of loop over ID

    #  Build Gametic Relationship
    ngam <- 2*n
    DHgam <- rep(DH,each=2)
    G <- diag(ngam)
    dimnames(G) <- list(paste(rep(ped$ID,each=2),rep(1:2,times=n),sep="_"), paste(rep(ped$ID,each=2),rep(1:2,times=n),sep="_"))

    # set inbreed coefficients of DHs on 1
    G[cbind((1:ngam)*DHgam,((1:ngam)+c(1,-1))*DHgam)] <- 1

    # caluclate gametic relationship
    # loop over gamets
    for(i in 1:(ngam-1-DHgam[2*n])){
      ip <- i+1 + (DHgam* rep(c(1,0),ngam))[i]

      for(j in ip:ngam){
          if(gamMat[j,2] > 0) {

            x <- 0.5*(G[i,gamMat[j,2]]+G[i,gamMat[j,3]])

            G[i,j] <- G[j,i] <- x
            }
    }
    } # end of loop over gamets

    # prepare return value
    what <- match.arg(ret,choices=c("add","kin","dom","gam"),several.ok = FALSE)

    # calculate addiditive and dominance relationship
    if(any(c("add","dom","kin") %in% what)){
      A <- D <- matrix(data=NA,nrow=n,ncol=n)
      dimnames(A) <-  dimnames(D) <- list(ped$ID, ped$ID)

     
      # initialize inbreeding coefficients

   # loop over individuals
      for(i in 1:n){
         ka <- (i-1)*2 + 1
         for(j in i:n){
            kb <- (j-1)*2 + 1
            fab <- 0.25*(G[ka,kb]+G[ka,kb+1]+G[ka+1,kb]+G[ka+1,kb+1])
            # set up numerator relationship matrix
            # (Oakley et al, 2007)
            # A[j,k] = 1 + F_j, j=k
            # A[j,k] = 2f_jk  , j!=k
            A[i,j] <- A[j,i] <- 2*fab
       
            # set up dominance relationship matrix
            # pedigree (Y,Z) -> j, (U,V) -> k
            # (Oakley et al, 2007)
            # D[j,k] = 1-F_j, j=k
            # D[j,k] = (f_YU f_ZV + f_YV f_ZU (1-F_j)(1-F_k), j!=k 
            dab <- (G[ka,kb]*G[ka+1,kb+1] + G[ka+1,kb]*G[ka,kb+1])#*(1-G[ka,ka+1])*(1-G[kb,kb+1])
            # acoount for inbreeding
            # dominance = 0 if Fi=1
            D[i,j] <- D[j,i] <- dab
            
        }
      } # end of loop over individuals

      diag(D) <- 1 - (diag(A)-1)
      
    }  # end of if

    # set return matrices
    if(what == "add") kmat <- A
    if(what == "dom") kmat <- D
    if(what == "kin") kmat <- A/2
    if(what == "gam") kmat <- G

    class(kmat) <- "relationshipMatrix"

    return(kmat)
}


