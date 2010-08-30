kinship <- function(ped,DH=NULL,ret=c("add","kin","dom","gam")){

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
        if(par1gam > 0){
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

      # loop over individuals
      for(i in 1:n){
         ka <- (i-1)*2 + 1
         for(j in i:n){
            kb <- (j-1)*2 + 1
            fab <- 0.5*(G[ka,kb]+G[ka,kb+1]+G[ka+1,kb]+G[ka+1,kb+1])
            A[i,j] <- A[j,i] <- fab
            dab <- G[ka,kb]*G[ka+1,kb+1] + G[ka+1,kb]*G[ka,kb+1]
            D[i,j] <- D[j,i] <- dab
        }
      } # end of loop over individuals
    }
    
    # set return matrices
    if(what == "add") kmat <- A
    if(what == "dom") kmat <- D
    if(what == "kin") kmat <- A/2
    if(what == "gam") kmat <- G

    class(kmat) <- "relationshipMatrix"

    return(kmat)
}

