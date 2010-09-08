simul.phenotype <- function(pedigree=NULL,A=NULL,mu=100,vc=NULL,Nloc=1,Nrepl=1,seed=123){

        if(is.null(A) & is.null(pedigree)) step("Either 'pedigree' or 'A' must be given")
        # create A matrix if missing
        if (is.null(A)) A <- kinship(pedigree,ret="add")
       
        
        # read data out of arguments
        N <- nrow(A)         
        n <- N*Nloc*Nrepl
        sigmae <- sqrt(vc$sigma2e)
        sigmaa <- sqrt(vc$sigma2a)
        sigmal <- sqrt(vc$sigma2l)
        sigmab <- sqrt(vc$sigma2b)
        namesA <- rownames(A) 
        if(is.null(rownames(A))) namesA <- 1:N  
                                     
        # initialize data

        ID <- rep(namesA,each=Nrepl*Nloc)
        Loc <- rep(1:Nloc,length.out=n,each=Nrepl)
        Block <- rep(1:Nrepl,lengt.out=n)
        
         # as matrix
        A <- matrix(A,nrow=N)  
                       
        # set starting value for simulation
        set.seed(seed)
        
        # simulate data for contribution parts of phenotypic values
        
        # true breeding values
        #tbv <- rmvnorm(1,rep(0,N),A*sigmaa^2)
        tbv <- chol(A)%*%rnorm(N,0,sigmaa)
        tbv <- rep(tbv,each=Nloc*Nrepl)
        # location effect
        locEff <- rnorm(Nloc,0,sigmal)
        locEff <- rep(locEff,length.out=n,each=Nrepl)
        # block effects 
        blockEff <- rnorm(Nrepl,0,sigmab)
        blockEff <- rep(blockEff,length.out=n)
        # residual
        residual <- rnorm(n,0,sigmae)
        
        # simulate phenotypic value
        Trait <- mu + tbv + locEff + blockEff + residual

        # combine to a data.frame
        ret <- data.frame(ID=factor(ID),Loc=factor(Loc),Block=factor(Block),Trait=Trait,TBV=tbv)
        return(ret)


}

