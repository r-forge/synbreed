simul.phenotype <- function(pedigree,h2,loc=1,repl=1,seed=123,genotypeEff=NULL,locEff=NULL){

        # read information from data
        gener <- pedigree$gener
        N <- nrow(pedigree)
        n <- N*loc*repl
                                       
        
        # initialize data
        Trait <- rep(0,n)
        ID <- rep(pedigree$ID,length.out=n,each=repl)
        Loc <- rep(1:loc,length.out=n,each=N*repl)
                       
        # set values for first generation
        gener1 <- pedigree$ID[gener==min(gener)]

        
        # set starting value
        set.seed(seed)
        
        
        
        # simulate data by
        
        # adding effect of genotype
        if (is.null(genotypeEff)){ 
          genotypeEff <- rep(100,times=n)
        }
        else genotypeEff <- rep(genotypeEff,length.out=n,each=repl)
        
        Trait <-  rnorm(n,mean=rep(genotypeEff),sd=0.2*N)
        
        
         # adding effect of pedigree              
        for (i in 1:n){
            if(ID[i] %in% gener1)    Trait[i] <- Trait[i]
            
            else {
              p1 <- pedigree[pedigree$ID==ID[i],"Par1"]
              p2 <- pedigree[pedigree$ID==ID[i],"Par2"]  
              # mixing infromation: h(0.5*Par1 + 0.5*Par2) + 1-h(N(mean(gener1),sd(gener1))                      
              Trait[i] <-  h2*(0.5 * mean(Trait[ID == p1]) + 0.5 * mean(Trait[ID == p2])) + (1-h2)*rnorm(1,mean(Trait[ID %in% gener1],na.rm=TRUE),sd=sd(Trait[ID %in% gener1],na.rm=TRUE))
            } 
        } 
        
        #  adding effect of location
        Trait <- Trait + rnorm(length(Trait),mean=Loc*2,sd=0.5*loc)
        
        ret <- data.frame(ID=ID,Trait=Trait,Loc=Loc)
        # sort by ID                                         
        ret <- ret[order(ID),]
        class(ret) <- c("phenoData","data.frame")
        return(ret)


}