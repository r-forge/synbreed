# Cross validation with different sampling and variance components estimation methods

crossVal <- function (gpData,trait=1,cov.matrix=NULL, k=2,Rep=1,Seed=NULL,sampling=c("random","within popStruc","across popStruc"), varComp=NULL,popStruc=NULL, VC.est=c("commit","ASReml","BRR","BL"),priorBLR=NULL,verbose=TRUE,nIter=5000,burnIn=1000,thin=10) 
{
    VC.est <- match.arg(VC.est)
    sampling <- match.arg(sampling)
    if(!gpData$info$codeGeno) stop("use function 'codeGeno' before using 'crossVal'") 

    # individuals with genotypes and phenotypes
    dataSet <- as.character(gpData$covar$id[gpData$covar$genotyped & gpData$covar$phenotyped])

    # number of individuals
    n <- length(dataSet)

    # constructing design matrices
    y <- data.frame(rownames(gpData$pheno),gpData$pheno[ ,trait])
    colnames(y) <- c("ID","TRAIT")
    y <- y[y$ID %in% dataSet, ]
    X <- matrix(rep(1,n,ncol=1))
    rownames(X) <- y[,1]
    Z <- diag(n)
    rownames(Z) <- y[,1]

    # checking covariance matrices
    if (is.null(cov.matrix) ){
	Z <- gpData$geno[rownames(gpData$geno) %in% dataSet, ]
	if (VC.est %in% c("commit","ASReml")) cov.matrix <- list(kin=diag(ncol(Z)))
    }

    if(is.null(colnames(X)) & !is.null(colnames(Z))) names.eff <- c(paste("X",1:ncol(X),sep=""),colnames(Z))
    if(is.null(colnames(X)) & is.null(colnames(Z))) names.eff <- c(paste("X",1:ncol(X),sep=""),paste("Z",1:ncol(Z),sep=""))

    # catch errors	
    if(is.null(varComp) & VC.est=="commit") stop("Variance components have to be specified")
    if(VC.est=="commit" & length(varComp)<2) stop("Variance components should be at least two, one for the random effect and one residual variance")
    if(sampling!="random" & is.null(popStruc) & is.null(gpData$covar$family)) stop("no popStruc was given")
    if(sampling!="random" & is.null(popStruc)){
	popStruc <- gpData$covar$family[gpData$covar$genotyped & gpData$covar$phenotyped]
    }
    if(sampling!="random" & !is.null(popStruc)){
      if(length(popStruc)!=n) stop("population structure must have equal length as obsersvations in data")
      if(any(is.na(popStruc))) stop("no missing values allowed in popStruc")
    }
    if ( k < 2) stop("folds should be equal or greater than 2")
    if ( k > n) stop("folds should be equal or less than the number of observations")
    if (VC.est=="commit" & !is.null(cov.matrix) & length(cov.matrix)!=length(varComp)-1) stop("number of variance components does not match given covariance matrices")
    if(VC.est=="BL" & is.null(priorBLR)) stop("prior for varE has to be specified")
    if(VC.est=="BRR" & is.null(priorBLR)) stop("prior for varBR and varE have to be specified")

    # prepare covariance matrices
    if (!is.null(cov.matrix)){
	m <- length(cov.matrix)
	cat("Model with ",m," covariance matrix/ces \n")
	if (VC.est=="commit"){
	   # function for constructing GI
	   rmat<-NULL
   	   for( i in 1:length(cov.matrix)){
       	   m <- solve(as.matrix(cov.matrix[[i]]))*(varComp[length(varComp)]/varComp[i])
       	     if(i==1) rmat <- m
       	     else
             {
              nr <- dim(m)[1]
              nc <- dim(m)[2]
              aa <- cbind(matrix(0,nr,dim(rmat)[2]),m)
              rmat <- cbind(rmat,matrix(0,dim(rmat)[1],nc))
              rmat <- rbind(rmat,aa)
             }
       	   }
         GI <- rmat
         }
	# covariance matrices for ASReml
	else { 
	   if(VC.est=="ASReml"){
      	   for ( i in 1:length(cov.matrix)){
	   write.relationshipMatrix(as.matrix(cov.matrix[[i]]),file=paste("ID",i,".giv",sep=""),type="inv",sorting="ASReml",digits=10)
	   }
	   ID1 <- paste("ID",1:length(cov.matrix),".giv \n",sep="",collapse="")
	   ID2 <- paste("giv(ID,",1:length(cov.matrix),") ",sep="",collapse="")
	   cat(paste("Model \n ID     	  !A \n Yield  	  !D* \n",ID1,"Pheno.txt !skip 1 !AISING !maxit 11\n!MVINCLUDE \n \nYield  ~ mu !r ",ID2,sep=""),file="Model.as")
	   cat("",file="Model.pin")
	  }
	}
    }
    # set seed for replications
    if(!is.null(Seed)) set.seed(Seed)
    seed2<-round(runif(Rep,1,100000),0)

    # begin replications
    COR3 <- NULL
    bu3 <- NULL
    lm.coeff2 <- NULL
    y.TS2 <- NULL
    n.TS2<-NULL
    id.TS2 <- list()   
    for (i in 1:Rep){ 

	# sampling of k sets
	# random sampling
	if(verbose) cat(sampling," sampling \n")
	if(sampling=="random"){
	  y.u <- unique(y[,1])
	  set.seed(seed2[i])
	  modu<-n%%k
	  val.samp2<-sample(c(rep(1:k,each=(n-modu)/k),sample(1:k,modu)),n,replace=FALSE)
	  val.samp3 <- data.frame(y.u,val.samp2)
	  }

	# within family sampling
	if(sampling=="within popStruc"){
	   which.pop <- unique(popStruc)# nr of families
	   y.u <- unique(y[,1])
	   val.samp3<- NULL
	   for (j in 1:length(which.pop)){
		y2<-matrix(y.u[popStruc==which.pop[j]],ncol=1)# select each family
		set.seed(seed2[i]+j) # in each family a different seed is used to result in more equal TS sizes
		modu<-nrow(y2)%%k
		if(!modu==0) val.samp<-sample(c(rep(1:k,each=(nrow(y2)-modu)/k),sample(1:k,modu)),nrow(y2),replace=FALSE)
		if(modu==0) val.samp<-sample(rep(1:k,each=(nrow(y2))/k),nrow(y2),replace=FALSE)
		val.samp2 <- data.frame(y2,val.samp)		
		val.samp3 <- as.data.frame(rbind(val.samp3,val.samp2))
	   }
	   val.samp3 <- val.samp3[order(as.character(val.samp3[,1])),]
	}

	# across family sampling
	if(sampling=="across popStruc"){
	  which.pop <- unique(popStruc)
	  y.u <- unique(y[,1])
	  y2 <- matrix(y.u[order(popStruc)],ncol=1)
	  b <- table(popStruc)
	  modu<-length(which.pop)%%k
	  set.seed(seed2[i])
	  val.samp<-sample(c(rep(1:k,each=(length(which.pop)-modu)/k),sample(1:k,modu)),length(which.pop),replace=FALSE)
	  val.samp2<- rep(val.samp,b)
	  val.samp3 <- data.frame(y2,val.samp2)
	  val.samp3 <- 	as.data.frame(val.samp3[order(as.character(val.samp3[,1])),])
	 #print(head(val.samp3))
	 }

     # start k folds
     COR2 <- NULL
     bu2 <- NULL
     lm.coeff <- NULL
     y.TS <- NULL
     n.TS <- NULL
     id.TS <- list()
     for (ii in 1:k){
	if (verbose) cat("Replication: ",i,"\t Fold: ",ii," \n")

	   # CV in R with comitting variance components
	   if (VC.est=="commit"){
		# select ES, k-times
		samp.es <- val.samp3[val.samp3[,2]!=ii,] 
	
		# vectors and matrices for MME
		y1 <- y[y[,1] %in% samp.es[,1],2]
		Z1 <-Z[rownames(Z) %in% samp.es[,1],]
		X1 <-X[rownames(X) %in% samp.es[,1],]
		# crossproducts
		XX <- crossprod(X1)
		XZ <- crossprod(X1,Z1)
		ZX <- crossprod(Z1,X1)
		ZZGI <-  crossprod(Z1)+ GI
		Xy <- crossprod(X1,y1)
		Zy <- crossprod(Z1,y1)
		# Left hand side	
		LHS <- rbind(cbind(XX, XZ),cbind(ZX,ZZGI))
		# Right hand side
		RHS <- rbind(Xy,Zy)
		
		# solve MME
		bu <-  as.vector(ginv(LHS)%*%RHS)
		}

	   # estimation of variance components with ASReml for every ES
   	   if (VC.est=="ASReml"){
		samp.kf<-val.samp3[,2]==ii
		y.samp<-y
		y.samp[samp.kf,2]<-NA # set values of TS to NA
		#print(head(y.samp))

		# for unix
		if(.Platform$OS.type == "unix"){

			# checking directories for ASReml
			ASTest <- system(paste("ls"),intern=TRUE)
			if (!any(ASTest %in% "ASReml")) system(paste("mkdir ASReml"))

			# data output for ASReml
			write.table(y.samp,'Pheno.txt',col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')

			# ASReml function
			asreml <- system(paste('asreml -ns10000 Model.as',sep=''),TRUE)
			system(paste('asreml -p Model.pin',sep='')) # for variance components in an extra file

			system(paste('mv Model.asr ','ASReml/Model_rep',i,'_fold',ii,'.asr',sep=''))
			system(paste('mv Model.sln ','ASReml/Model_rep',i,'_fold',ii,'.sln',sep=''))
			system(paste('mv Model.vvp ','ASReml/Model_rep',i,'_fold',ii,'.vvp',sep=''))
			system(paste('mv Model.yht ','ASReml/Model_rep',i,'_fold',ii,'.vht',sep=''))
			system(paste('mv Model.pvc ','ASReml/Model_rep',i,'_fold',ii,'.pvc',sep=''))				
		}

		# for windows
		#if(.Platform$OS.type == "windows"){

			# checking directories for ASReml
		#	ASTest <- shell(paste("dir /b"),intern=TRUE)
		#	if (!any(ASTest %in% "ASReml")) shell(paste("md ASReml"))
			# data output for ASReml
		#	write.table(y.samp,'Pheno.txt',col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')
			# ASReml function
		#	system(paste('ASReml.exe -ns10000 Model.as',sep=''),wait=TRUE,show.output.on.console=FALSE)
		##	system(paste('ASReml.exe -p Model.pin',sep=''),wait=TRUE,show.output.on.console=FALSE)
		#	shell(paste('move Model.asr ','ASReml/Model_rep',i,'_fold',ii,'.asr',sep=''),wait=TRUE,translate=TRUE)
		#	shell(paste('move Model.sln ','ASReml/Model_rep',i,'_fold',ii,'.sln',sep=''),wait=TRUE,translate=TRUE)
		#	shell(paste('move Model.vvp ','ASReml/Model_rep',i,'_fold',ii,'.vvp',sep=''),wait=TRUE,translate=TRUE)
		#	shell(paste('move Model.yht ','ASReml/Model_rep',i,'_fold',ii,'.vht',sep=''),wait=TRUE,translate=TRUE)
		#	shell(paste('move Model.pvc ','ASReml/Model_rep',i,'_fold',ii,'.pvc',sep=''),wait=TRUE,translate=TRUE)				
		#}

		samp.es <- val.samp3[val.samp3[,2]!=ii,]
		# read in ASReml solutions
		asreml.sln<-matrix(scan(paste('ASReml/Model_rep',i,'_fold',ii,'.sln',sep=''),what='character'),ncol=4,byrow=TRUE)
		# solve MME
		bu <-  as.numeric(asreml.sln[,3])
	   }

	   # estimation of variance components with Bayesian ridge regression for every ES
   	   if (VC.est=="BRR"){
		samp.kf<-val.samp3[,2]==ii
		y.samp<-y
		y.samp[samp.kf,2]<-NA # set values of TS to NA
		#print(head(y.samp))

		# checking directories for BRR
		if(.Platform$OS.type == "unix"){
			BRRTest <- system(paste("ls"),intern=TRUE)
			if (!any(BRRTest %in% "BRR")) system(paste("mkdir BRR"))
		}
		if(.Platform$OS.type == "windows"){
			BRRTest <- shell(paste("dir /b"),intern=TRUE)
			if (!any(BRRTest %in% "BRR")) shell(paste("md BRR"))
		}

		# BRR function
		if(is.null(cov.matrix)) capture.output(mod50k <- BLR(y=y.samp[,2],XR=Z,prior=priorBLR,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=paste("BRR/50k_rep",i,"_fold",ii,sep="")),file=paste("BRR/BRRout_rep",i,"_fold",ii,".txt",sep=""))
		if(!is.null(cov.matrix)) capture.output(mod50k <- BLR(y=y.samp[,2],GF=list(ID=1:n,A=cov.matrix[[1]]),prior=priorBLR,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=paste("BRR/50k_rep",i,"_fold",ii,sep="")),file=paste("BRR/BRRout_rep",i,"_fold",ii,".txt",sep=""))

		samp.es <- val.samp3[val.samp3[,2]!=ii,]

		# solution
		if(is.null(cov.matrix)) bu <-  as.numeric(c(mod50k$mu,mod50k$bR))
		if(!is.null(cov.matrix)) bu <-  as.numeric(c(mod50k$mu,mod50k$u))
	  }

	  # estimation of variance components with Baysian Lasso for every ES
   	  if (VC.est=="BL"){
		samp.kf<-val.samp3[,2]==ii
		y.samp<-y
		y.samp[samp.kf,2]<-NA # set values of TS to NA
		#print(head(y.samp))

		# checking directory for BL
		if(.Platform$OS.type == "unix"){
			BLTest <- system(paste("ls"),intern=TRUE)
			if (!any(BLTest %in% "BL")) system(paste("mkdir BL"))
		}
		if(.Platform$OS.type == "windows"){
			BLTest <- shell(paste("dir /b"),intern=TRUE)
			if (!any(BLTest %in% "BL")) shell(paste("md BL"))
		}

		# BL function
		if(is.null(cov.matrix)) capture.output(mod50k <- BLR(y=y.samp[,2],XL=Z,prior=priorBLR,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=paste("BL/50k_rep",i,"_fold",ii,sep="")),file=paste("BL/BLout_rep",i,"_fold",ii,".txt",sep=""))
		if(!is.null(cov.matrix)) capture.output(mod50k <- BLR(y=y.samp[,2],GF=list(ID=1:n,A=cov.matrix[[1]]),prior=priorBLR,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=paste("BL/50k_rep",i,"_fold",ii,sep="")),file=paste("BL/BLout_rep",i,"_fold",ii,".txt",sep=""))

		samp.es <- val.samp3[val.samp3[,2]!=ii,]

		# solutions of BL

		if(is.null(cov.matrix)) bu <-  as.numeric(c(mod50k$mu,mod50k$bL))
		if(!is.null(cov.matrix)) bu <-  as.numeric(c(mod50k$mu,mod50k$u))
		#print(length(bu))
	  }

	  # solution vector
	  bu2 <- cbind(bu2,bu)
	  colnames(bu2)[ii]<-paste("rep",i,"_fold",ii,sep="")
	  # TS
	  Z2 <- Z[!(rownames(Z) %in% samp.es[,1]),]
	  X2 <- X[!(rownames(X) %in% samp.es[,1]),]
	  XZ2 <- cbind(X2,Z2)
	  #print(dim(XZ2))
	  y2 <- y[!(y[,1] %in% samp.es[,1]),2]
	  y.dach <- XZ2%*%bu
	  #print(head(y.dach))
	  #print(dim(y.dach))
	  n.TS <- rbind(n.TS,nrow(y.dach))
	  rownames(n.TS)[ii]<-paste("fold",ii,sep="")
          # Predicted breeding/testcross values
          y.TS <- rbind(y.TS,y.dach)
	  # predictive ability
	  COR <- round(cor(y2,y.dach),digits=4)
	  COR2 <- rbind(COR2,COR)
	  rownames(COR2)[ii]<-paste("fold",ii,sep="")
	  # regression = bias
	  lm1 <- lm(y2~as.numeric(y.dach))
	  #print(lm1)
	  #print(y.dach)
 	  lm.coeff <- rbind(lm.coeff,lm1$coefficients[2])
	  rownames(lm.coeff)[ii]<-paste("fold",ii,sep="")
	  # save IDs of TS
	  id.TS[[ii]] <- rownames(Z2)
	  names(id.TS)[[ii]] <- paste("fold",ii,sep="")
   	}  # end loop for k-folds

	n.TS2<-cbind(n.TS2,n.TS)
    	colnames(n.TS2)[i] <- paste("rep",i,sep="")
	y.TS <- y.TS[sort.list(rownames(y.TS)),]
    	y.TS2 <- cbind(y.TS2,y.TS)
	#print(dim(y.TS2))
    	colnames(y.TS2)[i] <- paste("rep",i,sep="")
    	COR3 <- cbind(COR3,COR2)
    	colnames(COR3)[i] <- paste("rep",i,sep="")
    	bu3 <- cbind(bu3,bu2)
	rownames(bu3)<-names.eff
    	lm.coeff2 <- cbind(lm.coeff2,lm.coeff)
    	colnames(lm.coeff2)[i] <- paste("rep",i,sep="")
	# save IDs of TS
	id.TS2[[i]] <- id.TS
	names(id.TS2)[[i]] <- paste("rep",i,sep="")
    }  # end loop for replication

    # return object
    if(VC.est=="commit") est.method <- "committed" else est.method <- paste("reestimated with ",VC.est,sep="")
    obj <- list( n.TS=n.TS2,id.TS=id.TS2,bu=bu3,y.TS=y.TS2,PredAbi=COR3,bias=lm.coeff2,k=k, Rep=Rep, sampling=sampling,Seed=Seed, rep.seed=seed2,nr.ranEff = length(cov.matrix),VC.est.method=est.method)
    class(obj) <- "cvData"
    return(obj)
}

