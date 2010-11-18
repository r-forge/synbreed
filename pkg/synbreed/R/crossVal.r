# Cross validation with different sampling and variance components estimation methods

crossVal <- function (y,X,Z,cov.matrix=NULL, k=2,Rep=1,Seed=NULL,sampling=c("random","within family","across family"), varComp=NULL,popStruc=NULL, VC.est=FALSE) 
{
    # number of observations 
    n <- length(unique(y[,1]))
    names.eff <- c(colnames(X),colnames(Z))

    # catch errors	
    if(is.null(varComp) & !VC.est) stop("Variance components have to be specified")
    if(!VC.est & length(varComp)<2) stop("Variance components should be at least two, one for the random effect and one residual variance")
    if(sampling!="random" & is.null(popStruc)) stop("popStruc has to be given")
    if(sampling!="random" & !is.null(popStruc)){
      if(length(popStruc)!=n) stop("population structure must have equal length as obsersvations in data")
      if(any(is.na(popStruc))) stop("no missing values allowed in popStruc")
    }
    if ( k < 2) stop("folds should be equal or greater than 2")
    if ( k > n) stop("folds should be equal or less than the number of observations")
    if (is.null(cov.matrix)){
	warning("no covariance matrix is given, assuming one iid random effect")
	cov.matrix <- list(diag(ncol(Z)))
    }
    if (!VC.est & length(cov.matrix)!=(length(varComp)-1)) stop("number of variance components does not match given covariance matrices")

    # prepare X,Z design matrices
    X <- as.matrix(X)
    rownames(X) <- y[,1]
    Z <- as.matrix(Z)
    rownames(Z) <- y[,1]

    # prepare covariance matrices
    m <- length(cov.matrix)
    cat("Model with ",m," covariance matrix/ces \n")
    if (VC.est==FALSE){
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
      	for ( i in 1:length(cov.matrix)){
	relMatASReml(as.matrix(cov.matrix[[i]]),ginv=FALSE,file=paste("ID",i,".giv",sep=""),digits=10)
	}
	cat(paste("Model \n"," ID      	  !A \n"," Yield  	  !D* \n","ID1.giv \n","Pheno.txt !skip 1 !AISING !maxit 11\n","!MVINCLUDE \n \n","Yield  ~ mu !r giv(ID,1)",sep=""),file="cov1.as")
	cat(paste("Model \n"," ID      	  !A \n"," Yield  	  !D* \n","ID1.giv \n","ID2.giv \n","Pheno.txt !skip 1 !AISING !maxit 11\n","!MVINCLUDE \n \n","Yield  ~ mu !r giv(ID,1) giv(ID,2)",sep=""),file="cov2.as")
	cat(paste("Model \n"," ID      	  !A \n"," Yield  	  !D* \n","ID1.giv \n","ID2.giv \n","ID3.giv \n","Pheno.txt !skip 1 !AISING !maxit 11\n","!MVINCLUDE \n \n","Yield  ~ mu !r giv(ID,1) giv(ID,2) giv(ID,3)",sep=""),file="cov3.as")
	cat("",file="cov1.pin")
	cat("",file="cov2.pin")
	cat("",file="cov3.pin")
     }

    # set seed for replications
    if(!is.null(Seed)) set.seed(Seed)
    seed2<-round(runif(Rep,1,100000),0)

    # begin replications
    COR3 <- NULL
    bu3 <- NULL
    lm.coeff2 <- NULL
    y.TS2 <- NULL
    for (i in 1:Rep){ 

	# sampling of k sets
	# random sampling
	cat(sampling," sampling \n")
	if(sampling=="random"){
	  y.u <- unique(y[,1])
	  set.seed(seed2[i])
	  modu<-n%%k
	  val.samp2<-sample(c(rep(1:k,each=(n-modu)/k),sample(1:k,modu)),n,replace=FALSE)
	  val.samp3 <- as.data.frame(cbind(y.u,val.samp2))
	  }

	# within family sampling
	if(sampling=="within family"){
	   which.pop <- unique(popStruc)# nr of families
	   y.u <- unique(y[,1])
	   val.samp3<- NULL
	   for (j in 1:length(which.pop)){
		y2<-matrix(y.u[popStruc==which.pop[j]],ncol=1)# select each family
		set.seed(seed2[i])
		modu<-nrow(y2)%%k
		if(!modu==0) val.samp<-sample(c(rep(1:k,each=(nrow(y2)-modu)/k),sample(1:k,modu)),nrow(y2),replace=FALSE)
		if(modu==0) val.samp<-sample(rep(1:k,each=(nrow(y2))/k),nrow(y2),replace=FALSE)
		val.samp2 <- cbind(y2,val.samp)		
		val.samp3 <- as.data.frame(rbind(val.samp3,val.samp2))
	   }
	   val.samp3 <- val.samp3[order(as.character(val.samp3[,1])),]
	}

	# across family sampling
	if(sampling=="across family"){
	  which.pop <- unique(popStruc)
	  y.u <- unique(y[,1])
	  y2 <- matrix(y.u[order(popStruc),],ncol=1)
	  b <- table(popStruc)
	  modu<-length(which.pop)%%k
	  val.samp<-sample(c(rep(1:k,each=(length(which.pop)-modu)/k),sample(1:k,modu)),length(which.pop),replace=FALSE)
	  val.samp2<- rep(val.samp,b)
	  val.samp3 <- cbind(y2,val.samp2)
	  val.samp3 <- 	as.data.frame(val.samp3[order(as.character(val.samp3[,1])),])
	 }

   # CV in R with comitting variance components
   if (!VC.est){
     # start k folds
     COR2 <- NULL
     bu2 <- NULL
     lm.coeff <- NULL
     y.TS <- NULL
     for (ii in 1:k){
	cat("Replication: ",i,"\t Fold: ",ii," \n")
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
	bu2 <- cbind(bu2,bu)
	colnames(bu2)[ii]<-paste("rep",i,"_fold",ii,sep="")

	# TS
	Z2 <- Z[!(rownames(Z) %in% samp.es[,1]),]
	X2 <- X[!(rownames(X) %in% samp.es[,1]),]
	XZ2 <- cbind(X2,Z2)
	y2 <- y[!(y[,1] %in% samp.es[,1]),2]
	y.dach <- XZ2%*%bu
	# Predicted breeding/testcross values
	y.TS <- rbind(y.TS,y.dach)
	# predictive ability
	COR <- round(cor(y2,y.dach),digits=4)
	COR2 <- rbind(COR2,COR)
	rownames(COR2)[ii]<-paste("fold",ii,sep="")
	# regression = bias
	lm1 <- lm(y2~y.dach)
	lm.coeff <- rbind(lm.coeff,lm1$coefficients[2])
	rownames(lm.coeff)[ii]<-paste("fold",ii,sep="")
	}
    }
    if (VC.est){# estimation of variance components with ASReml for every ES
		# start k folds
     		COR2 <- NULL
     		bu2 <- NULL
     		lm.coeff <- NULL
     		y.TS <- NULL
		for (ii in 1:k){
			cat('\n number of cov-matrix:',m,'Replication: ',i,'\t Fold: ',ii,'\n \n')
			samp.kf<-val.samp3==ii
			y.samp<-y
			y.samp[y.samp[,1] %in% samp.kf[,2],2]<-NA # set values of TS to NA

			# for unix
			if(.Platform$OS.type == "unix"){
				write.table(y.samp,'Pheno.txt',col.names=TRUE,row.names=TRUE,quote=FALSE,sep='\t')
				asreml <- system(paste('asreml -ns10000 cov',m,'.as',sep=''),TRUE)
				system(paste('asreml -p cov',m,'.pin',sep=''))
				system(paste('mv cov',m,'.asr ','cov',m,'_rep',i,'_fold',ii,'.asr',sep=''))
				system(paste('mv cov',m,'.sln ','cov',m,'_rep',i,'_fold',ii,'.sln',sep=''))
				system(paste('mv cov',m,'.vvp ','cov',m,'_rep',i,'_fold',ii,'.vvp',sep=''))
				system(paste('mv cov',m,'.yht ','cov',m,'_rep',i,'_fold',ii,'.vht',sep=''))
				system(paste('mv cov',m,'.pvc ','cov',m,'_rep',i,'_fold',ii,'.pvc',sep=''))				
			}

			# for windows
			if(.Platform$OS.type == "windows"){
				write.table(y.samp,'Pheno.txt',col.names=TRUE,row.names=TRUE,quote=FALSE,sep='\t')
				system(paste('ASReml.exe -ns10000 cov',m,'.as',sep=''),wait=TRUE,show.output.on.console=FALSE)
				system(paste('ASReml.exe -p cov',m,'.pin',sep=''),wait=TRUE,show.output.on.console=FALSE)
				shell(paste('move cov',m,'.asr ','cov',m,'_rep',i,'_fold',ii,'.asr',sep=''),wait=TRUE,translate=TRUE)
				shell(paste('move cov',m,'.sln ','cov',m,'_rep',i,'_fold',ii,'.sln',sep=''),wait=TRUE,translate=TRUE)
				shell(paste('move cov',m,'.vvp ','cov',m,'_rep',i,'_fold',ii,'.vvp',sep=''),wait=TRUE,translate=TRUE)
				shell(paste('move cov',m,'.yht ','cov',m,'_rep',i,'_fold',ii,'.vht',sep=''),wait=TRUE,translate=TRUE)
				shell(paste('move cov',m,'.pvc ','cov',m,'_rep',i,'_fold',ii,'.pvc',sep=''),wait=TRUE,translate=TRUE)				
			}

			samp.es <- val.samp3[val.samp3[,2]!=ii,]
			# read in ASReml solutions
			asreml.sln<-matrix(scan(paste('cov',m,'_rep',i,'_fold',ii,'.sln',sep=''),what='character'),ncol=4,byrow=TRUE)
			# solve MME
			bu <-  as.numeric(asreml.sln[,3])
			bu2 <- cbind(bu2,bu)
			colnames(bu2)[ii]<-paste("rep",i,"_fold",ii,sep="")

			# TS
			Z2 <- Z[!(rownames(Z) %in% samp.es[,1]),]
			X2 <- X[!(rownames(X) %in% samp.es[,1]),]
			XZ2 <- cbind(X2,Z2)
			y2 <- y[!(y[,1] %in% samp.es[,1]),2]
			y.dach <- XZ2%*%bu
			# Predicted breeding/testcross values of TS
			y.TS <- rbind(y.TS,y.dach)
			# predictive ability
			COR <- round(cor(y2,y.dach),digits=4)
			COR2 <- rbind(COR2,COR)
			rownames(COR2)[ii]<-paste("fold",ii,sep="")
			# regression = bias
			lm1 <- lm(y2~y.dach)
			lm.coeff <- rbind(lm.coeff,lm1$coefficients[2])
			rownames(lm.coeff)[ii]<-paste("fold",ii,sep="")
		}
	}
	y.TS <- y.TS[order(rownames(y.TS)),]
    	y.TS2 <- cbind(y.TS2,y.TS)
    	colnames(y.TS2)[i] <- paste("rep",i,sep="")
    	COR3 <- cbind(COR3,COR2)
    	colnames(COR3)[i] <- paste("rep",i,sep="")
    	bu3 <- cbind(bu3,bu2)
	rownames(bu3)<-names.eff
    	lm.coeff2 <- cbind(lm.coeff2,lm.coeff)
    	colnames(lm.coeff2)[i] <- paste("rep",i,sep="")
    }
    # return object
    obj <- list(bu=bu3,y.TS=y.TS2,PredAbi=COR3,bias=lm.coeff2)
    class(obj) <- "cvData"
    return(obj)
}

