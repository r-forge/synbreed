# wrapper to apply genomic prediction models to an object of class gpData
# models used: regress and BLR

# author: Valentin Wimmer
# date: 2011 - 05 - 03

gpMod <- function(gpData,model=c("BLUP","BL","BRR","RR BLUP"),kin=NULL,trait=1,...){
    model <- match.arg(model)
    
    # inidividuls with genotypes and phenotypes
    trainSet <- as.character(gpData$covar$id[gpData$covar$genotyped & gpData$covar$phenotyped])
    n <- length(trainSet)

    # take data from gpData object
    y <- gpData$pheno[rownames(gpData$pheno) %in% trainSet,trait]

    if(model=="BLUP"){
      if (is.null(kin)) stop("Missing object 'kin'")
      if (n!=nrow(kin)){
         kin <- kin[rownames(kin) %in% trainSet,rownames(kin) %in% trainSet]
      }
      res <- regress(y~1,Vformula=~kin,...)
      genVal <- res$predicted
      m <- NULL
    }
    if(model=="BL"){
      X <- gpData$geno[rownames(gpData$geno) %in% trainSet,]
      capture.output(res <- BLR(y=y,XL=X,...),file="BLRout.txt")
      if(!is.null(kin)) res <- BLR(y=y,XL=X,GF=list(ID=1:n,A=kin),...)
      genVal <- res$yHat
      m <- res$bL
    }
    if(model=="BRR"){
      X <- gpData$geno[rownames(gpData$geno) %in% trainSet,]
      capture.output(res <- BLR(y=y,XR=X,...),file="BLRout.txt")
      if(!is.null(kin)) res <- BLR(y=y,XR=X,GF=list(ID=1:n,A=kin),...)
      genVal <- res$yHat
      m <- res$bR
    }
    if(model=="RR BLUP"){
      if(!gpData$info$codeGeno) stop("use function 'codeGeno' before")
      if (is.null(kin)) stop("Missing object 'kin'")
      if (n!=nrow(kin)){
         kin <- kin[rownames(kin) %in% trainSet,rownames(kin) %in% trainSet]
      }
      # first run BLUP model
      res <- regress(y~1,Vformula=~kin,...)
      genVal <- res$predicted
      sigma2u <- res$sigma[1]
      sigma2  <- res$sigma[2]
      p <- colMeans(gpData$geno)/2
      sumP <- 2*sum(p*(1-p))
      # use transformation rule for vc (Albrecht et al. 2011)
      sigma2m <- sigma2u/sumP
      # set up design matrices
      X <- matrix(1,nrow=n)
      Z <- gpData$geno
      GI <- diag(rep(sigma2/sigma2m,ncol(Z)))
      RI <- diag(n)
      sol <- MME(X, Z, GI, RI, y)
      m <- sol$u 
    }

    ret <- list(fit=res,model=model,y=y,g=genVal,m=m,kin=kin)
    class(ret) = "gpMod"
    return(ret)

}

summary.gpMod <- function(object,...){
    ans <- list()
    ans$model <- object$model
    if(object$model=="BLUP") ans$summaryFit <- summary(object$fit)
    if(object$model=="BL") ans$summaryFit <- list(mu = object$fit$mu, varE=object$fit$varE, lambda=object$fit$lambda, nIter = object$fit$nIter,burnIn = object$fit$burnIn,thin=object$fit$thin)
    if(object$model=="BRR") ans$summaryFit <- list(mu = object$fit$mu, varE=object$fit$varE, varBr=object$fit$varBr, nIter = object$fit$nIter,burnIn = object$fit$burnIn,thin=object$fit$thin)
    ans$n <- length(object$y)
    ans$sumNA <- sum(is.na(object$y))
    ans$summaryG <- summary(as.numeric(object$g))
    class(ans) <- "summary.gpMod"
    ans
}

print.summary.gpMod <- function(x,...){
    cat("object of class 'gpMod' \n")
    cat("Model used:",x$model,"\n")
    cat("nr. observations ",x$n," (NA = ",x$sumNA,") \n",sep="")
    cat("Genetic values: \n")
    cat("Min. 1st Qu. Median Mean 3rd Qu. Max \n")
    cat(x$summaryG, " \n")
    cat("--\n")
    cat("Model fit \n")
    if(x$model=="BLUP") cat(print(x$summaryFit),"\n")
    else{
    cat("MCMC options: nIter = ",x$summaryFit$nIter,", burnIn = ",x$summaryFit$burnIn,", thin = ",x$summaryFit$thin,"\n",sep="")
    cat("             Posterior mean \n")
    cat("(Intercept) ",x$summaryFit$mu,"\n")
    cat("VarE        ",x$summaryFit$varE,"\n")
    if(x$model=="BL"){
    cat("lambda      ",x$summaryFit$lambda,"\n")
    }
    if(x$model=="BRR"){
    cat("varBr       ",x$summaryFit$varBr,"\n")
    }
    }
}