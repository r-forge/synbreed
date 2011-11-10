# wrapper to apply genomic prediction models to an object of class gpData
# models used: regress and BLR

# author: Valentin Wimmer
# date: 2011 - 05 - 03

gpMod <- function(gpData,model=c("modA","modU","modRR","modBL","modBRR"),kin=NULL,trait=1,...){
    model <- match.arg(model)
    if(model == "modA") {
      trainSet <- as.character(gpData$covar$id[gpData$covar$phenotyped ])
    } else {
      trainSet <- as.character(gpData$covar$id[gpData$covar$phenotyped & gpData$covar$genotyped ])
    }
      # remove missing observations
      trainSet <- trainSet[trainSet %in% rownames(gpData$pheno)[!is.na(gpData$pheno[,trait])]]

    if(model %in% c("modA","modU")){
      n <- length(trainSet)

      # take data from gpData object
      y <- gpData$pheno[rownames(gpData$pheno) %in% trainSet,trait]
    
      if (is.null(kin)) stop("Missing object 'kin'")
      if (n!=nrow(kin)){
         kinTS <- kin[rownames(kin) %in% trainSet,rownames(kin) %in% trainSet]
      } else kinTS <- kin

      res <- regress(y~1,Vformula=~kinTS,...)
      genVal <- res$predicted
      names(genVal) <- trainSet
      m <- NULL
    } else {
      n <- length(trainSet)

      # take data from gpData object
      y <- gpData$pheno[rownames(gpData$pheno) %in% trainSet,trait]

    if(model=="modBL"){
      X <- gpData$geno[rownames(gpData$geno) %in% trainSet,]
      capture.output(res <- BLR(y=y,XL=X,...),file="BLRout.txt")
      if(!is.null(kin)) res <- BLR(y=y,XL=X,GF=list(ID=1:n,A=kin),...)
      genVal <- res$yHat
      names(genVal) <- rownames(X)
      m <- res$bL
    }
    if(model=="modBRR"){
      X <- gpData$geno[rownames(gpData$geno) %in% trainSet,]
      capture.output(res <- BLR(y=y,XR=X,...),file="BLRout.txt")
      if(!is.null(kin)) res <- BLR(y=y,XR=X,GF=list(ID=1:n,A=kin),...)
      genVal <- res$yHat
      names(genVal) <- rownames(X)
      m <- res$bR
    }
    if(model=="modRR"){
      if(!gpData$info$codeGeno) stop("use function 'codeGeno' before")
      if (is.null(kin)) stop("Missing object 'kin'")
      if (n!=nrow(kin)){
         kin <- kin[rownames(kin) %in% trainSet,rownames(kin) %in% trainSet]
      }
      # first run BLUP model
      res <- regress(y~1,Vformula=~kin,...)
      genVal <- res$predicted
      names(genVal) <- trainSet
      sigma2u <- res$sigma[1]
      sigma2  <- res$sigma[2]
      p <- colMeans(gpData$geno)/2
      sumP <- 2*sum(p*(1-p))
      # use transformation rule for vc (Albrecht et al. 2011)
      sigma2m <- sigma2u/sumP
      # set up design matrices
      X <- matrix(1,nrow=n)
      Z <- gpData$geno[rownames(gpData$geno) %in% trainSet,]
      GI <- diag(rep(sigma2/sigma2m,ncol(Z)))
      RI <- diag(n)
      sol <- MME(X, Z, GI, RI, y)
      m <- sol$u 
      names(m) <- colnames(Z)
    }
    }
    
    ret <- list(fit=res,model=model,trainingSet=trainSet,y=y,g=genVal,m=m,kin=kin)
    class(ret) = "gpMod"
    return(ret)

}

summary.gpMod <- function(object,...){
    ans <- list()
    ans$model <- object$model
    if(object$model %in% c("modA","modU","modRR")) ans$summaryFit <- summary(object$fit)
    if(object$model=="modBL") ans$summaryFit <- list(mu = object$fit$mu, varE=object$fit$varE, lambda=object$fit$lambda, nIter = object$fit$nIter,burnIn = object$fit$burnIn,thin=object$fit$thin)
    if(object$model=="modBRR") ans$summaryFit <- list(mu = object$fit$mu, varE=object$fit$varE, varBr=object$fit$varBr, nIter = object$fit$nIter,burnIn = object$fit$burnIn,thin=object$fit$thin)
    ans$n <- sum(!is.na(object$y))
    ans$sumNA <- sum(is.na(object$y))
    ans$summaryG <- summary(as.numeric(object$g))
    class(ans) <- "summary.gpMod"
    ans
}

print.summary.gpMod <- function(x,...){
    cat("Object of class 'gpMod' \n")
    cat("Model used:",x$model,"\n")
    cat("Nr. observations ",x$n," \n",sep="")
    cat("Genetic performances: \n")
    cat("  Min.    1st Qu. Median  Mean    3rd Qu. Max    \n")
    cat(format(x$summaryG,width=7,trim=TRUE), "\n",sep=" ")
    cat("--\n")
    cat("Model fit \n")
    if(x$model %in% c("modA","modU","modRR")) cat(print(x$summaryFit),"\n")
    else{
    cat("MCMC options: nIter = ",x$summaryFit$nIter,", burnIn = ",x$summaryFit$burnIn,", thin = ",x$summaryFit$thin,"\n",sep="")
    cat("             Posterior mean \n")
    cat("(Intercept) ",x$summaryFit$mu,"\n")
    cat("VarE        ",x$summaryFit$varE,"\n")
    if(x$model=="modBL"){
    cat("lambda      ",x$summaryFit$lambda,"\n")
    }
    if(x$model=="modBRR"){
    cat("varBr       ",x$summaryFit$varBr,"\n")
    }
    }
}
