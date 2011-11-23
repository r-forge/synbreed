# wrapper to apply genomic prediction models to an object of class gpData
# models used: regress and BLR

# author: Valentin Wimmer
# date: 2011 - 05 - 03
# changes: Hans-J�rgen Auinger
# date: 2011 - 11 - 21

gpMod <- function(gpData,model=c("BLUP","modBL","modBRR"),kin=NULL,trait=1,repl=NULL,markerEffects=FALSE,fixed=NULL,random=NULL,...){
  model <- match.arg(model)
  m <- NULL
  for(i in trait){
    df.trait <- gpData2data.frame(gpData, i, onlyPheno=TRUE, repl=repl)
    # take data from gpData object
    if (is.null(kin)){
      if(!gpData$info$codeGeno) stop("Missing object 'kin'")
      kin <- kin(gpData, ret="realized")
    }
    vec.bool <- colnames(df.trait) == "ID" | colnames(df.trait) %in% unlist(strsplit(paste(fixed), " ")) | colnames(df.trait) %in% unlist(strsplit(paste(random), " "))
    if(i %in% 1:ncol(df.trait)) {
      yName <- dimnames(gpData$pheno)[[2]][as.numeric(i)]
      vec.bool[colnames(df.trait) %in% yName] <- TRUE 
    } else {
      vec.bool <- vec.bool | colnames(df.trait) == i
      yName <- i
    }
    df.trait <- df.trait[, vec.bool]
    df.trait <- df.trait[!apply(is.na(df.trait), 1, sum), ]
    kinNames <- unique(df.trait$ID[!df.trait$ID %in% rownames(kin)])
    if(length(kinNames) != 0){
      df.trait <- df.trait[!df.trait$ID %in% kinNames,]
      warning("Some phenotyped IDs are not in the kinship matrix!\nThese are removed from the analysis")
    }
    kinTS <- kin[df.trait$ID, df.trait$ID]

    if(model == "BLUP"){
      if(is.null(fixed)) fixed <- " ~ 1"
      if(is.null(random)) random <- "~ " else random <- paste(paste(random, collapse=" "), " + ")
      res <- regress(as.formula(paste(yName, paste(fixed, collapse=" "))), Vformula=as.formula(paste(paste(random, collapse=" "), "kinTS")),data=df.trait,...)
      genVal <- res$predicted
      if(markerEffects){
        genVal <- res$predicted
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

    if(model=="modBL"){
      if(dim(gpData$pheno)[3] > 1) stop("This method is not developed for a one-stage analysis yet. \nA phenotypic analysis have to be done fist.")
      X <- gpData$geno[df.trait$ID, ]
      capture.output(res <- BLR(y=df.trait[, yName],XL=X,...),file="BLRout.txt")
      if(!is.null(kin)) res <- BLR(y=df.trait[, yName],XL=X,GF=list(ID=1:n,A=kinTS),...)
      genVal <- res$yHat
      names(genVal) <- rownames(X)
      m <- res$bL
    }
    if(model=="modBRR"){
      if(dim(gpData$pheno)[3] > 1) stop("This method is not developed for a one-stage analysis yet. \nA phenotypic analysis have to be done fist.")
      X <- gpData$geno[rownames(gpData$geno) %in% trainSet,]
      capture.output(res <- BLR(y=df.trait[, yName],XR=X,...),file="BLRout.txt")
      if(!is.null(kin)) res <- BLR(y=df.trait[, yName],XR=X,GF=list(ID=1:n,A=kin),...)
      genVal <- res$yHat
      names(genVal) <- rownames(X)
      m <- res$bR
    }

    
    ret <- list(fit=res,model=model,trainingSet=trainSet,y=y,g=genVal,m=m,kin=kin)
    class(ret) = "gpMod"
    return(ret)
  }
}

summary.gpMod <- function(object,...){
    ans <- list()
    ans$model <- object$model
    if(object$model %in% c("BLUP")) ans$summaryFit <- summary(object$fit)
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
    if(x$model %in% c("BLUP")) cat(print(x$summaryFit),"\n")
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
