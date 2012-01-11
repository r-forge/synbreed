# wrapper to apply genomic prediction models to an object of class gpData
# models used: regress and BLR

# author: Valentin Wimmer
# date: 2011 - 05 - 03
# changes: Hans-Jürgen Auinger
# date: 2011 - 11 - 21
# changes: return argument by Valentin Wimmer
# date: 2011 - 11 - 30

gpMod <- function(gpData,model=c("BLUP","BL","BRR"),kin=NULL,trait=1,repl=NULL,markerEffects=FALSE,fixed=NULL,random=NULL,...){
  ans <- list()
  model <- match.arg(model)
  m <- NULL
  if(length(trait) > 1) warning("\n   The return will be a list of gpMod-objects!\n")
  if(is.null(fixed)) fixed <- " ~ 1"
  if(is.null(random)) random <- "~ " else random <- paste(paste(random, collapse=" "), "+ ")
  if(model == "BLUP"){
    if (is.null(kin)){
      if(!gpData$info$codeGeno) stop("Missing object 'kin', or use function codeGeno first!")
      kin <- kin(gpData, ret="realized")
      # transposed crossproduct of the genotype matrix is used as relationship to obtain the variance components and mean of RR-BLUP 
      if(markerEffects) kin <- tcrossprod(gpData$geno) 
    } else { 
      if(markerEffects) kin <- tcrossprod(gpData$geno)
    }
  }
  for(i in trait){
    df.trait <- gpData2data.frame(gpData, i, onlyPheno=TRUE, repl=repl)
    # take data from gpData object
    vec.bool <- colnames(df.trait) == "ID" | colnames(df.trait) %in% unlist(strsplit(paste(fixed), " ")) | colnames(df.trait) %in% unlist(strsplit(paste(random), " "))
    if(i %in% 1:(dim(gpData$pheno)[2]+ dim(gpData$phenoCovars)[2])) {
      yName <- dimnames(gpData$pheno)[[2]][as.numeric(i)]
      vec.bool[colnames(df.trait) %in% yName] <- TRUE 
    } else {
      vec.bool <- vec.bool | colnames(df.trait) == i
      yName <- i
    }
    df.trait <- df.trait[, vec.bool]
    df.trait <- df.trait[!apply(is.na(df.trait), 1, sum), ]
    kinNames <- unique(df.trait$ID[!df.trait$ID %in% rownames(kin)])
    if(length(kinNames) != 0 & model=="BLUP"){
      df.trait <- df.trait[!df.trait$ID %in% kinNames,]
      warning("Some phenotyped IDs are not in the kinship matrix!\nThese are removed from the analysis")
    }
    kinTS <- kin[df.trait$ID, df.trait$ID]# expand the matrix to what is needed
    if(model == "BLUP"){
      res <- regress(as.formula(paste(yName, paste(fixed, collapse=" "))), Vformula=as.formula(paste(paste(random, collapse=" "), "kinTS")),data=df.trait, identity=TRUE,...)
      us <- BLUP(res)$Mean
      genVal <- us[grep("kinTS", names(us))]
      genVal <- genVal[!duplicated(names(genVal))]
      names(genVal) <-  unlist(strsplit(names(genVal), "kinTS."))[(1:length(genVal))*2]
      # genVal <- NULL
      if(markerEffects){
        sigma2m <- res$sigma["kinTS"]
        # set up design matrices
        Z <- NULL; X <- NULL
        X <- model.matrix(as.formula(fixed), data=df.trait)# fixed part of the model
        if(substr(random, nchar(random)-1, nchar(random)-1) == "+"){
          randomM <- substr(random, 1, nchar(random)-3)
          term <- labels(terms(as.formula(randomM)))
          if(!all(term %in% colnames(df.trait))) stop("for markerEffects = TRUE only factors or regressors as random covariables are allowed!")
          Z <-  res$Z[[term[1]]]
          colnames(Z) <- paste(term[1], colnames(Z), sep=".")
          if(length(term) > 1) for(ii in 2:length(term)){
            Z1 <- res$Z[[term[ii]]]
            colnames(Z1) <- paste(term[ii], colnames(Z1), sep=".")
            Z <- cbind(Z, Z1)  
          }
          Z <- cbind(Z, gpData$geno[df.trait$ID, ]) 
          GI <- diag(ncol(Z))
          cnt <- 1
          for(ii in term){
            cnt1 <- nlevels(df.trait[, ii])-1
            if(cnt1 == 0 ) cnt1 <- 0
            GI[cnt:(cnt+cnt1), cnt:(cnt+cnt1)] <- GI[cnt:(cnt+cnt1), cnt:(cnt+cnt1)] / res$sigma[ii] * res$sigma["In"]
            cnt <- cnt + cnt1 + 1
          }
        } else {
          Z <- gpData$geno[df.trait$ID, ]
          GI <- diag(ncol(Z))
        }
        GI[(nrow(GI)-ncol(gpData$geno)+1):nrow(GI), (nrow(GI)-ncol(gpData$geno)+1):nrow(GI)] <- 
                GI[(nrow(GI)-ncol(gpData$geno)+1):nrow(GI), (nrow(GI)-ncol(gpData$geno)+1):nrow(GI)] / sigma2m * res$sigma["In"]
        RI <- res$Z$In
        sol <- MME(X, Z, GI, RI, df.trait[, yName])
        m <- sol$u 
        names(m) <- colnames(Z)
        m <- m[colnames(gpData$geno)]
      }
    }

    if(model=="BL"){
      if(dim(gpData$pheno)[3] > 1) stop("This method is not developed for a one-stage analysis yet. \nA phenotypic analysis have to be done fist.")
       X <- gpData$geno[rownames(gpData$geno) %in% df.trait$ID,]
       y <- df.trait[df.trait$ID %in% rownames(gpData$geno), yName]
      capture.output(res <- BLR(y=y,XL=X,...),file="BLRout.txt")
      if(!is.null(kin)) capture.output(res <- BLR(y=y,XL=X,GF=list(ID=df.trait$ID,A=kinTS),...),file="BLRout.txt")
      else kin <- NULL
      genVal <- res$yHat
      names(genVal) <- rownames(X)
      m <- res$bL
    }
    if(model=="BRR"){
      if(dim(gpData$pheno)[3] > 1) stop("This method is not developed for a one-stage analysis yet. \nA phenotypic analysis have to be done fist.")
        X <-  gpData$geno[rownames(gpData$geno) %in% df.trait$ID,]
        y <- df.trait[df.trait$ID %in% rownames(gpData$geno), yName]
        capture.output(res <- BLR(y=y,XR=X,...),file="BLRout.txt")
        if(!is.null(kin))  capture.output(res <- BLR(y=y,XR=X,GF=list(ID=df.trait$ID,A=kin),...),file="BLRout.txt")
        else kin <- NULL
        genVal <- res$yHat
        names(genVal) <- rownames(X)
        m <- res$bR
      }

    y <- df.trait[,yName]
    names(y) <- df.trait[,"ID"]
    ret <- list(fit=res,model=model,y=y,g=genVal,m=m,kin=kin)
    class(ret) = "gpMod"
    ans[[i]] <- ret
    names(ans)[length(ans)] <- yName
  }
  if(length(trait) > 1){
   class(ans) <- "gpModList"
   return(ans)
  } else {
    return(ret)
  }
}

summary.gpMod <- function(object,...){
    ans <- list()
    ans$model <- object$model
    if(object$model %in% c("BLUP")) ans$summaryFit <- summary(object$fit)
    if(object$model=="BL") ans$summaryFit <- list(mu = object$fit$mu, varE=object$fit$varE, varU=object$fit$varU, lambda=object$fit$lambda, nIter = object$fit$nIter,burnIn = object$fit$burnIn,thin=object$fit$thin)
    if(object$model=="BRR") ans$summaryFit <- list(mu = object$fit$mu, varE=object$fit$varE, varBr=object$fit$varBr, varU=object$fit$varU, nIter = object$fit$nIter,burnIn = object$fit$burnIn,thin=object$fit$thin)
    ans$n <- sum(!is.na(object$y))
    ans$sumNA <- sum(is.na(object$y))
    ans$summaryG <- summary(as.numeric(object$g))
    class(ans) <- "summary.gpMod"
    ans
}

summary.gpModList <- function(object,...){
  ret <- list()
  for(i in 1:length(object)){
    ret[[i]] <- summary(object[[i]])
  }
  class(ret) <- "summary.gpModList"
  names(ret) <- names(object)
  ret
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
    if(x$model %in% c("BLUP")) 
      cat(print(x$summaryFit),"\n")
    else {
      cat("MCMC options: nIter = ",x$summaryFit$nIter,", burnIn = ",x$summaryFit$burnIn,", thin = ",x$summaryFit$thin,"\n",sep="")
      cat("             Posterior mean \n")
      cat("(Intercept) ",x$summaryFit$mu,"\n")
      cat("VarE        ",x$summaryFit$varE,"\n")
      if(!is.null(x$summaryFit$varU)) cat("VarU        ",x$summaryFit$varU,"\n")
      if(x$model=="BL"){
        cat("lambda      ",x$summaryFit$lambda,"\n")
      }
      if(x$model=="BRR"){
        cat("varBr       ",x$summaryFit$varBr,"\n")
      }
    }
}

print.summary.gpModList <- function(x,...) {
  for(i in 1:length(x)){
    cat(paste("\n\tTrait ", names(x)[i], "\n\n\n"))
    print(x[[i]]) 
    if(i != length(x)) cat("-------------------------------------------------\n") 
  }
}
