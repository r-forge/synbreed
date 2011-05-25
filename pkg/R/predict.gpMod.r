# predictions for objects of class gpMod

predict.gpMod <- function(object,newdata,...){
  if(class(object)!="gpMod") stop("'object' must be of class 'gpMod'")
  model <- object$model
  if(model %in% c("BRR","BL")){
      if(!is.null(object$kin)) ("including a polygenic effect is not yet implemented")
      if(class(newdata)!="gpData") stop("object 'newdata' must be of class 'gpData'")
      X <- newdata$geno
      m <- object$m
      mu <- object$fit$mu
      prediction <- mu + X %*% m
  }
  if(model=="RR BLUP"){
     if(class(newdata)!="gpData") stop("object 'newdata' must be of class 'gpData'")
     X <- newdata$geno
     m <- object$m
     mu <- object$fit$beta
     prediction <- mu + X %*% m
  }
  if(model=="BLUP"){
      warning("valid results only for a pedigree-based relationship matrix, prediction for a realized relationship matrix is not yet implemented")
      G <- object$kin[c(object$trainingSet,newdata),c(object$trainingSet,newdata)]
      y <- object$y
      n <- length(y) # size of the training set
      Z <- cbind(diag(n),matrix(data=0,nrow=n,ncol=length(newdata)))
      X <- matrix(1,nrow=n)
      sigma2g <- object$fit$sigma[1]
      sigma2  <- object$fit$sigma[2]
      GI <- ginv(G)*sigma2/sigma2g
      RI <- diag(n)
      sol <- MME(X, Z, GI, RI, y)
      prediction <- sol$b + sol$u[-(1:n)]
      names(prediction) <- newdata
  }
  return(prediction)
}                                         
                             