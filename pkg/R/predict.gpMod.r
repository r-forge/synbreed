# predictions for objects of class gpMod

predict.gpMod <- function(object,newdata=NULL,...){
  if(class(object)!="gpMod") stop("'object' must be of class 'gpMod'")
  if (is.null(newdata)) prediction <- object$g
  else{
  model <- object$model
  
  if(model %in% c("modBRR","modBL")){
      if(!is.null(object$kin)) ("including a polygenic effect is not yet implemented")
      if(class(newdata)!="gpData") stop("object 'newdata' must be of class 'gpData'")
      X <- newdata$geno
      m <- object$m
      mu <- object$fit$mu
      prediction <- mu + X %*% m
  }
  if(model=="modRR"){
     if(class(newdata)!="gpData") stop("object 'newdata' must be of class 'gpData'")
     X <- newdata$geno
     m <- object$m
     mu <- object$fit$beta
     prediction <- mu + X %*% m
  }
  if(model %in% c("modA","modU")){
      G <- object$kin[c(object$trainingSet,newdata),c(object$trainingSet,newdata)]
      y <- object$y
      n <- length(y) # size of the training set
      Z <- cbind(diag(n),matrix(data=0,nrow=n,ncol=length(newdata)))
      X <- matrix(1,nrow=n)
      sigma2g <- object$fit$sigma[1]
      sigma2  <- object$fit$sigma[2]
      diag(G) <- diag(G) + 0.00001
      GI <- solve(G)*sigma2/sigma2g    # to ensure a solution for the inverse
      RI <- diag(n)
      sol <- MME(X, Z, GI, RI, y)
      prediction <- sol$b + sol$u[-(1:n)]
      names(prediction) <- newdata
  }

  }
  return(prediction)
}                                         
                             