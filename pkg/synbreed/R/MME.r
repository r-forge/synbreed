MME <- function(X,Z,GI,RI,y){
  SST <- t(y) %*% RI %*% y
  XX <- t(X) %*% RI %*% X
  XZ <- t(X) %*% RI %*% Z
  ZZ <- t(Z) %*% RI %*% Z
  Xy <- t(X) %*% RI %*% y
  Zy <- t(Z) %*% RI %*% y
  R1 <- cbind(XX,XZ)
  R2 <- cbind(t(XZ),(ZZ+GI))
  LHS <- rbind(R1,R2)
  RHS <- rbind(Xy,Zy)
  C <- ginv(LHS)
  bhat <- C %*% RHS
  SSR <- t(bhat) %*% RHS
  SSE <- SST - SSR
  N <- nrow(X)
  rX <- sum(diag(X%*%ginv(X)))
  sigma2e <- SSE/(N-rX)
  SEP <- sqrt(diag(C)*sigma2e)
  b <- bhat[1:ncol(X)]
  u <- bhat[-c(1:ncol(X))]
  return(list(b=b,u=u,LHS=LHS,RHS=RHS,C=C,SEP=SEP,SST=SST,SSR=SSR,residuals=y-X%*%b-Z%*%u))
  }
  
  set.seed(123)
n <- 10
X <- matrix(1,nrow=n,ncol=1)
Z <- cbind(matrix(0,10,10),diag(n))
ped <- simul.pedigree(3,c(5,5,10))
A <- kinship(ped,ret="add")
AI <- solve(A)
GI <- AI *  0.5
RI <- diag(n)
y <- c(10,12,13,15,8,14,7,8,9,9)
mod1 <- MME(X,Z,GI,RI,y)
mod1$u

library(regress)
mod2 <- regress(y~1,Vformula=~A[11:20,11:20])
summary(mod2)
mod2$predicted-fitted(mod2)


A <- A[11:20,11:20]
Z <- diag(n)
AI <- solve(A)
GI <- AI *  0.93/7.865
mod2 <- regress(y~1,Vformula=~A)
mod2$predicted-fitted(mod2)
mod1 <- MME(X,Z,GI,RI,y)
mod1$u