# summary method for cross validation results
summary.cvData <- function(object,...){

     obj <- object
     ans <- list()

     # summary of method
     ans$k <- obj$k
     ans$Rep <- obj$Rep
     ans$est.method <- obj$VC.est.method
     ans$sampling <- obj$sampling
     ans$nr.ranEff <- obj$nr.ranEff

     # number of individuals
     ans$n <- nrow(obj$y.TS)

     # size of TS
     ans$nmin.TS <- min(obj$n.TS)
     ans$nmax.TS <- max(obj$n.TS)

     # Results
     # Predictive ability
     colmean.pa <- colMeans(obj$PredAbi)
     ans$se.pa <- format(round(sd(colmean.pa)/sqrt(length(colmean.pa)),digits=4),digits=4,nsmall=4)
     ans$min.pa <- format(min(obj$PredAbi),digits=4,nsmall=4)
     ans$mean.pa <- format(mean(obj$PredAbi),digits=4,nsmall=4)
     ans$max.pa <- format(max(obj$PredAbi),digits=4,nsmall=4)
     # Bias
     colmean.b <- colMeans(obj$bias)
     ans$se.b <- format(round(sd(colmean.b)/sqrt(length(colmean.b)),digits=4),digits=4,nsmall=4)
     ans$min.b<- format(min(obj$bias),digits=4,nsmall=4)
     ans$mean.b <- format(mean(obj$bias),digits=4,nsmall=4)
     ans$max.b <-format(max(obj$bias),nsmall=4,digits=4)

     # Seed
     ans$Seed <- obj$Seed
     ans$rep.seed <- obj$rep.seed

     class(ans) <- "summary.cvData"
     ans
}

# print method for summary.pedigree
print.summary.cvData <- function(x,...){
    cat("Object of class 'cvData' \n")
    cat("\n",x$k,"-fold cross validation with",x$Rep,"replications \n")
    cat("     Sampling:                ",x$sampling,"\n")
    cat("     Variance components:     ",x$est.method,"\n")
    cat("     Number of random effects:",x$nr.ranEff,"\n")
    cat("     Number of individuals:   ",x$n," \n")
    cat("     Size of the TS:          ",x$nmin.TS,"--", x$nmax.TS," \n")
    cat("\nResults: \n")
    cat("                      Min         Mean +- pooled SE     Max \n")
    cat(" Predictive ability: ",x$min.pa,"    ",x$mean.pa,"+-",x$se.pa,"    ",x$max.pa,"\n")
    cat(" Bias:               ",x$min.b,"    ",x$mean.b,"+-",x$se.b,"    ",x$max.b,"\n")

    cat("\nSeed start: ",x$Seed,"\n")
    cat("Seed replications: \n")
    print(x$rep.seed)
} 



