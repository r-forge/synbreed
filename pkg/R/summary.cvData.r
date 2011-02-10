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
     ans$se.pa <- round(sd(colmean.pa)/sqrt(length(colmean.pa)),digits=4)
     ans$min.pa <- round(min(obj$PredAbi),digits=4)
     ans$mean.pa <- round(mean(obj$PredAbi),digits=4)
     ans$max.pa <- round(max(obj$PredAbi),digits=4)
     # Bias
     colmean.b <- colMeans(obj$bias)
     ans$se.b <- round(sd(colmean.b)/sqrt(length(colmean.b)),digits=4)
     ans$min.b<- round(min(obj$bias),digits=4)
     ans$mean.b <- round(mean(obj$bias),digits=4)
     ans$max.b <- round(max(obj$bias),digits=4)

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
    cat("\tSampling:\t \t",x$sampling,"\n")
    cat("\tVariance components: \t",x$est.method,"\n")
    cat("\tNumber of random effects:",x$nr.ranEff,"\n")
    cat("\tNumber of individuals:\t",x$n," \n")
    cat("\tSize of the TS:\t \t",x$nmin.TS,"--", x$nmax.TS," \n")
    cat("\nResults: \n")
    cat(" \t \t      Min \t Mean +- pooled SE \t Max \n")
    cat(" Predictive ability: ",x$min.pa,"\t",x$mean.pa,"+-",x$se.pa,"\t",x$max.pa,"\n")
    cat(" Bias: \t \t     ",x$min.b,"\t",x$mean.b,"+-",x$se.b,"\t",x$max.b,"\n")

    cat("\nSeed start: ",x$Seed,"\n")
    cat("Seed replications: \n")
    print(x$rep.seed)
} 



