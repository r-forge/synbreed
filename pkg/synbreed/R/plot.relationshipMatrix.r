# heatmap for relationshipMatrix objects

plot.relationshipMatrix <- function(relationshipMatrix,...){

         class(relationshipMatrix) <- "matrix"
         n <- nrow(relationshipMatrix)
         if ( n < 35) levelplot(t(relationshipMatrix),axes=FALSE,col.regions=brewer.pal(9,"OrRd"),cuts=8,xlab="",ylab="",scales=list(cex=1-(n-20)/50,rot=c(40,0),abbreviate=TRUE,minlength=5),colorkey=list(labels=list(cex=1-(n-20)/50)),...)
         else levelplot(t(relationshipMatrix),axes=FALSE,col.regions=brewer.pal(9,"OrRd"),cuts=8,xlab="",ylab="",scales=list(at=c(1,n/2,n),labels=c(1,"...",paste(n)),tck=0),...)                                     
}