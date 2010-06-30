# heatmap for relationshipMatrix objects

plot.relationshipMatrix <- function(relationshipMatrix,...){

         class(relationshipMatrix) <- "matrix"
         levelplot(relationshipMatrix,axes=FALSE,col.regions=brewer.pal(9,"OrRd"),cuts=8,...)

}