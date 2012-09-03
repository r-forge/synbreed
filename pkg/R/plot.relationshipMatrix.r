# heatmap for relationshipMatrix objects

plot.relationshipMatrix <- function(x, file=NULL, fileFormat="pdf",...){
    # open device in file
    if(!is.null(file)){
      if(substr(file, nchar(file)-nchar(fileFormat)+1, nchar(file)) != fileFormat | nchar(file) < 5)
        file <- paste(file, ".", fileFormat, sep="")
      if(fileFormat == "pdf") pdf(file)
      else if (fileFormat == "png") png(file)
      else stop("not supported file format choosen!")
    }


    relationshipMatrix <- x
    class(relationshipMatrix) <- "matrix"
    n <- nrow(relationshipMatrix)
    color = c("#FFF7EC","#FEE8C8","#FDD49E","#FDBB84","#FC8D59","#EF6548","#D7301F","#B30000","#7F0000")
 
    if ( n < 35) levelplot(t(relationshipMatrix),axes=FALSE,col.regions=color,cuts=8,xlab="",ylab="",scales=list(cex=1-(n-20)/50,rot=c(40,0),abbreviate=TRUE,minlength=5),ylim=c(n+1,0),...)
    else levelplot(t(relationshipMatrix),axes=FALSE,col.regions=color,cuts=8,xlab="",ylab="",scales=list(at=c(1,n/2,n),labels=c(1,"...",paste(n)),tck=0),ylim=c(n+1,0),...)  
    if(!is.null(file)) dev.off()
                                   
}
