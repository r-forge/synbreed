manhattanPlot <- function(b,gpData=NULL,...){
    if(is.null(gpData)) plot(b,...)
    else{               
       cols <- rep(c(grey(0.3),grey(0.7)),times=length(unique(gpData$map$chr)))
       chr <- cumsum(as.numeric(gpData$map$pos))
       plot(chr,b,col=cols[as.numeric(gpData$map$chr)],type="p",axes=FALSE,...)
       axis(side=1,at=c(chr[!duplicated(gpData$map$chr)],max(chr)),labels=NA)
       axis(side=1,at=chr[!duplicated(gpData$map$chr)]+diff(c(chr[!duplicated(gpData$map$chr)],max(chr)))/2,tick=FALSE,labels=unique(gpData$map$chr),hadj=0)
       axis(side=2)
       box()
    }
}
