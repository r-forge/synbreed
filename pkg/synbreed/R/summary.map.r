# summary for maker maps

summaryGenMap <- function(map){

     # information from arguments
     if(class(map)=="gpData") map <- map$map
     chr <- map$chr
     pos <- map$pos

     # extract information
     # number of markers
     len <- tapply(pos,chr,length)
     rge <- tapply(pos,chr,max,na.rm=TRUE)-tapply(pos,chr,min,na.rm=TRUE)

     # differences of markers on each chr
     diffs <- tapply(pos[!is.na(pos)],chr[!is.na(pos)],diff,na.rm=TRUE)
     avDist <- as.numeric(lapply(diffs,mean,na.rm=TRUE))
     maxDist <- as.numeric(lapply(diffs,max,na.rm=TRUE))
     minDist <- as.numeric(lapply(diffs,min,na.rm=TRUE))

     # return data.frame
     ret <- data.frame(noM = len, range=rge, avDist = avDist, maxDist = maxDist, minDist = minDist,row.names=names(len))
     return(ret)

}
