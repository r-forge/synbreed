LDMap <- function(LDmat,gpData,file=NULL,...){

    # catch (possible) errors
    if(class(LDmat)!="LDmat") stop("'LDmat' must be of class 'LDmat'")
    lg <- 1:length(LDmat$LD)
    pos <- gpData$map$pos
    names(pos) <- rownames(gpData$map)

     # use LD from input arguement
    ret <- LDmat
     
     if(!is.null(file)) pdf(file)
      for (i in 1:length(lg)){
     
        color = c("#7F0000","#B30000","#D7301F","#EF6548","#FC8D59","#FDBB84","#FDD49E","#FEE8C8","#FFF7EC")
        #   using function LDheatmap
        MapUnit <- ifelse(gpData$info$map.unit=="cM","genetics","physical")
        LDheatmap(LDmat$LD[[i]], LDmeasure="r", color=color, genetic.distances=pos[rownames(LDmat$LD[[i]])],distances=MapUnit,   geneMapLabelY=0.12, geneMapLabelX=0.35,...)
         
      }
     if(!is.null(file)) dev.off()
 
}