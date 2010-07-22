LDMap <- function(marker,linkageGroup,pos,file=NULL,...){
     # Calculating LD as r2
     ld.r <- cor(marker,method="spearman")^2
     lg <- unique(linkageGroup)
     for (i in 1:length(lg)){
      LDheatmap(ld.r[linkageGroup == lg[i],linkageGroup == lg[i]], LDmeasure="r",title=paste("Pairwise LD r2, Linkage Group ",i,sep=""), color=brewer.pal(9,"OrRd")[9:1], genetic.distances=pos[linkageGroup == lg[i]], geneMapLabelY=0.12, geneMapLabelX=0.35, distances = "genetic",...)
     }
     if(!is.null(file))
}

data(maize)
LDMap(maize.geno[,-1],maize.marker.pos$chr,maize.marker.pos$pos)
