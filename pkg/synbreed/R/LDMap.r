LDMap <- function(gpData,chr=NULL,file=NULL,...){

    # catch (possible) errors
    if(is.null(gpData$geno)) stop("no genotypic data available")
    if(!gpData$info$codeGeno) stop("use function 'codeGeno' before")
    if(is.null(gpData$map))  stop("no map information available")


     # extract information from gpData 
     mapped <- !(is.na(gpData$map$chr) & is.na(gpData$map$pos)) 
     marker <- gpData$geno[,mapped]
     linkageGroup <- gpData$map$chr[mapped]
     pos <- gpData$map$pos[mapped]

    # select chromosomes if 'chr' is specified
    lg <- unique(linkageGroup)
    if(!is.null(chr)){ 
        lg <- chr
        if(chr=="all") linkageGroup <- rep("all",length(linkageGroup))

    }

     # Calculating LD as r2
     ld.r2 <- cor(marker,method="spearman")^2

     if(!is.null(file)) pdf(file)
     for (i in 1:length(lg)){
      LDheatmap(ld.r2[linkageGroup == lg[i],linkageGroup == lg[i]], LDmeasure="r",title=paste("Pairwise LD r2, Linkage Group ",i,sep=""), color=brewer.pal(9,"OrRd")[9:1], genetic.distances=pos[linkageGroup == lg[i]], geneMapLabelY=0.12, geneMapLabelX=0.35, distances = "genetic",...)
     }
     if(!is.null(file)) dev.off()
     
     # return matrix of r2
     colnames(ld.r2) <- rownames(ld.r2) <- colnames(marker)
     invisible(ld.r2)
}