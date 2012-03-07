plotGenMap <- function (map, dense = FALSE, nMarker = TRUE, bw=1, centr=NULL, ...){
    oldPar <- par()
    if (class(map) == "gpData"){
       map.unit <- map$info$map.unit
       map <- map$map
    } else map.unit <- "unit"
    
    chr <- unique(map$chr)
    chr <- chr[!is.na(chr)]
    if(class(map$chr) == "factor") bord <- "transparent" else bord <- NULL
    map <- map[!is.na(map$chr), ]

    # centromere positions of maize
    if(!is.null(centr)) if(centr == "maize") centr <- c(133,90,95,104.6,105.5,50,55.3,46.5,68.8,59.9)

    # norm pos
    if (!is.null(centr)) map$pos <- map$pos - centr[map$chr]

    # add legend to the left side
    if (dense)  layout(matrix(2:1, ncol = 2), widths = c(0.82, 0.25))

    # colors from RColorBrewer red - green
    # display.brewer.pal(11, "RdYlGn")
    #cols <- c( "#FFFFBF","#FEE08B","#FDAE61","#F46D43","#D73027","#A50026")
    # display.brewer.pal(7, "Reds")
    cols <- c("#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#99000D")

    # compute density for a grid of values
    # cpmpute in advanve to use maxDens for legend
    if (dense) {
    x.grid <- y.grid <- list()
    maxDens <- 0
    for (i in seq(along = chr)) {
        start <- min(map$pos[map$chr == chr[i]], na.rm = TRUE)
        end <- max(map$pos[map$chr == chr[i]], na.rm = TRUE)
        x.grid[[i]] <- seq(from=start,to=end,by=bw)
        if(length(x.grid[[i]]) > 10000) warning("large discrepancy between map length and bandwith, maybe choose a larger value for 'bw'")
        y.grid[[i]] <- rep(NA,length(x.grid))
        for(j in seq(along=x.grid[[i]])){
           y.grid[[i]][j] <- sum(map$pos[map$chr == chr[i]] >= x.grid[[i]][j]-bw/2 & map$pos[map$chr == chr[i]] <= x.grid[[i]][j]+bw/2)
           if (y.grid[[i]][j]>maxDens) maxDens <- y.grid[[i]][j]
        }
      }
    }

    # add legend to the left margin of the plot
    if (dense) {

        par(mar = c(5, 1, 4, 3.8) + 0.1)
        image(seq(-0.3, 0.3, length = 20), seq(from = 0, to =  maxDens,
            length = 6), matrix(rep(seq(from = 0, to = maxDens, length = 6),
            20), nrow = 20, byrow = TRUE), col = cols, breaks=round(seq(0,maxDens,length=7)), axes = FALSE,
            xlab = "",main=paste("Nr. of SNPs \n within",bw,map.unit),xlim=c(-0.6,0.6))
        axis(side = 4, at = round(seq(from = 0, to = maxDens, length = 6))
            , labels=round(seq(from = 0, to = maxDens, length = 6)),las = 1)
        par(mar = c(5, 4, 4, 1) + 0.1)
    }

   # make an empty plot 
    if(!is.null(centr)) {
        plot(map, type = "n", xaxt = "n", xlim = c(0.5, length(chr) + 0.5), border = bord,
             ylim = c( max(map$pos,na.rm = TRUE) * 1.1, min(map$pos, na.rm = TRUE)),axes=FALSE, ...)
    } else{
        plot(map, type = "n", xaxt = "n", xlim = c(0.5, length(chr) + 0.5), border = bord,
             ylim = c( max(map$pos,na.rm = TRUE) * 1.1, min(map$pos, na.rm = TRUE)), ...)
    }    
   # x-axis     
    axis(side = 1, at = seq(along = chr), labels = chr)
   # y-axis
    if(!is.null(centr)){
        box()
        axis(side=2,at=-seq(-round(max(map$pos, na.rm = TRUE),-2),round(max(map$pos, na.rm = TRUE),-2),by=25),labels=abs(-seq(-round(max(map$pos, na.rm = TRUE),-2),round(max(map$pos, na.rm = TRUE),-2),by=25)),las=1)
    }

   # plot each chromosome
    for (i in seq(along = chr)) {

        if (dense) {
            image(seq(i - 0.35, i + 0.35, length = 20), x.grid[[i]],
                 matrix(rep(y.grid[[i]], 20), nrow = 20, byrow = TRUE),
                 col = cols, breaks=round(seq(0,maxDens,length=7)), add = TRUE)
            if(!is.null(centr)){
                # centromere
                polygon(x=c(i-0.4,i-0.1,i-0.1,i-0.4,i-0.4),y=c(-10,-1,1,10,-10),col="white",border="white")
                polygon(x=c(i+0.4,i+0.1,i+0.1,i+0.4,i+0.4),y=c(-10,-1,1,10,-10),col="white",border="white")
            }
        } else {
            n <- sum(map$chr == chr[i], na.rm = TRUE)
            start <- min(map$pos[map$chr == chr[i]], na.rm = TRUE)
            end <- max(map$pos[map$chr == chr[i]], na.rm = TRUE)
            lines(x = c(i, i), y = c(start, end))
            for (j in 1:n) {
                lines(x = c(i - 0.4, i + 0.4), y = rep(map$pos[map$chr == chr[i]][j], 2))
            }
        }
        # add nr. of markers
        if (nMarker)
            text(i, max(map$pos,na.rm=TRUE) * 1.05, sum(map$chr == chr[i],na.rm=TRUE))
    }
  oldPar$cin <- oldPar$cra <- oldPar$csi <- oldPar$cxy <- oldPar$din <- NULL
  par(oldPar)
}
