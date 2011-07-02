plotNeighbourLD <- function(LD,map,nMarker=TRUE,dense=FALSE,...){

    if(class(map)=="gpData") map <- map$map
   
    # extract chromosomes
    chr <-  (1:length(LD$LD))[!as.logical(lapply(LD$LD,is.null))]
    
    
    # initialize map
    layout(matrix(1:2,ncol=2),width=c(0.82,0.18))
    plot(map,type="n",xaxt="n",xlim=range(chr,na.rm=TRUE)+c(-0.5,0.5),ylim= c(max(map$pos,na.rm = TRUE) * 1.1, min(map$pos, na.rm = TRUE)),...) 
    axis(side=1,at=chr,labels=chr)
    
    # loop over chromosomes
    for (i in chr){
    
        n <- sum(map$chr==i,na.rm=TRUE) 
        start <- min(map$pos[map$chr==i],na.rm=TRUE)
        end <- max(map$pos[map$chr==i],na.rm=TRUE)

        # from RColorBrewer, brewer.pal(11,"RdYlGn")
        #cols <- c("#A50026","#D73027","#F46D43","#FDAE61","#FEE08B","#FFFFBF","#D9EF8B","#A6D96A","#66BD63","#1A9850","#006837")[11:1] 
        
        # display.brewer.pal(7, "Reds")s
        cols <- c("#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#99000D")

	if(dense){  # calculating averaged LD
		# map posittions
		mapPos <- smooth(map$pos[map$chr==i])
		# matrices with 10 rows
		avLD <- matrix(c(diag(LD$LD[[i]][-nrow(LD$LD[[i]]),-1]),rep(NA,(ceiling((n-1)/10)*10-(n-1)))),nrow=10)
		avPos <- matrix(c(mapPos[-length(mapPos)],rep(NA,(ceiling((n-1)/10)*10-(n-1)))),nrow=10)

        	# visualisation of map
        	image(seq(i-0.4,i+0.4,length=20), colMeans(avPos,na.rm=TRUE),matrix(rep(colMeans(avLD,na.rm=TRUE),20),nrow=20,byrow=TRUE),col=cols,add=TRUE)

	}  # end of if,  smooth LD calculation
	else{  # using LD directly
		# visualisation of map
        	image(seq(i-0.4,i+0.4,length=20), map$pos[map$chr==i],matrix(rep(diag(LD$LD[[i]][-nrow(LD$LD[[i]]),-1]),20),nrow=20,byrow=TRUE),col=cols,add=TRUE)
	}
      if(nMarker) text(i,max(map$pos)*1.05,sum(map$chr==i))
 
      } # end chromosome loop
    
    # add legend
      par(mar=c(5,0,4,3.8)+0.1)
      image(seq(-0.4,0.4,length=20),seq(from=0,to=1,length=11),matrix(rep(seq(from=0,to=1,length=11),20),nrow=20,byrow=TRUE),col=cols,axes=FALSE,xlab="")
      axis(side=4,at=round(seq(from=0,to=1,length=11),4),las=1)
      par(mar=c(5,4,4,1)+0.1)
}

