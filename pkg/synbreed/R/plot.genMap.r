plotGenMap <- function(map,dense=FALSE,nMarker=TRUE,...){
   
    if(class(map)=="gpData") map <- map$map
   
    # extract chromosomes
    chr <- unique(map$chr)
    # without NA
    chr <- chr[!is.na(chr)]
    
    
    # initialize map
    if(dense)layout(matrix(1:2,ncol=2),width=c(0.82,0.18))
    plot(map,type="n",xaxt="n",xlim=c(0.5,length(chr)+0.5),ylim=c(min(map$pos,na.rm=TRUE),max(map$pos,na.rm=TRUE)*1.1),...) 
    axis(side=1,at=seq(along=chr),labels=chr)
    
    # loop over chromosomes
    for (i in seq(along=chr)){
    
        n <- sum(map$chr==chr[i],na.rm=TRUE) 
        start <- min(map$pos[map$chr==chr[i]],na.rm=TRUE)
        end <- max(map$pos[map$chr==chr[i]],na.rm=TRUE)
  
      if(dense){ # kernel density estimation

        # fixed bandwith for density estimation
        bw <- (end-start)/n

        # density estimation
        densEst <- density(map$pos[map$chr==chr[i]],kernel="rectangular",from=min(map$pos[map$chr==chr[i]],na.rm=TRUE),to=max(map$pos[map$chr==chr[i]],na.rm=TRUE),cut=TRUE,bw=bw)
             
        # norm to 1
        densEst$y <- densEst$y/max(densEst$y)
        
        # from RColorBrewer, brewer.pal(11,"RdYlGn")
        cols <- c("#A50026","#D73027","#F46D43","#FDAE61","#FEE08B","#FFFFBF","#D9EF8B","#A6D96A","#66BD63","#1A9850","#006837")[11:1]

        # visualisation of map
        image(seq(i-0.4,i+0.4,length=20),densEst$x,matrix(rep(densEst$y,20),nrow=20,byrow=TRUE),col=cols,add=TRUE)
        
        } # end of if,  kernel denisty estimation
      else{
        lines(x=c(i,i),y=c(start,end))
        for(j in 1:n){
          lines(x=c(i-0.4,i+0.4),y=rep(map$pos[map$chr==chr[i]][j],2))
        }
      } # end of else
 
      if(nMarker) text(i,max(map$pos)*1.05,sum(map$chr==chr[i]))
 
      } # end chromosome loop
    
    # add legend
    if (dense){
      par(mar=c(5,0,4,3.8)+0.1)
      image(seq(-0.4,0.4,length=20),seq(from=0,to=1,length=11),matrix(rep(seq(from=0,to=1,length=11),20),nrow=20,byrow=TRUE),col=cols,axes=FALSE,xlab="")
      axis(side=4,at=round(seq(from=0,to=1,length=11),4),las=1)
      par(mar=c(5,4,4,1)+0.1)
    }
}

