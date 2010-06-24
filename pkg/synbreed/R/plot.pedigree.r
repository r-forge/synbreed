plot.pedigree <- function(pedigree,effect=NULL,...){
  
  library(igraph)
  # gener has to start with 0
  
  if(min(pedigree$gener) != 0) pedigree$gener <- pedigree$gener-min(pedigree$gener)
  
  if (!is.null(effect) & length(effect)!=nrow(pedigree)) stop("length of effect does not equal nrow(pedigree)")
  

  gener <- pedigree$gener
  
  # set parents which are not in pedigree to unknown
  pedigree[!(pedigree$Par1 %in%  pedigree$ID),]$Par1 <- 0 
  pedigree[!(pedigree$Par2 %in%  pedigree$ID),]$Par2 <- 0 
    
  relations <- rbind(as.matrix(pedigree[pedigree$Par1!=0,c("Par1","ID")],ncol=2),as.matrix(pedigree[pedigree$Par2!=0,c("Par2","ID")],ncol=2))
  ped.graph <- graph.data.frame(relations,directed=TRUE,pedigree)
  # define x-y positions
  pos <- matrix(data=NA,nrow=nrow(pedigree),ncol=2)
 
  
  pos[,2] <- max(gener)-gener
  if (is.null(effect)) pos[,1] <- order(gener,partial=order(pedigree$Par1)) - cumsum(c(0,table(gener)))[gener+1]   
  else pos[,1] <- effect



  # plot
  plot(ped.graph,rescale=TRUE,vertex.label=pedigree$ID,layout=pos,edge.color=1,edge.width=0.5,edge.arrow.size=0.5,...)
  if (!is.null(effect)) axis(side=1,at=seq(-1,1,length=10),labels=round(seq(min(pos[,1]),max(pos[,1]),length=10),0))
}


