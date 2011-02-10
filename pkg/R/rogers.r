# similarity coefficient based on rogers distance and transformed to
# relationship coefficient  according to formula of Melchinger (1991) or
# Hayes and Goddard (2008)

rogers <- function(marker,correction=c("none","Hayes","Melchinger")){
          
          # extract information from arguments
          if(any(class(marker)=="gpData")){
             if(!marker$info$codeGeno) stop("use function 'codeGeno' before using 'rogers'") 
             marker <- marker$geno
          }
          else marker <- marker
          
          # should be a matrix
          if(!is.matrix(marker)) marker <- as.matrix(marker)
          
          # no correction if not specified
          if(is.null(correction)) correction <- "none"
          else correction <- match.arg(correction,c("none","Hayes","Melchinger"))
          
          # code marker to -1/0/1 from 0,1,2
          marker <- marker - (max(marker,na.rm=TRUE)-1)
          m <- ncol(marker)
          
          # rogers distance
          d <- 1- (tcrossprod(marker) + m)/(2*m)
          
          # different corrections
          if(correction=="Melchinger") f <- 1-d/mean(d,na.rm=TRUE)
          if(correction=="Hayes"){
            s <- 1-d
            smin <- min(s,na.rm=TRUE)
            f <- (s-smin)/(1-smin)
          }
          else f <- 1-d
          
          # mulitply with 2 to get relationship and not kinship
          class(f) <- "relationshipMatrix"
          return(2*f)
}


