# relationship coefficient based on rogers distance and transformed to
# relationship coefficient  according to formula of Melchinger (1991) or
# Hayes and Goddard (2008)


rogers <- function(marker,correction=c("none","Hayes","Melchinger")){
          if(!is.matrix(marker)) marker <- as.matrix(marker)
          if(is.null(correction)) correction <- "none"
          else correction <- match.arg(correction,c("none","Hayes","Melchinger"))
          # code marker to -1/0/1
          marker <- marker - (max(marker,na.rm=TRUE)-1)
          m <- ncol(marker)
          d <- 1- (tcrossprod(marker) + m)/(2*m)
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


