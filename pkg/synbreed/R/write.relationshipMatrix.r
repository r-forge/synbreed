write.relationshipMatrix <- function(relationshipMatrix,file=NULL,sorting="WOMBAT",ginv=TRUE,digits=10){
        
        if(sorting=="WOMBAT" & ginv==FALSE) stop("GINV must be specified for WOMBAT")
        
        # compute inverse relationship matrix
        if(ginv) rMinv <- ginv(relationshipMatrix)
        
        rMinv <- round(rMinv,digits)
        
        # add rownames and colnames
        res <- data.frame(coeff = as.numeric(rMinv),
                        rowname = rep(1:nrow(rMinv), nrow(rMinv)),
                        colname = rep(1:nrow(rMinv), each = nrow(rMinv)),
                        lower = as.logical(lower.tri(rMinv, diag = TRUE)))
                        
      
    
        # only use lower triangle
        res <- res[res$lower == TRUE, c("rowname", "colname", "coeff")]
        
        if (sorting=="WOMBAT"){
          res <- res[,c(2,1,3)]
          res <-  res[order( res$rowname,  res$colname), ] 
        }
        else res <-  res[order( res$colname,  res$rowname), ]  
        
        # write to given file
        if (!is.null(file)) write.table(res, file, sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

        else return(res)


}