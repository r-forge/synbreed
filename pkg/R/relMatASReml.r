# function to compute inverse relationshipmatrix and save it in format suitable for ASReml

relMatASReml <- function(relationshipMatrix,ginv=FALSE,file=NULL,digits=10){
        # compute inverse relationship matrix
        if(ginv) rMinv <- ginv(relationshipMatrix)
        else  rMinv <- solve(relationshipMatrix)
        rMinv <- round(rMinv,digits)
        
        # add rownames and colnames
        res <- data.frame(coeff = as.numeric(rMinv),
                        rowname = rep(1:nrow(rMinv), nrow(rMinv)),
                        colname = rep(1:nrow(rMinv), each = nrow(rMinv)),
                        lower = as.logical(lower.tri(rMinv, diag = TRUE)))
                        
        # order res first by rowname and then by colnames
         res <-  res[order( res$rowname,  res$colname), ]
        # only use lower triangle
        res <- res[res$lower == TRUE, c("rowname", "colname", "coeff")]
        
        # write to given file
        if (!is.null(file)) write.table(res, file, sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

        else return(res)
}