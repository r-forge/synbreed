cross2gpData <- function(cross,...){
    # check for class
    # sometimes multiple classes in package 'qtl'
    if(all(class(cross)!="cross")) stop("object '",substitute(cross),"' not of class 'cross'")  
    # phenotypic data
    pheno <- cross$pheno
    # names of chromosomes
    chr <- names(cross$geno)
    # extract informations from first chromosome
    geno <- cross$geno[[chr[1]]]$data
    map <- data.frame(pos=cross$geno[[chr[1]]]$map,chr=chr[1])
    # loop over chromosomes 2:n
    if (length(chr)>1){
       for(i in chr[-1]){
          # add new genotypes
          geno <- cbind(geno,cross$geno[[i]]$data)
          # add new map info
          map <- rbind(map,data.frame(pos=cross$geno[[i]]$map,chr=i))
       }
    }
    # combine pheno, geno, and map  in class 'gpData'
    ret <- create.gpData(pheno=pheno,geno=geno,map=map)
    # use codeGeno
    # cat("using function 'codeGeno'\n")
    ret <- codeGeno(ret,label.heter="1",...)
    return(ret)
} 