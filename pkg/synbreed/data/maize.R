maize.geno <- read.table("maize.geno.tab",header=TRUE,sep=" ")
maize.pheno <- read.csv2("maize.pheno.csv",header=TRUE,sep=",",dec=".")
maize.ped <- read.table("maize.ped.tab",header=TRUE,sep="\t")
maize.marker.pos <- read.table("maize.marker.pos.tab",header=TRUE,sep="\t")