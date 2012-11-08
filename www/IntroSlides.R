### R code from vignette source 'IntroSlides.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: IntroSlides.Rnw:182-186
###################################################
set.seed(1234)
library(synbreed)
data(maize)
data(mice)


###################################################
### code chunk number 2: IntroSlides.Rnw:200-201 (eval = FALSE)
###################################################
## install.packages("synbreed")


###################################################
### code chunk number 3: IntroSlides.Rnw:204-205 (eval = FALSE)
###################################################
## install.packages("synbreed",repos="http://r-forge.r-project.org")


###################################################
### code chunk number 4: IntroSlides.Rnw:210-211 (eval = FALSE)
###################################################
## library(synbreed)


###################################################
### code chunk number 5: IntroSlides.Rnw:220-221 (eval = FALSE)
###################################################
## vignette("IntroSyn")	


###################################################
### code chunk number 6: IntroSlides.Rnw:225-226 (eval = FALSE)
###################################################
## help(package="synbreed")


###################################################
### code chunk number 7: IntroSlides.Rnw:229-230 (eval = FALSE)
###################################################
## demo(package="synbreed")


###################################################
### code chunk number 8: IntroSlides.Rnw:241-242
###################################################
citation(package="synbreed")


###################################################
### code chunk number 9: IntroSlides.Rnw:300-301 (eval = FALSE)
###################################################
## gp <- create.gpData(pheno,geno,map,pedigree,covar,map.unit="cM")


###################################################
### code chunk number 10: IntroSlides.Rnw:325-329
###################################################
id <- c("A","B","C","D","E")
par1 <- c(0,0,"A","A","D")
par2 <- c(0,0,"B","C","B") 
(ped <- create.pedigree(id,par1,par2))


###################################################
### code chunk number 11: IntroSlides.Rnw:350-360 (eval = FALSE)
###################################################
## # Read file TrueEBV.txt with pedigree, trait, and tbv
## dat <- read.table("TrueEBV.txt",header=TRUE,stringsAsFactors=FALSE)
## # Create object of class 'pedigree'
## ped <- with(dat,create.pedigree(ID=id,Par1=sire,Par2=dam,gener=gen,sex=abs(sex-2)))
## # Phenotypic data
## pheno <- data.frame(trait=dat$Phenotype,row.names=dat$id) 
## # covar = tbv
## covar <- data.frame(tbv=dat$GeneticValue,row.names=dat$id)
## # genotypic data
## geno <- read.table("genotype_cor.txt",header=FALSE,stringsAsFactors=FALSE)


###################################################
### code chunk number 12: IntroSlides.Rnw:368-385 (eval = FALSE)
###################################################
## # gametes to genotypes
## geno2 <- matrix(data=NA,nrow=nrow(geno),ncol=(ncol(geno)-1)/2)
## for (j in 1:ncol(geno2)){
##   # combine phased data to a genotype
##   geno2[,j] <- paste(as.character(geno[,2*j]),as.character(geno[,2*j+1]),sep="")
## }
## # create map
## # 6 chromosomes with 1000 markers
## # dist between adjacent markers = 0.1cM
## chr <- rep(1:6,each=1000)
## pos <- rep(seq(from=0,to=99.9,by=.1),times=6)
## map <- data.frame(chr=chr,pos=pos)
## # create gpData object
## qtlMASdata <- create.gpData(pheno=pheno,geno=geno2,map=map,pedigree=ped,covar=covar,map.unit="cM")
## # save data as object of class gpData in Rdata-format
## save("qtlMASdata",file="qtlMASdata.Rdata")
## # for loading data, function load() and ls() might be useful


###################################################
### code chunk number 13: IntroSlides.Rnw:433-436
###################################################
library(synbreed)
data(maize)
summary(maize)


###################################################
### code chunk number 14: IntroSlides.Rnw:443-444 (eval = FALSE)
###################################################
## str(maize)


###################################################
### code chunk number 15: IntroSlides.Rnw:448-449
###################################################
head(maize$pheno[,1,])


###################################################
### code chunk number 16: IntroSlides.Rnw:452-453
###################################################
maize$geno[10:13,20:25]


###################################################
### code chunk number 17: IntroSlides.Rnw:460-461
###################################################
head(maize$covar,n=4)


###################################################
### code chunk number 18: IntroSlides.Rnw:469-470 (eval = FALSE)
###################################################
## maize$covar$id[maize$covar$phenotyped]


###################################################
### code chunk number 19: IntroSlides.Rnw:482-484
###################################################
maizeChr1to5 <- discard.markers(maize,rownames(maize$map)[maize$map$chr > 5])
summary(maizeChr1to5$map)


###################################################
### code chunk number 20: IntroSlides.Rnw:492-493
###################################################
plotGenMap(maize)


###################################################
### code chunk number 21: IntroSlides.Rnw:501-502
###################################################
plotGenMap(mice,dense=TRUE,nMarker = FALSE, bw=1)


###################################################
### code chunk number 22: IntroSlides.Rnw:509-510
###################################################
summaryGenMap(maize)


###################################################
### code chunk number 23: IntroSlides.Rnw:567-569
###################################################
maizeC <- codeGeno(maize,maf=0.05,nmiss=0.1,
verbose=TRUE)


###################################################
### code chunk number 24: IntroSlides.Rnw:572-573 (eval = FALSE)
###################################################
## maizeLD <- pairwiseLD(maizeC,chr=1,type="data.frame")


###################################################
### code chunk number 25: IntroSlides.Rnw:581-585 (eval = FALSE)
###################################################
## codeGeno(gpData, impute = FALSE, impute.type = c("fix", 
##      "random", "family", "Beagle", "BeagleAfterFamily"), 
##      replace.value = NULL, maf = NULL, nmiss = NULL, label.heter = "AB", 
##       keep.identical = TRUE, verbose = FALSE)


###################################################
### code chunk number 26: IntroSlides.Rnw:611-612 (eval = FALSE)
###################################################
## plot(maizeLD); plot(maizeLD,type="bars")


###################################################
### code chunk number 27: IntroSlides.Rnw:683-684 (eval = FALSE)
###################################################
## kin(gpData,ret="add")


###################################################
### code chunk number 28: IntroSlides.Rnw:687-688 (eval = FALSE)
###################################################
## kin(gpData,ret="dom")


###################################################
### code chunk number 29: IntroSlides.Rnw:691-692 (eval = FALSE)
###################################################
## kin(gpData,ret="kin")


###################################################
### code chunk number 30: IntroSlides.Rnw:695-696 (eval = FALSE)
###################################################
## kin(gpData,ret="gam")


###################################################
### code chunk number 31: IntroSlides.Rnw:701-702 (eval = FALSE)
###################################################
## A <- kin(maizeC,ret="kin",DH=maize$covar$DH)


###################################################
### code chunk number 32: IntroSlides.Rnw:721-723 (eval = FALSE)
###################################################
## U <- kin(maizeC,ret="realized")/2
## plot(A[maize$covar$genotyped,maize$covar$genotyped]); plot(U)


###################################################
### code chunk number 33: IntroSlides.Rnw:854-857
###################################################
library(xtable)
dat <- data.frame(y=c(132,147,156,172),time=c(1,2,1,2),animal=c(1,2,3,4))
print(xtable(dat,digits=0),include.rownames=FALSE)


###################################################
### code chunk number 34: IntroSlides.Rnw:861-865
###################################################
library(synbreed)
ped <- create.pedigree(ID=c(6,5,1,2,3,4),Par1=c(0,0,5,5,1,6),Par2=c(0,0,0,0,6,2))
par(mar=c(0,0,0,0))
plot(ped,c(-.5,.5,2,-2,1,-1))


###################################################
### code chunk number 35: IntroSlides.Rnw:1026-1032
###################################################
dat <- data.frame(y=c(132,147,156,172),time=c(1,2,1,2),animal=c(1,2,3,4))
ped <- create.pedigree(ID=c(6,5,1,2,3,4),Par1=c(0,0,5,5,1,6),Par2=c(0,0,0,0,6,2))
gp <- create.gpData(pheno=dat,pedigree=ped)
A <- kin(gp,ret="add")
(X <- matrix(c(1,0,1,0,0,1,0,1),ncol=2))
(Z <- diag(6)[-c(1,2),])


###################################################
### code chunk number 36: IntroSlides.Rnw:1037-1042
###################################################
(AI <- solve(A))
RI <- diag(4)

res <- MME(X,Z,AI*3,RI,dat$y)
res$b; res$u


###################################################
### code chunk number 37: IntroSlides.Rnw:1051-1053 (eval = FALSE)
###################################################
## modA <- gpMod(maizeC,model="BLUP",kin=A)
## modU <- gpMod(maizeC,model="BLUP",kin=U)


###################################################
### code chunk number 38: IntroSlides.Rnw:1057-1059 (eval = FALSE)
###################################################
## gA <- predict(modA)
## gU <- predict(modU)


###################################################
### code chunk number 39: IntroSlides.Rnw:1062-1063 (eval = FALSE)
###################################################
## tbv <- maizeC$covar$tbv[maizeC$covar$phenotyped]


###################################################
### code chunk number 40: IntroSlides.Rnw:1132-1134 (eval = FALSE)
###################################################
## last50 <- rownames(maizeC$pheno)[1201:1250]
## maizeC2 <- discard.individuals(maizeC,last50)


###################################################
### code chunk number 41: IntroSlides.Rnw:1137-1138 (eval = FALSE)
###################################################
## modU24 <- gpMod(maizeC2,model="modU",kin=U)


###################################################
### code chunk number 42: IntroSlides.Rnw:1142-1143 (eval = FALSE)
###################################################
## g <- predict(modU24,rownames(maizeC$pheno)[1201:1250])


###################################################
### code chunk number 43: IntroSlides.Rnw:1165-1166 (eval = FALSE)
###################################################
##  cv.maize <- crossVal(maizeC,cov.matrix=list(U),k=5,Rep=2,Seed=123,sampling="random",varComp=modU$fit$sigma)


###################################################
### code chunk number 44: IntroSlides.Rnw:1205-1210
###################################################
y <- maize$pheno[,1,]
X <- maize$geno
sX2 <- sum(X^2)
h2 <- 0.5 # priori expectation
(lambdaStart <- sqrt(2*sum(X^2)*(1-h2)/h2/nrow(X)))


###################################################
### code chunk number 45: IntroSlides.Rnw:1216-1220
###################################################
lambda <- seq(from=0,to=100,by=1)
dens <- dgamma(x=lambda^2,shape=.52,rate=3e-5)*lambda*2 # distribution for lambda 
plot(dens~lambda,type='l',ylab="density")
abline(v=lambdaStart)


###################################################
### code chunk number 46: IntroSlides.Rnw:1226-1228 (eval = FALSE)
###################################################
## prior <- list(varE=list(df=3,S=35),lambda = list(shape=0.52,rate=1e-4,value=lambdaStart,type='random'))
## modBL <- gpMod(maizeC,model="BL",prior=prior,nIter=6000,burnIn=1000,thin=5)


###################################################
### code chunk number 47: IntroSlides.Rnw:1231-1233 (eval = FALSE)
###################################################
## cv.BL <- crossVal(maizeC,k=5,Rep=2,Seed=123,sampling="random",VC.est="BL",prior=prior)
## summary(cv.BL)


###################################################
### code chunk number 48: IntroSlides.Rnw:1243-1244 (eval = FALSE)
###################################################
## gpData2cross(gpDataObj)


###################################################
### code chunk number 49: IntroSlides.Rnw:1247-1248 (eval = FALSE)
###################################################
## cross2gpData(crossObj)


