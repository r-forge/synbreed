#############################################
## Create a genomic prediction Data object
## from the wheat data in the BLR package
##
##
## author : Valentin Wimmer
## date : 2011 - 11 - 30
##
##############################################

data(wheat)
# X = genotypes
# Y = phenotypes

# assign names
rownames(X) <- paste("ID",1000+(1:nrow(X)),sep="")
colnames(Y) <- paste("Env",1:4,sep="")
rownames(Y) <- paste("ID",1000+(1:nrow(X)),sep="")

# create a gpData object
gpWheat <- create.gpData(pheno=Y,geno=X)
gpWheat <- codeGeno(gpWheat)



# predictive ability using Bayesian Lasso
# use prior values from Crossa et al. (2010)
priorCrossa <- list(varE=list(df=4,S=1),lambda = list(shape=0.6,rate=1e-4,value=20,type='random'))

modBL <- gpMod(gpWheat,trait=1,model="BL",prior=priorCrossa)

cv <- crossVal(gpWheat,trait=1,VC.est="BL",priorBLR=priorCrossa,k=10,Rep=1)
#Object of class 'cvData'
#
# 10 -fold cross validation with 1 replications
#     Sampling:                 random
#     Variance components:      reestimated with BL
#     Number of random effects: 0
#     Number of individuals:    599 -- 599
#     Size of the TS:           59 -- 60
#
#Results:
#                      Min         Mean +- pooled SE       Max
# Predictive ability:  0.2897      0.5212 +- NA            0.6788
# Rank correlation:    0.2989      0.4760 +- NA            0.5985
# Bias:                0.5041      1.0505 +- NA            1.3903
# 10% best predicted:  0.56        0.56 +- NA      0.56
#
#Seed start:
#Seed replications:
#[1] 4178
