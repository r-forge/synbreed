\name{crossVal}
\alias{crossVal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross validation of different prediction models
}
\description{Function for the application of the cross validation procedure on prediction 
models with fixed and random effects. Covariance matrices must be committed to the function 
and variance components can be committed or reestimated with ASReml or the BLR function.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
crossVal(y, X, Z, cov.matrix = NULL,  k = 2, Rep = 1, Seed = NULL, 
         sampling = c("random", "within family", "across family"), 
         varComp = NULL, popStruc = NULL, VC.est = c("commit",
         "ASReml","BRR","BL"),priorBLR=NULL,verbose=TRUE,nIter=5000,
	 burnIn=1000,thin=10)
}
\arguments{
  \item{y}{Phenotypic records of class \code{matrix} with two columns, where the first column contains IDs of individuals and the second column contains the phenotypic values.
}
  \item{X}{Design matrix for the fixed effects, where the colnames represents the names of the fixed effects.
}
  \item{Z}{Design matrix for the random effects, where the colnames represents the names of the random effects.
}
  \item{cov.matrix}{A object of class \code{list} including covariance matrices for the random effects. Size and order of rows and columns should be equal to rownames of \code{y}. 
If no covariance is given, the identity matrix is used.
}
  \item{k}{Number of folds for k-fold cross validation, thus \code{k} should be in [2,\code{nrow(y)}].
}
  \item{Rep}{Number of replications.
}
  \item{Seed}{Number for \code{set.seed()} to make results reproducable.
}
  \item{sampling}{Different sampling strategies can be "\code{random}", "\code{within family}" or "\code{across family}".
}
  \item{varComp}{A \code{vector} of variance components for the random effects, which has to be specified if \code{VC.est="commit"}. The first variance components should be 
the same order as the given covariance matrices, the last given variance component is for the residuals.
}
  \item{popStruc}{Vector of length \code{nrow(y)} assigning individuals to families.
}
  \item{VC.est}{Should variance components be reestimated with "\code{"ASReml"}" or with a Bayesian approach "\code{"BRR"}" and "\code{"BL"}" within the estimation set of each fold in the cross validation? If
\code{VC.est="commit"}, the variance components have to be defined in \code{varComp}. For "\code{ASReml}", ASReml software have to be installed on the system. 
}
  \item{priorBLR}{A vector with priors for varBR and varE within the BLR function. For \code{VC.est="BRR"}, varBR and varE have to be specified, for \code{VC.est="BL"}, varE has to be specified.
}
  \item{verbose}{Logical. Wether output shows replications and folds.
}
  \item{nIter, burnIn, thin}{Number of iterations, burn-in and thinning for \code{VC.est="BRR"} and \code{VC.est="BL"} in the BLR-function
}
}
\details{
In cross validation the data set is splitted into an estimation (ES) and a test set (TS). The effects are estimated with the ES to predict observations in the TS.
For sampling into ES and TS, k-fold cross validation is applied, where the data set is splitted into k subsets and k-1 comprising the ES and 1 is the TS, repeated for each subset.

To account for the family structure, \code{sampling} can be defined as:
\describe{
\item{random}{ Does not account for family structure, random sampling within the complete data set}
\item{within family}{ Accounts for within family information, each family is splitted into k subsets}
\item{across family}{ Accounts for across family information, ES and TS contains a set of complete families}}
The following mixed model equation is used for \code{VC.est="commit"}:

\deqn{\bf y=\bf{Xb}+\bf{Zu}+\bf e}{y=Xb+Zu+e} with \deqn{\bf u \sim N(0,G\sigma^2_u)}{u=N(0,Gsigma2u)}

\deqn{\left(\begin{array}{cc} \bf X'\bf X & \bf X'\bf Z \\ \bf Z'\bf X & \bf Z'\bf Z + \bf G^{-1}\frac{\sigma^2_e}{\sigma^2_u} \end{array}\right) \left(\begin{array}{c} \bf b \\ \bf  u \end{array}\right) = \left(\begin{array}{c}\bf  X'\bf y \\ \bf Z'\bf y \end{array}\right)}{(X'X,X'Z,Z'X,ZZ+GIsigma2e/sigma2u)(b,u)=(X'y,Z'y)}
}
\value{
A object of class \code{list} with follwoing items:
\item{bu }{Estimated fixed and random effects of each fold within each replication.}
\item{y.TS}{Predicted values of all test sets within each replication.}
\item{PredAbi}{Predictive ability of each fold within each replication calculated as \eqn{r(y_{TS},\hat y_{TS})} .}
\item{bias}{Regression coefficients of a regression of the observed values on the predicted values in the TS. 
A regression coefficient \eqn{< 1} implies extreme high (low) values of the predicted values over- (under-) estimated the observed values, and vice versa for a regression coefficient \eqn{> 1}.}
\item{k}{Number of folds}
\item{Rep}{Replications} 
\item{sampling}{Sampling method}
\item{Seed}{Seed for \code{set.seed()}}
\item{rep.seed}{Calculated seeds for each replication}
\item{nr.ranEff}{Number of random effects} 
\item{VC.est.method}{Method for the variance components ("\code{committed}" or "\code{reestimated with ASReml/BRR/BL}")}
}
\references{
Legarra A, Robert-Granie C, Manfredi E, Elsen J (2008) Performance of genomic selection in mice. Genetics 180:611-618

Luan T, Wooliams JA, Lien, S, Kent M, Svendsen M, Meuwissen THE (2009) The accuracy of genomic selection in Norwegia red cattle assessed by crosss-validation. Genetics 183:1119-1126

Mosier CI (1951) I. Problems and design of cross-validation 1. Educ Psychol Measurement 11:5-11

Crossa J, de los Campos G, Perez P, Gianola D, Burgueno J, et al. (2010) Prediction of genetic values of quantitative traits in plant breeding using pedigree and molecular markers, Genetics 186:713-724

}
\author{Theresa Albrecht
}
\seealso{
\code{\link{summary.cvData}}}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# loading the maize data set
data(maize)
# prepare matrices for cross validation
maize2 <-codeGeno(maize)
maize2$pheno <- data.frame(rownames(maize2$pheno),maize2$pheno[,1])
rad<-vanRaden(maize2$geno)
diag(rad)<-diag(rad)+0.00001 # to avoid singularities
X <- matrix(rep(1,nrow(maize2$pheno)),ncol=1)
Z <- diag(nrow(maize2$pheno))
# cross validation
cv.maize <- crossVal(maize2$pheno,X,Z,cov.matrix=list(rad),k=5,Rep=1,
            Seed=123,sampling="random",varComp=c(26.5282,48.5785),VC.est="commit")
cv.maize2 <- crossVal(maize2$pheno,X,Z=maize2$geno,k=5,Rep=1,
             Seed=123,sampling="random",varComp=c(0.0704447,48.5785),VC.est="commit")
# comparing results, both are equal!
cv.maize$PredAbi
cv.maize2$PredAbi
summary(cv.maize) 
summary(cv.maize2) 
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

