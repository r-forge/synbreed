\name{create.pedigree}
\alias{create.pedigree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to create an object of class pedigree
}
\description{
This function can be used to create a \code{pedigree} object. This function allows inlcuding inbreeds in the pedigree.
}
\usage{
create.pedigree(ID, Par1, Par2, gener)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ID}{
vector of unique IDs identifying e.g. each genotype. Should be numeric.
}
  \item{Par1}{
vector of IDs identifying parent 1 (father)
}
  \item{Par2}{
vector of IDs identifying parent 1 (mother)
}
  \item{gener}{
vector identifying the generation
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class \code{pedigree} with the same contents as above. Note that gener starts from 0.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Valentin Wimmer
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{plot.pedigree}
}
\examples{
ped.raw <- matrix(c(1,  0,   0,   1, 
                 2,  0,   0,   1,
                 3,  0,   0,   1,
                 4,  0,   0,   1,
                5,  1,   2,   2, 
                 6,  1,   3,   2,
                 7,  1,   2,   2,
                 8,  1,   2,   2) ,
                 byrow = TRUE, ncol = 4)
ped <- create.pedigree(ped.raw[,1],ped.raw[,2],ped.raw[,3],ped.raw[,4])
ped
}