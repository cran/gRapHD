\name{calcStat}
\alias{calcStat}
\title{Pairwise weights}
\description{
  Calculates pairwise statistics (-2*log-LR, AIC, or BIC) for each variable pair (edge) in the dataset.
}

\usage{
calcStat(dataset,homog=TRUE,forbEdges=NULL,stat="LR")
}

\arguments{
  \item{dataset}{matrix or data frame (\code{nrow(dataset)} observations and
                 \code{ncol(dataset)} variables).}
  \item{homog}{\code{TRUE} for homogeneous covariance structure, \code{FALSE}
               for heterogeneous. This is only meaningful with mixed models.
               Default is homogeneous (\code{TRUE}).}
  \item{forbEdges}{list with edges that should not be considered. Matrix with 2
                   columns, each row representing one edge, and each column one
                   of the vertices in the edge. Default is \code{NULL}.}
  \item{stat}{measure to be minimized: LR (-2*log-likelihood), AIC, or BIC.
              Default is LR. It can also be a user defined function with 
              format: \code{FUN(newEdge,}\code{numCat,}\code{dataset)}; where 
              \code{numCat} is a vector with number of levels for each variable 
              (0 if continuous); \code{newEdge} is a vector with length two; 
              and \code{dataset} is a matrix (n by p).}
}

\value{
  A matrix with \code{p(p-1)/2} lines and \code{4} columns, where each line
  refers to a possible edge, and the columns are: vertex 1, vertex 2, value of
  the statistic, and number of estimated parameters (degrees of freedom) 
  associated with the edge.
}

\details{
  Calculates pairwise statistics (-2*log-LR, AIC, or BIC) for all possible
  edges, returning the values sorted in descending order.
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk}) \cr
Rodrigo Labouriau (\email{Rodrigo.Labouriau@agrsci.dk}) \cr
David Edwards (\email{David.Edwards@agrsci.dk})
}

\examples{
set.seed(7,kind="Mersenne-Twister")
dataset <- matrix(rnorm(1000),nrow=100,ncol=10)
m <- calcStat(dataset,stat="BIC")

data(dsCont)
# m1 <- calcStat(dataset,homog=TRUE,forbEdges=NULL,stat="LR")
#          1. in this case, there is no use for homog
#          2. no forbidden edges
#          3. the measure used is the LR (the result is a tree)
v <- calcStat(dsCont,homog=TRUE,forbEdges=NULL,stat="LR")

# result
head(v)
# column 1: first vertex of the edge
# column 2: second vertex of the edge
# column 3: in this case, -LR
# column 4: number of parameters for the edge
#         [,1] [,2]     [,3] [,4]
#    [1,]   17   27 393.0072    1
#    [2,]   21   27 343.5780    1
#    [3,]   22   25 306.0097    1
#    [4,]   17   21 302.9414    1
#    [5,]   27   32 300.0275    1
#    [6,]   21   32 289.4179    1
}
\keyword{graphs}
