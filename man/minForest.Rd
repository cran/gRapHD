\name{minForest}
\alias{minForest}
\title{Minimum forest}
\description{
  Returns the forest that minimises the -2*log-likelihood, AIC, or BIC.
}

\usage{
  minForest(dataset,homog=TRUE,forbEdges=NULL,stat="BIC")
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
              Default is BIC. It can also be a user
              defined function with format: \code{FUN(newEdge,varType,numCat,
              dataset)}; where the parameters \code{varType} and \code{numCat}
              are as defined in the Value section; \code{newEdge} is a vector
              with length two; and \code{dataset} is a matrix (n by p).}
}

\value{
  A list containing:
  \item{edges}{matrix with 2 columns, each row representing one edge, and each
               column one of the vertices in the edge. Column 1 contains the
               vertex with lower index.}
  \item{p}{number of variables (vertices) in the model.}
  \item{stat.minForest}{measure used (LR, AIC, or BIC).}
  \item{statSeq}{vector with value of \code{stat.minForest} for each edge.}
  \item{varType}{vector indicating the type of each variable: 0 if continuous,
                 or 1 if discrete.}
  \item{numCat}{vector with number of levels for each variable (0 if continuous).}
  \item{homog}{\code{TRUE} if the covariance is homogeneous.}
  \item{numP}{vector with number of estimated parameters for each edge.}
  \item{minForest}{first and last edges found with \code{minForest}.}
}

\details{
Returns for the tree or forest that minimizes the -2*log-likelihood, AIC, or
BIC. If the log-likelihood is used, the result is a tree, if AIC or BIC is used,
the result is a tree or forest.The \code{dataset} contains variables
(vertices) in the columns, and observations in the rows. The result has vertices
numbered according to the column indexes in \code{dataset}.\cr
All discrete variables must be factors. All factor levels must be represented in
the data. Missing values are not allowed.
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk}) \cr
Rodrigo Labouriau (\email{Rodrigo.Labouriau@agrsci.dk}) \cr
David Edwards (\email{David.Edwards@agrsci.dk})
}

\references{
Chow, C.K. and Liu, C.N. (1968) Approximating discrete probability distributions
with dependence trees. \emph{IEEE Transactions on Information Theory},
Vol. IT-14, 3:462-7. \cr
Edwards, D., de Abreu, G.C.G. and Labouriau, R. (2009). High-dimensional Mixed
Graphical Models Using Minimal AIC and BIC forests. \emph{BMC Bioinformatics}.
(submitted).
}

\examples{
  set.seed(7,kind="Mersenne-Twister")
  dataset <- matrix(rnorm(1000),nrow=100,ncol=10)
  m <- minForest(dataset,stat="BIC")

  ##############################################################################
  # Example with continuous variables
  data(dsCont)
  # m1 <- minForest(dataset,varType=0,homog=TRUE,forbEdges=NULL,stat="LR")
  #          1. in this case, there is no use for homog
  #          2. no forbidden edges
  #          3. the measure used is the LR (the result is a tree)
  m1 <- minForest(dsCont,homog=TRUE,forbEdges=NULL,stat="LR")
  plotG(model=m1,numIter=1000)

  ##############################################################################
  # Example with discrete variables
  data(dsDiscr)
  # m1 <- minForest(dataset,varType=1,homog=TRUE,forbEdges=NULL,stat="LR")
  #          1. in this case, there is no use for homog
  #          2. no forbidden edges
  #          3. the measure used is the LR (the result is a tree)
  m1 <- minForest(dsDiscr,homog=TRUE,forbEdges=NULL,stat="LR")
  plotG(model=m1,numIter=1000)

  ##############################################################################
  # Example with mixed variables
  data(dsMixed)
  # m1 <- minForest(dataset,varType=1,homog=TRUE,forbEdges=NULL,stat="LR")
  #          1. it is to be considered homogeneous
  #          2. no forbidden edges
  #          3. the measure used is the LR (the result is a tree)
  m1 <- minForest(dsMixed,homog=TRUE,forbEdges=NULL,stat="LR")
  plotG(model=m1,numIter=1000)
  
  ##############################################################################
  # Example using a user defined function
  #   The function userFun calculates the same edges weigths as the option
  # stat="LR". It means that the final result, using either option, is the
  # same.
  userFun <- function(newEdge,varType,numCat,dataset)
  {
    sigma <- var(dataset[,newEdge])
    v <- nrow(dataset)*log(prod(diag(sigma))/det(sigma))
    return(c(v,1))
  }
  
  data(dsCont)
  m <- minForest(dsCont,stat="LR")
  m1 <- minForest(dsCont,stat=userFun)
  identical(m$edges,m1$edges)

}
\keyword{graphs}
