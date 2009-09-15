\name{stepw}
\alias{stepw}
\title{Stepwise forward selection}
\description{
  Stepwise forward selection.
}

\usage{
  stepw(model,dataset,stat="BIC",saveCH=NULL,forbEdges=NULL,exact=FALSE,
        initial=NULL,threshold=0,join=FALSE)
}

\arguments{
  \item{model}{\code{gRapHD} object.}
  \item{dataset}{matrix (\code{nrow(dataset)} observations and
                 \code{ncol(dataset)} variables).}
  \item{stat}{measure to be minimized: LR, AIC, or BIC. Default is BIC.
              It can also be a user defined function with
              format: \code{FUN(model,dataset,previous,forbEdges)}; where the
              parameters are defined as in
              \code{\link[gRapHD:chStat]{chStat}}.
              The function must return a structure as in
              \code{\link[gRapHD:chStat]{chStat}}.}
  \item{saveCH}{pattern of a file name to save each iteration results.
                \code{NULL} not to save (default).}
  \item{forbEdges}{list with edges that should not be considered. Matrix with 2
                   columns, each row representing one edge, and each column a
                   vertex. Default is \code{NULL}.}
  \item{exact}{logical indicting whether the exact algorithm for finding add-eligible
               edges is to be used. Default is \code{FALSE}.}
  \item{initial}{continue the algorithm from a previous point. The parameter
                 must have the same structure as the result of
                 \code{\link[gRapHD:chStat]{chStat}}.}
  \item{threshold}{values greater than the threshold are not included in the
                   model. Default is \code{0}.}
  \item{join}{logical indicating whether disjoint components can be joined.
              Default is \code{FALSE}.}
}

\value{
  The same structure as \code{model}, but adding the new edges, and updating
  other relevant information. See \code{\link[gRapHD:minForest]{minForest}} for
  the description of the structure.
}

\details{
  Performs a stepwise forward selection of edges to be added to a graph. Only
  edges preserving the triangularity are considered
  (\code{\link[gRapHD:findEd]{findEd}}). At each step the edge giving the greatest improvement in the
  chosen statistic is added. The process ends when no further improvement is possible.\cr
  If \code{exact} is \code{FALSE}, the list of edges returned may contain
  a few edges that do not preserve triangularity, requiring further tests (for
  example \code{mcs}) before adding an edge. Otherwise, the list will only contain
  edges that preserve triangularity. That is, each edge that may be added
  to the graph such that the resulting graph is triangulated.\cr
  For graphs with both discrete and continuous vertices, the graph should be
  triangulated and contain no forbidden paths, and the edges that may be added
  must preserve both properties. See Lauritzen (1996), p. 11-13.
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk}) \cr
Rodrigo Labouriau (\email{Rodrigo.Labouriau@agrsci.dk}) \cr
David Edwards (\email{David.Edwards@agrsci.dk})
}

\examples{
  set.seed(7,kind="Mersenne-Twister")
  dataset <- matrix(rnorm(1000),nrow=100,ncol=10)
  m <- minForest(dataset,stat="BIC")

  M <- stepw(m,dataset,stat="LR",NULL,NULL)

  ##############################################################################
  # Example with continuous variables
  data(dsCont)
  # m1 <- minForest(dataset,homog=TRUE,forbEdges=NULL,stat="LR")
  #          1. in this case, there is no use for homog
  #          2. no forbidden edges
  #          3. the measure used is the LR (the result is a tree)
  m1 <- minForest(dsCont,homog=TRUE,forbEdges=NULL,stat="LR")
  plot(m1,numIter=1000)

  # m2 <- stepw(m1,dataset,stat="BIC",saveCH=NULL,forbEdges=NULL)
  #          1. m1 is the result of minForest
  #          2. the same dataset (this is not checked)
  #          3. the default is BIC
  #          4. if saveCH="XXX", a file "XXX_00000i.RData" is saved for each iter
  #          5. no forbidden edges
  m2 <- stepw(m1,dsCont,stat="BIC",saveCH=NULL,forbEdges=NULL)
  plot(m2,numIter=1000)

  ##############################################################################
  # Example with discrete variables
  data(dsDiscr)
  # m1 <- minForest(dataset,homog=TRUE,forbEdges=NULL,stat="LR")
  #          1. in this case, there is no use for homog
  #          2. no forbidden edges
  #          3. the measure used is the LR (the result is a tree)
  m1 <- minForest(dsDiscr,homog=TRUE,forbEdges=NULL,stat="LR")
  plot(m1,numIter=1000)

  # m2 <- stepw(m1,dataset,stat="BIC",saveCH=NULL,forbEdges=NULL)
  #          1. m1 is the result of minForest
  #          2. the same dataset (this is not checked)
  #          3. the default is BIC
  #          4. if saveCH="XXX", a file "XXX_00000i.RData" is saved for each iter
  #          5. no forbidden edges
  m2 <- stepw(m1,dsDiscr,stat="BIC",saveCH=NULL,forbEdges=NULL)
  plot(m2,numIter=1000)

  ##############################################################################
  # Example with mixed variables
  data(dsMixed)
  # m1 <- minForest(dataset,homog=TRUE,forbEdges=NULL,stat="LR")
  #          1. it is to be considered homogeneous
  #          2. no forbidden edges
  #          3. the measure used is the LR (the result is a tree)
  m1 <- minForest(dsMixed,homog=TRUE,forbEdges=NULL,stat="LR")
  plot(m1,numIter=1000)

  # m2 <- stepw(m1,dataset,stat="BIC",saveCH=NULL,forbEdges=NULL)
  #          1. m1 is the result of minForest
  #          2. the same dataset (this is not checked)
  #          3. the default is BIC
  #          4. if saveCH="XXX", a file "XXX_00000i.RData" is saved for each iter
  #          5. no forbidden edges
  m2 <- stepw(m1,dsMixed,stat="BIC",saveCH=NULL,forbEdges=NULL)
  plot(m2,numIter=1000)

  ##############################################################################
  # Example using a user defined function
  #   The function userFun calculates the same edges weigths as the option
  # stat="BIC". It means that the final result, using either option, is the
  # same.
  userFun <- function(model,dataset,previous=NULL,forbEdges=NULL)
  {
    p <- ncol(dataset) # number of variables (vertices)
    n <- nrow(dataset) # number of observations
    edges.to.test <- previous$edges.to.test
    SS <- previous$S # minimal separators
    rm(previous)

    num <- nrow(edges.to.test)

    if (num > 0)
      for (i in 1:num)
      {
        x <- edges.to.test[i,1]
        y <- edges.to.test[i,2]
        if ((edges.to.test[i,4]==0) &
            (!is.element((x-1)*p-(x-1)*x/2+y-x,forbEdges)))
        {
          S <- SS[[edges.to.test[i,3]]]
          clique <- c(edges.to.test[i,1:2],S)
          CM <- cov(dataset[,clique],use="pairwise.complete.obs")*(nrow(dataset)-1)
          a <- match(edges.to.test[i,1:2],clique)
          d <- log(det(CM)) + log(det(matrix(CM[-a,-a],length(clique)-2)))
          d <- d - (log(det(matrix(CM[-a[1],-a[1]],length(clique)-1))) +
                    log(det(matrix(CM[-a[2],-a[2]],length(clique)-1))))
          edges.to.test[i,4] <- nrow(dataset)*d + log(n)
          edges.to.test[i,5] <- 1
        }
      }
    return(list(edges.to.test=edges.to.test,S=SS))
  }

  data(dsCont)
  m <- minForest(dsCont,stat="BIC")
  m1 <- stepw(m,dsCont,stat="BIC")
  m2 <- stepw(m,dsCont,stat=userFun)
  identical(m1$edges,m2$edges)

  ##############################################################################
  # Example of joining disjoint components
  data(dsCont)
  m <- minForest(dsCont,stat="BIC")
  m1 <- stepw(m,dsCont)
  m2 <- m
  m2$edges <- matrix(integer(0),,2)
  m3 <- stepw(m2,dsCont,join=TRUE)
  identical(sortMat(m1$edges),sortMat(m3$edges))
  # TRUE

}
\keyword{graphs}
