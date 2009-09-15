\name{chStat}
\alias{chStat}
\title{Internal use}
\description{
  Calculates the deviance associated with the addition of each add-eligible edge.
  Called from \code{\link[gRapHD:stepw]{stepw}}.
}

\usage{
chStat(model,dataset,previous=NULL,forbEdges=NULL)
}

\arguments{
  \item{model}{a \code{gRapHD} object.}
  \item{dataset}{matrix (\code{nrow(dataset)} observations and
                 \code{ncol(dataset)} variables).}
  \item{previous}{result of a previous run of \code{chStat}.}
  \item{forbEdges}{list with edges that should not be considered. Matrix with 2
                   columns, each row representing an edge, and each column a
                   vertex. Default is \code{NULL}.}
}

\value{
  A list with:
  \item{edges.to.test}{matrix (k by 5), with columns:\cr
                      1 - first vertex of the tested edge\cr
                      2 - second vertex (the two values are ordered)\cr
                      3 - index for the separator in S\cr
                      4 - deviance statistic associated with adding this edge\cr
                      5 - degrees of freedom associated with the deviance}
  \item{S}{list of separators.}
}

\details{
  The deviance and degrees of freedom associated with
  each add-eligible edge is returned. If previous results are specified these are reused as appropriate, and a 
  concatenated result list is returned. This function is used by \code{stepw}.
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk}) \cr
Rodrigo Labouriau (\email{Rodrigo.Labouriau@agrsci.dk}) \cr
David Edwards (\email{David.Edwards@agrsci.dk})
}

\examples{
data(dsCont)
m1 <- minForest(dsCont,homog=TRUE,forbEdges=NULL,stat="LR")
ch <- findEd(m1$edges,m1$p,NULL,m1$varType,0)
ch <- chStat(m1,dsCont,ch,forbEdges=0)
str(ch)
# List of 2
#  $ edges.to.test: num [1:53, 1:5] 1 11 19 11 17 2 2 10 10 2 ...
#  $ S            :List of 53
#   ..$ : int 11
#   ..$ : int 21
#   ..$ : int 21
#   ..$ : int 21
#   ..$ : int 27
#   ..$ : int 17
#    ...
head(ch$edges.to.test)
#      [,1] [,2] [,3]         [,4] [,5]
# [1,]    1   21    1  -0.61733689    1
# [2,]   11   19    2  -0.24637623    1
# [3,]   19   27    3  -0.47194908    1
# [4,]   11   27    4  -7.00259895    1
# [5,]   17   21    5 -11.09310305    1
# [6,]    2   27    6  -0.04690911    1

# the columns in ch$edges.to.test
#    1: first vertex in the edge
#    2: second vertex in the edge
#    3: index os the separator in ch$S
#    4: change in the LR for the edge
#    5: number of parameters for the edge
}
\keyword{graphs}
