\name{findEd}
\alias{findEd}
\title{Internal use}
\description{
  Find all add-eligible edges for a given triangulated graph.
}

\usage{
  findEd(edges,p,previous=NULL,varType,from=0,exact=FALSE,join=FALSE)
}

\arguments{
  \item{edges}{matrix with 2 columns, each row representing one edge, and each
               column one of the vertices in the edge.}
  \item{p}{number of vertices.}
  \item{previous}{result of a previous run of \code{findEd}.}
  \item{varType}{vector indicating the type of each variable: 0 if continuous,
                 or 1 if discrete.}
  \item{from}{initial vertex to be used in \code{\link[gRapHD:MCS]{MCS}}.}
  \item{exact}{logical indicating if the exact algorithm for finding
               add-eligible edges is to be used. Default is
               \code{FALSE}.}
  \item{join}{logical indicating if the disjoint components can be joined.
              Default is \code{FALSE}.}
}

\value{
  A list with:
  \item{edges.to.test}{matrix (k by 5), with columns:\cr
                      1 - first vertex of the tested edge\cr
                      2 - second vertex (the two values are ordered)\cr
                      3 - index for the separator in S\cr
                      4 - change in the LR by adding this edge\cr
                      5 - number of parameters for that edge}
  \item{S}{list with the separators.}
}

\details{
  Returns all add-eligible edges for a given triangulated graph, that is, edges that preserve the triangulated
  property when added. In the case of a mixed graph, only edges that do not result in forbidden paths are returned.\cr
  The argument \code{from} can be used to indicate the initial vertex used in
  the MCS algorithm. If \code{0}, the first vertex is used.\cr
  If \code{exact} is \code{FALSE}, the edge list may contain a few
  edges that are not add-eligible. Further tests (for
  example \code{mcs}) will be required before adding edges. Otherwise, the list contains
  only edges that preserve triangularity. That is, each edge that may be added
  to the graph such that the resulting graph is triangulated.\cr
  For graphs with both discrete and continuous vertices, the graph should be
  triangulated and contain no forbidden paths, and the edges that may be added
  preserving both properties are returned. See Lauritzen (1996), p. 11-13.
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk}) \cr
Rodrigo Labouriau (\email{Rodrigo.Labouriau@agrsci.dk}) \cr
David Edwards (\email{David.Edwards@agrsci.dk})
}

\references{
Lauritzen, S.L. (1996) \emph{Graphical Models}, Oxford University Press.\cr
}

\examples{
  edges <- matrix(c(1,2,2,3,2,4,2,5,2,6,3,4,4,5,5,6),ncol=2,byrow=TRUE)
  addEligible <- findEd(edges=edges,p=6,previous=NULL,varType=rep(0,6),from=1)
  #    > str(addEligible)
  #    List of 2
  #     $ edges: num [1:7, 1:5] 1 1 3 1 4 3 1 3 4 5 ...
  #     $ S    :List of 6
  #      ..$ : int 2
  #      ..$ : int [1:2] 2 4
  #      ..$ : int 2
  #      ..$ : int [1:2] 2 5
  #      ..$ : int 2
  #      ..$ : int 2
  #    > addEligible$edges
  #         [,1] [,2] [,3] [,4] [,5]
  #    [1,]    1    3    1    0    0
  #    [2,]    1    4    1    0    0
  #    [3,]    3    5    2    0    0
  #    [4,]    1    5    3    0    0
  #    [5,]    4    6    4    0    0
  #    [6,]    3    6    5    0    0
  #    [7,]    1    6    6    0    0

  # the columns in addEligible$edges
  #    1: first vertex in the edge
  #    2: second vertex in the edge
  #    3: index os the separator in addEligible$S
  #    4: change in the LRT for the edge (used if previous != NULL)
  #    5: number of parameters for the edge (used if previous != NULL)

  # note that the edge 3-6 (row 6) is actually a "false positive". If it's used
  # from=3,4,5, or 6, this "error" doesn't happen.
}
\keyword{graphs}
