\name{DFS}
\alias{DFS}
\title{Depth-first search}
\description{
  Returns all vertices reachable from one specific vertex (assuming that there
  are no cycles).
}

\usage{
DFS(model=NULL,edges=NULL, v, p=NULL)
}

\arguments{
  \item{model}{\code{gRapHD} object.}
  \item{edges}{matrix with 2 columns, each row representing one edge, and each
               column one of the vertices in the edge.}
  \item{v}{initial vertex (\code{0<v<=p}).}
  \item{p}{number of vertices.}
}

\value{
  Vector with all vertices reachable from \code{v}, or \code{0} if \code{v} is
  an isolated vertex.
}

\details{
  Given a list of edges, and a specific vertex \code{v}, returns a vector with all 
  vertices in the connected component containing v. The function assumes that the graph is acyclic.
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{abreu_ga@yahoo.com.br})
}

\references{
Cormen, T.H., Leiserson, C.E., Rivest, R.L. and Stein, C.
\emph{Introduction to Algorithms}, 2nd Edition. MIT Press and McGraw-Hill, 2001,
pp.540:9.
}

\examples{
set.seed(7,kind="Mersenne-Twister")
dataset <- matrix(rnorm(1000),nrow=100,ncol=10)
m <- minForest(dataset,stat="BIC")

DFS(edges=m@edges,v=1,p=10)
# [1] 5 2 9 8
#######################################################################
data(dsDiscr)
m1 <- minForest(dsDiscr,homog=TRUE,forbEdges=NULL,stat="BIC")
vertices <- DFS(edges=m1@edges, v=1, p=m1@p)

# result
vertices
# numeric(0)
# meaning that 1 is an isolated vertex

# OR
m1 <- minForest(dsDiscr,homog=TRUE,forbEdges=NULL,stat="LR")
vertices <- DFS(edges=m1@edges, v=1, p=m1@p)

# result
vertices
# [1]  4  8 12 19 18 14  7 17  5  3 10 13 15  9  6 20 16 11  2
# meaning that 1 reachs all vertices (a tree)
}
\keyword{graphs}
