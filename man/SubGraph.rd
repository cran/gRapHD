\name{SubGraph}
\alias{SubGraph}
\title{Generates a subgraph}
\description{
  Generates a sub-graph.
}

\usage{
SubGraph(model=NULL,edges=NULL,v=NULL,p=0)
}

\arguments{
  \item{model}{\code{gRapHD} object.}
  \item{edges}{matrix with 2 columns, each row representing one edge, and each
               column one of the vertices in the edge.}
  \item{v}{list of vertices in the sub-graph.}
  \item{p}{Number of vartices (used only if \code{edges} is not \code{NULL}).}
}

\value{
  Returns a \code{gRapHD} object, in which the edge list contains only
  edges where both vertices are in \code{v}.
}

\details{
  Based on a list of vertices, generate a sub-graph.
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk})
}

\examples{
data(dsCont)
m1 <- minForest(dsCont,homog=TRUE,forbEdges=NULL,stat="LR")
plot(m1,numIter=1000)

v <- c(1,11,21,19,30,25,22,24,34,9,20,29)
subM1 <- SubGraph(model=m1,v=v)
plot(subM1,numIter=1500)
}
\keyword{graphs}
