\name{adjMat}
\alias{adjMat}
\title{Adjacency matrix.}
\description{
  Returns the adjacency matrix based on a list of edges, supplied in a \code{gRapHD} object or as a matrix.
}

\usage{
  adjMat(model=NULL,edges=NULL,p=NULL)
}

\arguments{
  \item{model}{\code{gRapHD} object.}
  \item{edges}{matrix with 2 columns, each row representing one edge, and each
               column one of the vertices in the edge. Column 1 contains the
               vertex with lower index.}
  \item{p}{number of vertices.}
}

\value{
  \item{matrix}{\code{p} by \code{p}.}
}

\details{
The dimension of the matrix is given by \code{model$p} or by the maximum of the number of 
vertices in \code{edges} and \code{p}.
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk}) \cr
Rodrigo Labouriau (\email{Rodrigo.Labouriau@agrsci.dk}) \cr
David Edwards (\email{David.Edwards@agrsci.dk})
}

\examples{
  data(dsCont)
  m1 <- minForest(dsCont,homog=TRUE,forbEdges=NULL,stat="LR")
  edges <- SubGraph(edges=m1$edges,v=1:10)$edges
  adjMat(edges=edges,p=10)
}
\keyword{graphs}
