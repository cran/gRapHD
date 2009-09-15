\name{neighbours}
\alias{neighbours}
\title{Finds all direct neighbours}
\description{
  Finds all direct neighbours of a given vertex in a given graph.
}

\usage{
neighbours(model=NULL,edges=NULL,v)
}

\arguments{
  \item{model}{\code{gRapHD} object.}
  \item{edges}{matrix with 2 columns, each row representing one edge, and each
               column one of the vertices in the edge.}
  \item{v}{reference vertex.}
}

\value{
  Vector with all neighbours of vertex \code{v}.
}

\details{
  Returns all vertices with a direct connection with vertex \code{v} in
  \code{edges}.
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{abreu_ga@yahoo.com.br})
}

\examples{
data(dsCont)
m1 <- minForest(dsCont,homog=TRUE,forbEdges=NULL,stat="LR")
neigh <- neighbours(edges=m1@edges, v=22)
#    > neigh
#    [1]  3  9 24 25
}
\keyword{graphs}
