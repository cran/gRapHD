\name{shortPath}
\alias{shortPath}
\title{Shortest path}
\description{
  Shortest paths between one vertex and all other vertices.
}

\usage{
shortPath(model=NULL,edges=NULL,v=NULL,p=NULL)
}

\arguments{
  \item{model}{\code{gRapHD} object.}
  \item{edges}{matrix with 2 columns, each row representing one edge, and each
               column one of the vertices in the edge. Column 1 contains the
               vertex with lower index.}
  \item{v}{vertex.}
  \item{p}{number of vertices. If \code{NULL}, the number of vertices will
           be considered as equal the maximum vertices' index, i.e.,
           \code{p=max(edges)}.}
}

\value{
  \item{vector}{with length equal \code{p}, with the shortest path length from
                \code{v} to each other vertex.}
}

\details{
Calculates the shortest path between the vertex \code{v} and all other
vertices.}

\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk}) \cr
Rodrigo Labouriau (\email{Rodrigo.Labouriau@agrsci.dk}) 
}

\examples{
data(dsCont)
m1 <- minForest(dsCont,homog=TRUE,forbEdges=NULL,stat="LR")
shortPath(edges=m1$edges,v=1)
}
\keyword{graphs}
