\name{Degree}
\alias{Degree}
\title{Degree}
\description{
  Returns the degree of a set of vertices.
}

\usage{
Degree(model=NULL,edges=NULL,v=NULL)
}

\arguments{
  \item{model}{\code{gRapHD} object.}
  \item{edges}{matrix with 2 columns, each row representing one edge, and each
               column one of the vertices in the edge. Column 1 contains the
               vertex with lower index.}
  \item{v}{set of vertices.}
}

\value{
  \item{vector}{\code{length(v)}.}
}

\details{
Calculates the degree of each vertex in \code{v}. If \code{v=NULL}, it returns
the degree of all vertices in \code{edges}. If \code{v} contains a vertex not in
edges, the corresponding value is \code{NA}.}

\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk}) \cr
Rodrigo Labouriau (\email{Rodrigo.Labouriau@agrsci.dk}) 
}

\examples{
data(dsCont)
m1 <- minForest(dsCont,homog=TRUE,forbEdges=NULL,stat="LR")
Degree(model=m1)
}
\keyword{graphs}
