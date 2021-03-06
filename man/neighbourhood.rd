\name{neighbourhood}
\alias{neighbourhood}
\title{Neighbourhood of a vertex}
\description{
  Finds the set of vertices with up to a given distance from a given vertex.
}

\usage{
neighbourhood(model=NULL,edges=NULL,orig=NULL,rad=1)
}

\arguments{
  \item{model}{\code{gRapHD} object.}
  \item{edges}{matrix with 2 columns, each row representing one edge, and each
               column one of the vertices in the edge. Column 1 contains the
               vertex with lower index.}
  \item{orig}{central vertex.}
  \item{rad}{distance.}
}

\value{
  Returns a list with:
  \item{subEdges}{matrix with 2 columns, each row representing one edge, and each
               column one of the vertices in the edge. Column 1 contains the
               vertex with lower index.}
  \item{v}{matrix with 2 columns, the first indicating the vertex index, and
           the second the distance to the \code{orig}.}
}

\details{
  Finds the set of vertices with up to a given distance from a given vertex.
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{abreu_ga@yahoo.com.br})
}

\examples{
data(dsCont)
m1 <- minForest(dsCont,homog=TRUE,forbEdges=NULL,stat="LR")
aux <- neighbourhood(model=m1,orig=27,rad=2)
plot(new("gRapHD",edges=aux$edges,p=m1@p),vert=aux$v[,1])
}
\keyword{graphs}
