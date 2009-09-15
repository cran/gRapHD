\name{ccoeff}
\alias{ccoeff}
\title{Clustering coefficient}
\description{
  Returns the clustering coefficients of the vertices in a graph.
}

\usage{
ccoeff(model=NULL,edges=NULL,p=NULL)
}

\arguments{
  \item{model}{\code{gRapHD} object.}
  \item{edges}{matrix with 2 columns, each row representing one edge, and each
               column one of the vertices in the edge. Column 1 contains the
               vertex with lower index.}
  \item{p}{number of vertices. If \code{NULL}, the \code{p=max(edges)}.}
}

\value{
  A vector with length \code{p} with the clustering coefficient of each vertex.
}

\details{
The clustering coefficient is given by \code{C_i=2*e_i/(k_i*(k_i-1))}, where
\code{k_i} is the number of neighbours the vertex \code{i} has, and \code{e_i}
is the number of edges between the neighbours of \code{i}.
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk})
}

\examples{
data(dsCont)
m <- minForest(dsCont,homog=TRUE,forbEdges=NULL,stat="BIC")
m1 <- stepw(m,dsCont)
cc <- ccoeff(edges=m1$edges,p=m1$p)
mean(cc)
}
\keyword{graphs}
