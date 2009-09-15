\name{MCS}
\alias{MCS}
\title{Maximum cardinality search}
\description{
  Returns a perfect ordering of the edges.
}

\usage{
  MCS(model=NULL,edges=NULL,v=0,p=NULL)
}

\arguments{
  \item{model}{\code{gRapHD} object.}
  \item{edges}{matrix with 2 columns, each row representing one edge, and each
               column one of the vertices in the edge.}
  \item{v}{initial vertex (\code{0<=v<=p}). If \code{v=0}, the algorithm
           starts from vertex 1.}
  \item{p}{number of vertices.}
}

\value{
  Zero if the graph is not triangulated, or a vector with the number of
  each vertex in the perfect ordering.
}

\details{
  Returns a perfect ordering of the vertices. For mixed graphs, the discrete vertices appear before the
  continuous vertices, in each connected component.
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk}) \cr
Rodrigo Labouriau (\email{Rodrigo.Labouriau@agrsci.dk}) \cr
David Edwards (\email{David.Edwards@agrsci.dk})
}

\references{
Tarjan, R.E., Yannakakis, M. Simple linear-time algorithms to test chordality of
graphs, test acyclicity of hypergraphs, and selectively reduce acyclic
hypergraphs. \emph{SIAM J. Comput.}, Vol 13, 3:566-79. \cr

Leimer, H. Triangulated graphs with marked vertices. \emph{Ann. Discr. Maths.}, Vol 41, 311-324.
}

\examples{
  set.seed(7,kind="Mersenne-Twister")
  dataset <- matrix(rnorm(1000),nrow=100,ncol=10)
  m <- minForest(dataset,stat="BIC")

  MCS(edges=m$edges,v=1,p=10)
}
\keyword{graphs}
