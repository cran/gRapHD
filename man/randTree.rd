\name{randTree}
\alias{randTree}
\title{Random tree}
\description{
  Generate a random tree.
}

\usage{
  randTree(p,seed=1)
}

\arguments{
  \item{p}{number of vertices.}
  \item{seed}{seed (for the random generator, see
              \code{\link[base:Random]{set.seed}}).}
}

\value{
  A list containing:
  \item{edges}{matrix with 2 columns, each row representing one edge, and each
               column one of the vertices in the edge.}
  \item{seed}{seed.}
  \item{p}{number of vertices.}
}



\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk}) \cr
Rodrigo Labouriau (\email{Rodrigo.Labouriau@agrsci.dk}) \cr
David Edwards (\email{David.Edwards@agrsci.dk})
}

\references{
Rodionov, A.S., Choo, H. On Generating Random Network Structures: Trees.
\emph{Springer-Verlag Lecture Notes in Computer Science}, vol. 2658,
pp. 879-887, June 2003.
}

\examples{
  tree <- randTree(p=10, seed=1)
  plotG(edges=tree$edges)
  #    > tree
  #    $edges
  #          [,1] [,2]
  #     [1,]    3    4
  #     [2,]    3    5
  #     [3,]    2    3
  #     [4,]    4    6
  #     [5,]    6    9
  #     [6,]    1    6
  #     [7,]    6   10
  #     [8,]    6    7
  #     [9,]    5    8
  #
  #    $seed
  #    [1] 1
  #
  #    $p
  #    [1] 10
}
\keyword{graphs}
