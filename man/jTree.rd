\name{jTree}
\alias{jTree}
\title{Junction tree}
\description{
  Finds a junction tree.
}

\usage{
jTree(model)
}

\arguments{
  \item{model}{object of \code{gRapHD} class.}
}

\value{
  A list with:
  \item{separators}{list with unique minimal separators.}
  \item{juncTree}{edges in the tree (each vertex is a clique in the list below).}
  \item{sepSubSetOfSep}{list in which each element gives all the separators
                        which contain this respective separator.}
  \item{indSepOrig}{index of the original separator (in the \code{MCS} result)
                    in the list above.}
  \item{cliques}{list with cliques.}
}

\details{
  Returns one possible junction tree. Note that each edge is associated to one
  separator in the list, and a separator may be contained in other(s)
  separator(s). To identify which separator is associated to each edge is enough
  to check \code{ind<-indSepOrig[which(indSepOrig!=1)]}. In this way, the edge
  \code{juncTree[i,]} is associated with separator \code{ind[i]}.
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{abreu_ga@yahoo.com.br})
}

\examples{
edges <- matrix(c(1,2,2,3,2,4,2,5,2,6,3,4,4,
                  5,5,6,7,8,7,9,8,9,8,10,9,10),ncol=2,byrow=TRUE)
m <- new("gRapHD",edges=edges)
jT <- jTree(m)
str(jT)
# List of 5
#  $ separators    :List of 5
#   ..$ : int(0)
#   ..$ : int 2
#   ..$ : int [1:2] 2 4
#   ..$ : int [1:2] 2 5
#   ..$ : int [1:2] 8 9
#  $ juncTree      : int [1:4, 1:2] 1 2 3 5 2 3 4 6
#  $ sepSubSetOfSep:List of 5
#   ..$ : int [1:4] 2 3 4 5
#   ..$ : int [1:2] 3 4
#   ..$ : int(0)
#   ..$ : int(0)
#   ..$ : int(0)
#  $ indSepOrig    : int [1:6] 1 2 3 4 1 5
#  $ cliques       :List of 6
#   ..$ : int [1:2] 1 2
#   ..$ : int [1:3] 2 3 4
#   ..$ : int [1:3] 2 4 5
#   ..$ : int [1:3] 2 5 6
#   ..$ : int [1:3] 7 8 9
#   ..$ : int [1:3] 8 9 10
}
\keyword{graphs}
