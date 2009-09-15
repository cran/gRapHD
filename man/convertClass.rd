\name{convertClass}
\alias{convertClass}
\title{Convert object between classes}
\description{
  Convert objects between gRapHD and graphNEL classes.
}

\usage{
  convertClass(object)
}

\arguments{
  \item{object}{an object of \code{gRapHD} or \code{graphNEL} class.}
}

\value{
  An object with the opposite class to the input object.
}

\details{
  As the gRapHD class does not use variables names, but only the column
  indexes in the dataset, the conversion may change the variables references. When
  converting a \code{gRapHD} object to a \code{graphNEL} object, the nodes
  names in the new object are only the original indexes converted to character.
  When doing the reverse conversion, the nodes indexes in the new
  \code{gRapHD} object are the respective indexes in the original element
  \code{nodes} in the \code{graphNEL} object. See the example below.
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk}) \cr
Rodrigo Labouriau (\email{Rodrigo.Labouriau@agrsci.dk}) \cr
David Edwards (\email{David.Edwards@agrsci.dk})
}

\references{
R. Gentleman, Elizabeth Whalen, W. Huber and S. Falcon. graph: A package to
handle graph data structures. R package version 1.22.2.
}

\examples{
  # convertion from gRapHD to graphNEL
  edges <- matrix(c(1,2,1,3,1,4),,2,byrow=TRUE)
  g <- as.gRapHD(edges)
  #List of 9
  # $ edges    : num [1:3, 1:2] 1 1 1 2 3 4
  # $ p        : int 4
  # $ stat.user: chr "LR"
  # $ statSeq  : num [1:3] NA NA NA
  # $ varType  : int [1:4] 0 0 0 0
  # $ numCat   : int [1:4] 0 0 0 0
  # $ homog    : logi TRUE
  # $ numP     : num [1:3] NA NA NA
  # $ userDef  : num [1:2] 1 3
  # - attr(*, "class")= chr "gRapHD"
  g1 <- convertClass(g)
  # A graphNEL graph with undirected edges
  # Number of Nodes = 4
  # Number of Edges = 3
  g1@nodes # the nodes names
}
\keyword{graphs}
