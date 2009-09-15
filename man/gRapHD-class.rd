\name{gRapHD-class}
\docType{class}
\alias{gRapHD-class}
\alias{initialize.gRapHD}
\alias{matrix.gRapHD}
\alias{summary.gRapHD}
\alias{print.gRapHD}
\alias{show.gRapHD}
\alias{gRapHD.graphNEL}
\alias{graphNEL.gRapHD}

\title{Class "gRapHD"}

\description{S4 class.}

\section{Objects from the Class}{
Objects can be created, for example, by \code{new("gRapHD", ...)}.
Where \code{...} may contain the initial attribute names and values, as
described in the Slots section.
}

\section{Slots}{
  \describe{
    \item{\code{edges}:}{integer, matrix with 2 columns, each row representing
                         one edges and each column one of the vertices in the
                         edge. Column 1 containsthe vertex with lower index.}
    \item{\code{homog}:}{logical, TRUE if the covariance is homogeneous.}
    \item{\code{minForest}:}{integer, first and last edges found with minForest.}
    \item{\code{numCat}:}{integer, vector with number of levels for each variable (0 if continuous).}
    \item{\code{numP}:}{integer, vector with number of estimated parameters for each edge.}
    \item{\code{p}:}{integer, number of variables (vertices) in the model.}
    \item{\code{stat.minForest}:}{character, measure used (LR, AIC, or BIC).}
    \item{\code{stat.stepw}:}{character, measure used (LR, AIC, or BIC).}
    \item{\code{stat.user}:}{character, user defined.}
    \item{\code{statSeq}:}{numeric, vector with value of stat.minForest for each edge.}
    \item{\code{stepw}:}{interger, first and last edges found with stepw.}
    \item{\code{userDef}:}{integer, first and last edges defined by the user.}
    \item{\code{vertNames}:}{character, vector with vertices' names.}
  }
}


\section{Methods}{
  \describe{
    \item{\code{setMethod("initialize","gRapHD",initialize.gRapHD)}}{}
    \item{\code{setMethod("summary",signature(object="gRapHD"),summary.gRapHD)}}{}
    \item{\code{setMethod("print",signature(x="gRapHD"),print.gRapHD)}}{}
    \item{\code{setMethod("show",signature(object="gRapHD"),show.gRapHD)}}{}
    \item{\code{setMethod("plot",signature(x="gRapHD"),plot.gRapHD)}}{}
    \item{\code{setAs(from="matrix",to="gRapHD",def=matrix.gRapHD)}}{}
    \item{\code{setAs(from="gRapHD",to="graphNEL",def=gRapHD.graphNEL)}}{}
    \item{\code{setAs(from="graphNEL",to="gRapHD",def=graphNEL.gRapHD)}}{}
  }
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
Gabriel Coelho Goncalves de Abreu (\email{abreu_ga@yahoo.com.br}) \cr
David Edwards (\email{David.Edwards@agrsci.dk})
}

\references{
R. Gentleman, Elizabeth Whalen, W. Huber and S. Falcon. graph: A package to
handle graph data structures. R package version 1.22.2.
}

\examples{
# convertion from gRapHD to graphNEL
edges <- matrix(c(1,2,1,3,1,4),,2,byrow=TRUE)
g <- as(edges,"gRapHD")
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
g1 <- as(g,"graphNEL")
# A graphNEL graph with undirected edges
# Number of Nodes = 4
# Number of Edges = 3
g1@nodes # the nodes names


as(matrix(integer(0),,2),"gRapHD")
# note that the vertices must be numbered consecutively from 1. In the
# following, vertex 2 is added as an isolated vertex.
m1 <- as(matrix(c(1,3,1,4),,2,byrow=TRUE),"gRapHD")
\dontrun{plot(m1)}

}

\keyword{classes}
