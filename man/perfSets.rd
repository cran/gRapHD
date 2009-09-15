\name{perfSets}
\alias{perfSets}
\title{Finds a perfect sequence}
\description{
  Finds a perfect sequence, returning the cliques, histories, residuals, and
  separators of a given triangulated graph.
}

\usage{
perfSets(model=NULL,edges=NULL,p=NULL,varType=0,from=0)
}

\arguments{
  \item{model}{\code{gRapHD} object.}
  \item{edges}{matrix with 2 columns, each row representing one edge, and each
               column one of the vertices in the edge.}
  \item{p}{number of vertices.}
  \item{varType}{vector indicating the type of each variable: 0 if continuous,
                 or 1 if discrete.}
  \item{from}{initial vertex to be used in \code{\link[gRapHD:MCS]{MCS}}.}
}

\value{
  A list containing:
  \item{cliques}{list.}
  \item{histories}{list.}
  \item{residuals}{list.}
  \item{separators}{list.}
}

\details{
  Based on the perfect numbering of \code{mcs}, returns the perfect sequnece.\cr
  The sequence is given by the cliques in the graph:
  \eqn{C_j=closure(\alpha_j)\cap\{\alpha_1,...,\alpha_j\}}{Cj=closure(\alpha(j)) intersection \{\alpha(1),...,\alpha(j)\}}, \eqn{j\ge1}{j>=1}.\cr
  The other sets are given by:\cr
  - Histories: \eqn{H_j=C_1\cup...\cup C_j}{Hj=C1 union ... union Cj} \cr
  - Residuals: \eqn{R_j=C_j\backslash H_{j-1}}{Rj=Cj \\ H(j-1)} \cr
  - Separators: \eqn{S_j=H_{j-1}\cap C_j}{Sj=H(j-1) intersection Cj}

}

\author{
Gabriel Coelho Goncalves de Abreu (\email{abreu_ga@yahoo.com.br})
}

\references{
Lauritzen, S.L. (1996) \emph{Graphical Models}, Oxford University Press.\cr
}

\examples{
edges <- matrix(c(1,2,2,3,2,4,2,5,2,6,3,4,4,5,5,6),ncol=2,byrow=TRUE)
setList <- perfSets(edges=edges, p=6, varType=0, from=1)
#    > str(setList)
#    List of 4
#     $ cliques   :List of 4
#      ..$ : int [1:2] 1 2
#      ..$ : int [1:3] 2 3 4
#      ..$ : int [1:3] 2 4 5
#      ..$ : int [1:3] 2 5 6
#     $ histories :List of 4
#      ..$ : int [1:2] 1 2
#      ..$ : int [1:4] 1 2 3 4
#      ..$ : int [1:5] 1 2 3 4 5
#      ..$ : int [1:6] 1 2 3 4 5 6
#     $ separators:List of 4
#      ..$ : NULL
#      ..$ : int 2
#      ..$ : int [1:2] 2 4
#      ..$ : int [1:2] 2 5
#     $ residuals :List of 4
#      ..$ : int [1:2] 1 2
#      ..$ : int [1:2] 3 4
#      ..$ : int 5
#      ..$ : int 6
}
\keyword{graphs}
