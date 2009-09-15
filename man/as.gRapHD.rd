\name{as.gRapHD}
\alias{as.gRapHD}
\title{Coerces to an object of type "gRapHD"}
\description{
  \code{\link{Methods}} for class \code{gRapHD}.
}

\usage{
\method{as}{gRapHD}(object,...)
}

\arguments{
  \item{object}{list of edges: NULL, or integer with dimension \code{(k,2)}.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\details{
  \dots can be:\cr
  \code{p} - number of variables (vertices) in the model.\cr
  \code{stat} - measure used (LR, AIC, or BIC).\cr
  \code{statSeq} - vector with \code{-2(log-likelihood)} for each edge.\cr
  \code{vertNames} - vector with the original vertices' names.\cr
  \code{numCat} - vector with number of levels for each variable (0 if continuous).\cr
  \code{homog} - \code{TRUE} if the covariance is homogeneous.\cr
  \code{numP} - vector with number of estimated parameters for each edge.\cr
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk})
}

\examples{
as.gRapHD(NULL)
as.gRapHD(matrix(integer(0),,2))
as.gRapHD(matrix(integer(0),,2),p=10,stat="BIC",homog=FALSE)

# note that the vertices must be numbered consecutively from 1. In the
# following, vertex 2 is added as an isolated vertex.
m1 <- as.gRapHD(matrix(c(1,3,1,4),,2,byrow=TRUE))
\dontrun{plot(m1)}
}
\keyword{graphs}
