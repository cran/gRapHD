\name{as.gRapHD}
\alias{as.gRapHD}
\title{Coerce to an object of type "gRapHD".}
\description{
  \code{\link{Methods}} for class \code{gRapHD}.
}

\usage{
  \method{as}{gRapHD}(object,...)
}

\arguments{
  \item{object}{list of edges: NULL, or integer with dimension (k,2).}
  \item{\dots}{further arguments passed to or from other methods.}
}

\details{
  \dots can have:\cr
  \code{p} - number of variables (vertices) in the model.\cr
  \code{stat} - measure used (LR, AIC, or BIC).\cr
  \code{statSeq} - vector with \code{-2(log-likelihood)} for each edge.\cr
  \code{vertNames} - vector with the original vertices' names.\cr
  \code{numCat} - vector with number of levels for each variable (0 if continuous).\cr
  \code{homog} - \code{TRUE} if the covariance is homogeneous.\cr
  \code{numP} - vector with number of estimated parameters for each edge.\cr
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk}) \cr
Rodrigo Labouriau (\email{Rodrigo.Labouriau@agrsci.dk}) \cr
David Edwards (\email{David.Edwards@agrsci.dk})
}

\examples{
  as.gRapHD(NULL)
  as.gRapHD(matrix(integer(0),,2))
  as.gRapHD(matrix(integer(0),,2),p=10,stat="BIC",homog=FALSE)
}
\keyword{graphs}
