\name{is.gRapHD}
\alias{is.gRapHD}
\title{Tests is a object is of "gRapHD" class.}
\description{
  \code{\link{Methods}} for class \code{gRapHD}.
}

\usage{
  \method{is}{gRapHD}(object)
}

\arguments{
  \item{object}{an \code{R} object.}
}

\details{
A logical constant, \code{TRUE} if the argument is of class \code{gRapHD}.}

\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk}) \cr
Rodrigo Labouriau (\email{Rodrigo.Labouriau@agrsci.dk}) \cr
David Edwards (\email{David.Edwards@agrsci.dk})
}

\examples{
  data(dsCont)
  m1 <- minForest(dsCont,homog=TRUE,forbEdges=NULL,stat="LR")
  is(m1,"gRapHD")
  is.gRapHD(m1)
}
\keyword{graphs}
