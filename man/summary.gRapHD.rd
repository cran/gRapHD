\name{summary.gRapHD}
\alias{summary.gRapHD}
\title{Summarizing model.}
\description{
  \code{\link{Methods}} for class \code{gRapHD}.
}

\usage{
  \method{summary}{gRapHD}(object,...)
}

\arguments{
  \item{object}{an object of class \code{"gRapHD"}, usually, a result of a
               call to \code{\link[gRapHD:minForest]{minForest}}.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\details{
Give details about the model structure.}

\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk}) \cr
Rodrigo Labouriau (\email{Rodrigo.Labouriau@agrsci.dk}) \cr
David Edwards (\email{David.Edwards@agrsci.dk})
}

\examples{
  data(dsCont)
  m1 <- minForest(dsCont,homog=TRUE,forbEdges=NULL,stat="LR")
  summary(m1)
}
\keyword{graphs}
