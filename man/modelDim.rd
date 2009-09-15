\name{modelDim}
\alias{modelDim}
\title{Calculate the dimension of a model.}
\description{
  Calculates the number of free parameters in the model.
}

\usage{
  modelDim(model)
}

\arguments{
  \item{model}{gRapHD class.}
}

\value{
  Number of free parameters in the model.
}

\details{
  See Lauritzen (1996), pages 202-203, and 215-216 for more details.
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk}) \cr
Rodrigo Labouriau (\email{Rodrigo.Labouriau@agrsci.dk}) \cr
David Edwards (\email{David.Edwards@agrsci.dk})
}

\references{
Lauritzen, S.L. (1996) \emph{Graphical Models}, Oxford University Press.\cr
}


\examples{
  data(dsMixed)
  m <- minForest(dsMixed,homog=TRUE,stat="LR")
  modelDim(m)
  # 71

  m <- minForest(dsMixed,homog=FALSE,stat="LR")
  modelDim(m)
  # 111
}
\keyword{graphs}
