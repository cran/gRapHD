\name{convData}
\alias{convData}
\title{Convert dataset}
\description{
  Convert a dataset to a structure required by other functions.
}

\usage{
  convData(dataset)
}

\arguments{
  \item{dataset}{matrix or data frame (discrete varibles must be factors).}
}

\value{
  List:
  \item{dataset}{Matrix.}
  \item{varType}{Vector: 0 for continuous, 1 for discrete.}
  \item{numCat}{Vector with the number of levels in each variable (0 if
                continuous).}
}

\details{
  Convert the dataset to the structure required by the other functions.
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk}) \cr
Rodrigo Labouriau (\email{Rodrigo.Labouriau@agrsci.dk}) \cr
David Edwards (\email{David.Edwards@agrsci.dk})
}


\examples{
  data(dsDiscr)
  ds <- convData(dsDiscr)
}
\keyword{graphs}
