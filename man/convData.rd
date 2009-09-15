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
  \item{numCat}{Vector with the number of levels in each variable (0 if
                continuous).}
  \code{vertNames} - vector with the original vertices' names.\cr
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
