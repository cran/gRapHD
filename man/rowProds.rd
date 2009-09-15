\name{rowProds}
\alias{rowProds}
\title{Row products}
\description{
  Form row products for numeric arrays.
}

\usage{
rowProds(x,na.rm=TRUE)
}

\arguments{
  \item{x}{matrix.}
  \item{na.rm}{logical. Whether missing values (including \code{NaN}) are
               omitted from the calculations.}
}

\value{
  Vector with length \code{nrow(x)}.
}

\details{
  Equivalent to use of apply with FUN = prod and MARGIN = 1, but is faster.
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk})
}

\examples{
set.seed(1,kind="Mersenne-Twister")
a <- matrix(rnorm(100),nrow=5)
x <- rowProds(x=a, na.rm=TRUE)
# x
# [1] -3.359208e-07 -2.861043e-10 -2.831108e-08 
# [4] -5.451996e-07  3.057436e-04
}
\keyword{array}
