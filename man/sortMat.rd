\name{sortMat}
\alias{sortMat}
\title{Sort matrix}
\description{
  Sorts the rows of a matrix by given columns.
}

\usage{
sortMat(mat,cols)
}

\arguments{
  \item{mat}{matrix.}
  \item{cols}{sequence os columns to sort by.}
}

\value{
  Matrix.
}

\details{
  It is just a interface to the function \code{\link[base:order]{order}}.
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{abreu_ga@yahoo.com.br})
}

\examples{
set.seed(1,kind="Mersenne-Twister")
a <- matrix(c(sample(2,6,TRUE),sample(3,6,TRUE),sample(2,6,TRUE)),
            nrow=6)
x <- sortMat(mat=a, cols=c(1:3))
a
#      [,1] [,2] [,3]
# [1,]    1    3    2
# [2,]    1    2    1
# [3,]    2    2    2
# [4,]    2    1    1
# [5,]    1    1    2
# [6,]    2    1    2
x
#      [,1] [,2] [,3]
# [1,]    1    1    2
# [2,]    1    2    1
# [3,]    1    3    2
# [4,]    2    1    1
# [5,]    2    1    2
# [6,]    2    2    2
}
\keyword{array}
