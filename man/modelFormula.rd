\name{modelFormula}
\alias{modelFormula}
\title{Model's formula}
\description{
  Returns the formula of the model.
}

\usage{
modelFormula(model)
}

\arguments{
  \item{model}{gRapHD class.}
}

\value{
  List with the generators of the model:
  \item{discrete}{terms \eqn{(d,\emptyset)}}
  \item{linear}{terms \eqn{(d,\gamma^2)}}
  \item{quadratic}{terms \eqn{(d,\gamma)} and \eqn{(d,\{\gamma,\mu\})}}
  \item{quadratic2}{terms \eqn{(d,c^2)}}
}

\details{
  See Lauritzen (1996), pages 202-203, and 215-216 for more details.
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{abreu_ga@yahoo.com.br})
}

\references{
Lauritzen, S.L. (1996) \emph{Graphical Models}, Oxford University Press.\cr
}

\examples{
data(dsMixed)
m <- minForest(dsMixed,homog=TRUE,stat="LR")
str(modelFormula(m))
#List of 4
# $ discrete  :List of 4
#  ..$ : int [1:2] 1 3
#  ..$ : int [1:2] 3 4
#  ..$ : int [1:2] 3 5
#  ..$ : int [1:2] 5 2
# $ linear    : list()
# $ quadratic :List of 5
#  ..$ : num [1:2] 5 8
#  ..$ : num [1:2] 5 9
#  ..$ : num [1:2] 4 11
#  ..$ : num [1:2] 5 14
#  ..$ : num [1:2] 2 15
# $ quadratic2:List of 6
#  ..$ : num [1:2] 6 11
#  ..$ : num [1:2] 7 8
#  ..$ : num [1:2] 9 10
#  ..$ : num [1:2] 9 13
#  ..$ : num [1:2] 12 15
#  ..$ : num 14

m <- minForest(dsMixed,homog=FALSE,stat="LR")
str(modelFormula(m))
#List of 4
# $ discrete  :List of 4
#  ..$ : int [1:2] 1 3
#  ..$ : int [1:2] 3 4
#  ..$ : int [1:2] 3 5
#  ..$ : int [1:2] 5 2
# $ linear    :List of 10
#  ..$ : num [1:2] 2 6
#  ..$ : num [1:2] 5 7
#  ..$ : num [1:2] 5 8
#  ..$ : num [1:2] 5 9
#  ..$ : num [1:2] 5 10
#  ..$ : num [1:2] 2 11
#  ..$ : num [1:2] 3 12
#  ..$ : num [1:2] 2 13
#  ..$ : num [1:2] 5 14
#  ..$ : num [1:2] 2 15
# $ quadratic : list()
# $ quadratic2: list()
}
\keyword{graphs}
