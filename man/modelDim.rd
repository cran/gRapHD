\name{modelDim}
\alias{modelDim}
\title{Model's dimension}
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
Gabriel Coelho Goncalves de Abreu (\email{abreu_ga@yahoo.com.br}) \cr
Rodrigo Labouriau (\email{Rodrigo.Labouriau@agrsci.dk}) 
}

\references{
Lauritzen, S.L. (1996) \emph{Graphical Models}, Oxford University Press.\cr
}


\examples{
data(dsCont)
m <- minForest(dsCont,stat="BIC")
modelDim(m)
# 102

m <- stepw(m,dsCont,stat="BIC")
modelDim(m)
# 149
}
\keyword{graphs}
