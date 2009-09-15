\name{fit}
\alias{fit}
\title{Log-likelihood, AIC, BIC}
\description{
  Calculate -2*log-likelihood, AIC, and BIC for a triangulated graph (decomposable model).
}

\usage{
  fit(model=NULL, edges=NULL, dataset, homog=TRUE)
}

\arguments{
  \item{model}{\code{gRapHD} object.}
  \item{edges}{matrix with 2 columns, each row representing one edge, and each
               column one of the vertices in the edge.}
  \item{dataset}{matrix or data frame (\code{nrow(dataset)} observations and
                 \code{ncol(dataset)} variables).}
  \item{homog}{\code{TRUE} if the mixed model is homogeneous.}
}

\value{
  Vector with: model dimension (no of free parameters), -2*log-likelihood, AIC, and BIC. Note that
  all parameters are assumed to be estimable in the dimension calculation. 
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
  data(dsCont)
  m1 <- minForest(dsCont,homog=TRUE,forbEdges=NULL,stat="LR")
  fit(edges=m1$edges,dataset=dsCont)
}
\keyword{graphs}
