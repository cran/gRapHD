\name{CI.test}
\alias{CI.test}
\title{Test of conditional independence}
\description{
  Test of conditional independence.
}

\usage{
CI.test(x,y,S,dataset,homog=TRUE)
}

\arguments{
  \item{x}{one of the variables.}
  \item{y}{the other variable.}
  \item{S}{separator (possibly \code{NULL}).}
  \item{dataset}{matrix or data frame (\code{nrow(dataset)} observations and
                 \code{ncol(dataset)} variables).}
  \item{homog}{\code{TRUE} for homogeneous covariance structure, \code{FALSE}
               for heterogeneous. This is only meaningful with mixed models.
               Default is homogeneous (\code{TRUE}).}
}

\value{
  A list with the deviance (\code{deviance}) and the adjusted degrees of freedom
  (\code{numP}).
}

\details{
Performs a test of conditional independence of x and y given a set of variables S. The variables are specified as
column numbers of the dataset. Under the alternative the variables are assumed to follow an unrestricted
(mixed) graphical model. If x and y are discrete, S must also be discrete.  
Note that the model dimension returned by the \code{\link[gRapHD:fit]{fit}}
function assumes that all parameters are estimable, which may not be the case for
high-dimensional sparse data. However, here and in the search functions we use
the adjusted degrees of freedom, which need no such assumptions and are believed to be correct.
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk}) \cr
Rodrigo Labouriau (\email{Rodrigo.Labouriau@agrsci.dk}) \cr
David Edwards (\email{David.Edwards@agrsci.dk})
}

\references{
Lauritzen, S.L. (1996) \emph{Graphical Models}, Oxford University Press.\cr
Edwards, D. (2000) \emph{Introduction to Graphical Modelling}, Springer-Verlag
New York Inc.\cr
}

\examples{
data(dsCont)
m1 <- minForest(dsCont,homog=TRUE,forbEdges=NULL,stat="BIC")
CI.test(20,29,c(9,11),dsCont)
#$deviance
#[1] 0.7617515263220724
#
#$numP
#[1] 1
}
\keyword{graphs}
