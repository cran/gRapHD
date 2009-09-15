\name{stepb}
\alias{stepb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stepwise backward selection
}
\description{
A function to perform stepwise backward selection to minimize AIC or BIC.
}
\usage{
stepb(G, dataset, fixed.edges = NULL, stat = "BIC")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{
A gRapHD object, the initial model.
}
  \item{dataset}{
A dataframe, whose variables correspond to the nodes of G.
}
  \item{fixed.edges}{
A boolean vector of length nrow(G@edges). Edges that are TRUE are not removed from the model.
}
  \item{stat}{
The measure to be minimized, either "BIC" (default) or "AIC".
}
}
\details{
Performs backward stepwise selection. The initial model must be decomposable. Only edges preserving decomposability are
eligible for removal. At each step, the edge resulting in the greatest reduction in BIC (or AIC) is removed. A pure graphical model
(i.e. either discrete or continuous) is decomposable iff its graph is triangulated. A mixed graphical models (i.e., with
both discrete and continuous variables) is decomposable iff its graph is triangulated and contains no forbidden paths.
}
\value{
A gRapHD object.
}
\author{
David Edwards (\email{David.Edwards@agrsci.dk})
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(attitude)
am <- NULL
for (i in 1:6) for (j in (i+1):7) am <- rbind(am, c(i,j))
satG <- new("gRapHD", edges=am, p=7, homog=TRUE, numCat=rep(0,7), vertNames=names(attitude))
G <- stepb(satG, attitude)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
