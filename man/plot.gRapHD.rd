\name{plot.gRapHD}
\alias{plot.gRapHD}
\title{Plot a gRapHD object}
\description{
  Plots a graph using the Fruchterman-Reingold algorithm. \code{\link{Methods}}
  for class \code{gRapHD}
}

\usage{
\method{plot}{gRapHD}(x,vert=NULL,numIter=50,main="",
     plotVert=TRUE,labelVert=TRUE,energy=FALSE,
     useWeights=FALSE,vert.hl=NULL,col.hl="red",
     vert.radii=0.01,coord=NULL,col.ed="darkgray",lty.ed=1,
     lwd.ed=1,lwd.vert=1,border=0,symbol.vert=1,
     cex.vert.label=.40,vert.labels=NULL,asp=NA,disp=TRUE,
     font=par("font"),...)
}

\arguments{
  \item{x}{a \code{gRapHD} object.}
  \item{vert}{list of vertices to be plotted. If \code{NULL}, all vertices will be
           plotted}
  \item{numIter}{number of iterations for the Fruchterman-Reingold algorithm.}
  \item{main}{main title.}
  \item{plotVert}{if \code{TRUE} the vertices are plotted.}
  \item{labelVert}{if \code{TRUE} the vertices are labelled.}
  \item{energy}{if \code{TRUE} use the minimum energy as initial values.}
  \item{useWeights}{if \code{TRUE} use the \code{model$statSeq} as edge length
                    weights).}
  \item{vert.hl}{list of vertices to be highlighted.}
  \item{col.hl}{colour to be used in the highlighted vertices.}
  \item{vert.radii}{radii of the edges (scalar or vector with length equal to the number
            of vertices). See 'Details'.}
  \item{coord}{initial coordinate values for the vertices.}
  \item{col.ed}{colour of the edges.}
  \item{lty.ed}{type of line for the edges.}
  \item{lwd.ed}{width of the line for the edges.}
  \item{lwd.vert}{width of the vertices' border.}
  \item{border}{colour of the vertices borders.}
  \item{symbol.vert}{symbol to be used in the vertices (length 1 or number of vertices):
            0 is a ellipse (\code{a=2*vs},\code{b=vs}), 1 a circle, 2 a square,
            3 or higher represents the number of sides.}
  \item{cex.vert.label}{numeric character expansion factor for the labels; multiplied by
              \code{\link[graphics:par]{par}} yields the final character size.
              \code{NULL} and \code{NA} are equivalent to \code{1.0}.}
  \item{vert.labels}{labels to be used in the vertices. If \code{NULL}, the vertices
               numbers are used.}
  \item{asp}{numeric, giving the aspect ratio y/x
             (see \code{\link[graphics:plot.window]{plot.window}} for more details).}
  \item{disp}{if \code{TRUE} (default), the graph is plotted.}
  \item{font}{an integer which specifies which font to use for the labels. If
              possible, device drivers let 1 correspond to plain
              text (the default), 2 to bold face, 3 to italic and 4 to bold
              italic.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
  If \code{ret} is \code{TRUE}, the coordinates for the vertices are returned.
}

\details{
  Plot a graph based on the list of edges.\cr
  Only one (\code{model} or \code{edges}) should be provided. If \code{model},
  the function uses also the information about the type of variables (discrete
  or continuous). If \code{edges}, then all variables are plotted as continuous
  (circles).\cr
  The plotting area is square, ranging from 0 to 1. The unit of parameter
  \code{vs} follow the axes. \cr
  The algorithm proposed by Fruchterman & Reingold (1991) is used to determine
  the position of each vertex. It is not initialised randomly, but in a regular
  grid.
}

\author{
Gabriel Coelho Goncalves de Abreu (\email{Gabriel.Abreu@agrsci.dk}) \cr
Rodrigo Labouriau (\email{Rodrigo.Labouriau@agrsci.dk}) \cr
David Edwards (\email{David.Edwards@agrsci.dk})
}

\references{
  Fruchterman, T.M.J. and Reingold, E.M. (1991) Graph Drawing by Force-directed
  Placement. \emph{SOFTWARE-PRACTICE AND EXPERIENCE}, VOL. 21(11), 1129-1164.\cr
  Csardi G, Nepusz T: The igraph software package for complex network
  research, InterJournal, Complex Systems 1695. 2006.
  http://igraph.sf.net
}

\examples{
  data(dsCont)
  m1 <- minForest(dsCont,homog=TRUE,forbEdges=NULL,stat="LR")
  plot(m1,numIter=1000)

  # or
  plot(m1,numIter=1000,plotVert=FALSE,labelVert=FALSE)
  
  #############
  r <- 3
  edges <- rep(1,r)
  x <- 2+r-1
  edges <- c(edges,sort(rep(2:x,r-1)))
  edges <- c(edges,sort(rep((x+1):(x+(x-1)*(r-1)),r-2)))
  edges <- c(edges,sort(rep((x+(x-1)*(r-1)+1):(x+(x-1)*(r-1)+(x-1)*(r-1)*(r-2)),
                             r-3)))
  edges <- cbind(edges,2:(length(edges)+1))
  a <- neighbourhood(edges=edges,orig=1,rad=r)
  vs <- a$v[,2]
  vs <- 1/vs
  vs[1] <- 2
  vs <- vs/30
  model <- as.gRapHD(edges)
  plot(model,numIter=200,vert.hl=a$v[,1],col.hl=colours()[386:383][a$v[,2]+1],
       vert.radii=vs,border="black",lwd.vert=2)
}
\keyword{graphs}
\keyword{dplot}
