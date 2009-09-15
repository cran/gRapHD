################################################################################
#   This file is part of gRapHD R package.
#
#   gRapHD R package
#   Copyright (C) 2009 Gabriel Coelho Goncalves de Abreu, Rodrigo Labouriau,
#   and David Edwards
#
#   gRapHD R package is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the Free
#   Software Foundation, either version 3 of the License, or any later version.
#
#   gRapHD R package program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

################################################################################
# Functions in this file:
#
#  SubGraph <- function(model=NULL,edges=NULL,v=NULL,p=0)
#  sortMat <- function(mat, cols)
#  MCS <- function(model=NULL,edges=NULL,v=0,p=NULL)
#  DFS <- function(model=NULL,edges=NULL, v, p=NULL)
#  minForest <- function(dataset,homog=TRUE,forbEdges=NULL,stat="BIC")
#  randTree <- function(p,seed=1)
#  chStat <- function(model,dataset,previous=NULL,forbEdges=NULL)
#  stepw <- function(model,dataset,stat="BIC",saveCH=NULL,forbEdges=NULL,
#                    exact=FALSE,initial=NULL,threshold=0,join=FALSE)
#  neighbours <- function(model=NULL,edges=NULL,v)
#  findEd <- function(edges,p,previous=NULL,numCat,from=0,exact=FALSE,join=F)
#  perfSets <- function(model=NULL,edges=NULL,p=NULL,varType=0,from=0)
#  rowProds <- function(x,na.rm=TRUE)
#  calcStat <- function(dataset,homog=TRUE,forbEdges=NULL,stat="LR")
#  plot.gRapHD <- function(x,vert=NULL,numIter=50,main="",
#                          plotVert=TRUE,labelVert=TRUE,energy=FALSE,
#                          useWeights=FALSE,vert.hl=NULL,col.hl="red",
#                          vert.radii=0.01,coord=NULL,col.ed="darkgray",lty.ed=1,
#                          lwd.ed=1,lwd.vert=1,border=0,symbol.vert=1,
#                          cex.vert.label=.40,vert.labels=NULL,asp=NA,disp=TRUE,
#                          font=par("font"),...)
#  neighbourhood <- function(model=NULL,edges=NULL,orig=NULL,rad=1)
#  adjMat <- function(model=NULL,edges=NULL,p=NULL)
#  Degree <- function(model=NULL,edges=NULL,v=NULL)
#  shortPath <- function(model=NULL, edges=NULL, v=NULL, p=NULL)
#  convData <- function(dataset)
#  fit <- function(model=NULL,edges=NULL,dataset,homog=TRUE)
#  summary.gRapHD <- function(object,...)
#  modelFormula <- function(model)
#  modelDim <- function(model)
#  as.gRapHD <- function(object,...)
#  print.gRapHD <- function(x,...)
#  is.gRapHD <- function(object)
#  jTree <- function(model)
#  ccoeff <- function(model=NULL,edges=NULL,p=NULL)
#  CI.test <- function(x,y,S,dataset,homog=TRUE)
#  seqLevels <- function(k,numCat)
#  normDens <- function (x, meanVec, covMat, logx = FALSE)
#  convertClass <- function(object)
################################################################################

################################################################################
# gRapHD object:
#   edges -  matrix with 2 columns, each row representing one edge, and each
#             column one of the vertices in the edge. Column 1 contains the
#             vertex with lower index.
#   p - number of variables (vertices) in the model.
#   stat.minForest - measure used (LR, AIC, or BIC).
#   stat.stepw - measure used (LR, AIC, or BIC).
#   statSeq - vector with value of stat.minForest for each edge.
#   numCat - vector with number of levels for each variable (0 if continuous).
#   homog - TRUE if the covariance is homogeneous.
#   numP - vector with number of estimated parameters for each edge.
#   minForest - first and last edges found with minForest.
#   stepw - first and last edges found with stepw.
#   userDef - first and last edges defined by the user.
################################################################################

################################################################################
# Description: Builds a sub-graph based on a list of vertices.
# In: model = gRapHD object
#     edges = matrix (n by 2, numeric).
#     v = vector k (vertices)
#     p = number of vertices
# Out: gRapHD object
################################################################################
SubGraph <- function(model=NULL,edges=NULL,v=NULL,p=0)
{
  if ((is.null(model) && is.null(edges)) || (!is.null(model) && !is.null(edges)))
    stop("Either model or edges must be provided.")
  if (is.null(edges))
  {
    if (class(model)!="gRapHD")
      stop("Model must be of class gRapHD.")
    result <- model
    edges <- model$edges
  }
  else
  {
    p <- ifelse(NROW(edges)==0,p,max(max(edges),p))
    result <- as.gRapHD(edges,p=p)
  }

  if (NROW(edges)==0)
    return(as.gRapHD(NULL,p=result$p,numCat=result$numCat,homog=result$homog))

  if (!is.null(v) & length(v)!=0)
  {
    v <- sort(unique(v))
    if (max(v)>result$p)
      stop(paste("v must be in [1,",result$p,"].",sep=""))
    ind <- apply(edges,1,function(x,y){length(intersect(x,y))==2},y=v)
    result$edges <- matrix(edges[ind,],ncol=2)
    result$statSeq <- result$statSeq[ind]
    result$numP <- result$numP[ind]
    result$vertNames <- result$vertNames[v]
    result$numCat <- result$numCat[v]
    aux <- rep(0,result$p)
    aux[v] <- 1:length(v)
    result$edges <- cbind(aux[result$edges[,1]],aux[result$edges[,2]])
    rm(aux)
    ind <- (1:nrow(edges))[ind]
    # if minForest was used, the resulting edges are the first in the list
    # and after come the edges from stepw
    if (!is.null(result$minForest))
    {
      if (length(ind)==0)
        result$minForest <- NULL
      else
        if (ind[1]>result$minForest[2])
          result$minForest <- NULL
        else
        {
          x <- which(ind<=result$minForest[2])
          x <- x[length(x)]
          result$minForest <- c(1,x)
        }
    }
    if (!is.null(result$stepw))
    {
      if (length(ind)==0)
        result$stepw <- NULL
      else
        if (ind[length(ind)]<result$stepw[1])
          result$stepw <- NULL
        else
        {
          x <- which(ind<=result$stepw[2])
          x <- x[length(x)]
          result$stepw <- c(x,length(ind))
        }
    }
    if (!is.null(result$userDef))
      if (length(ind)>0)
        result$userDef <- c(1,length(ind))
    result$p <- length(v)
  }
  return(result)
}
################################################################################

################################################################################
# Sorts the rows of a matrix, using the cols (in sequence).
# In: mat = matrix
#     cols = vector numeric (sequence of columns to be used in the sorting)
# Out: sorted matrix
################################################################################
sortMat <- function(mat, cols)
{
  m <- do.call("order", as.data.frame(mat[, cols]))
  mat[m, ]
}
################################################################################

################################################################################
# Maximum Cardinality Search (Tarjan & Yannakakis, in Jayson Rome, 2002).
# Searches for a perfect numbering.
# In: model = gRapHD object
#     edges = matrix (n by 2, numeric).
#     v     = starting vertex;
#     p     = number of vertices.
# Out: vector = (p by 1) where each position represents a
#               variable, containing the respective number in the perfect
#               numbering; os 0 is the graph is not triangulated.
################################################################################
MCS <- function(model=NULL,edges=NULL,v=0,p=NULL)
{
  if ((is.null(model) && is.null(edges)) || (!is.null(model) && !is.null(edges)))
    stop("Either model or edges must be provided.")
  if (is.null(edges))
  {
    if (class(model)!="gRapHD")
      stop("Model must be of class gRapHD.")
    edges <- model$edges
    p <- model$p
  }
  else
    v <- max(c(v,1))  

  if (NROW(edges)==0 && is.null(p))
    stop("p must not be NULL.")
  if (is.null(p))
    p <- max(edges)
  v <- ifelse(v==0,1,v)
  if (v > p)
    stop(paste("v must be in [1,",p,"].",sep=""))
  v1 <- edges[,1]
  v2 <- edges[,2]
  storage.mode(v1) <- storage.mode(v2) <- "integer"
  storage.mode(v) <- storage.mode(p) <- "integer"
  result <- .Call("mcs",v1,v2,v,p,PACKAGE="gRapHD")
  return(as.vector(result))
}
################################################################################

################################################################################
# Depth-first search, searches for all vertices reacheable from one specific
# vertex (considering that there is no cycle).
# In: model = gRapHD object
#     edges = matrix (n by 2, numeric).
#     v    = starting vertex;
#     p    = number of vertices.
# Out: K = vector containing all the vertices reacheable from v.
################################################################################
DFS <- function(model=NULL,edges=NULL, v, p=NULL)
{
  if ((is.null(model) && is.null(edges)) || (!is.null(model) && !is.null(edges)))
    stop("Either model or edges must be provided.")
  if (is.null(edges))
  {
    if (class(model)!="gRapHD")
      stop("Model must be of class gRapHD.")
    edges <- model$edges
    p <- model$p
  }
  else
    if (NROW(edges)==0 || NCOL(edges)!=2)
      stop("No edges.")
  if (is.null(p))
    p <- max(edges)
  if (v>p)
    stop(paste("v must be in [1,",p,"].",sep=""))

  v1 <- edges[,1]
  v2 <- edges[,2]
  storage.mode(v1) <- storage.mode(v2) <- "integer"
  storage.mode(v) <- storage.mode(p) <- "integer"
  result <- .Call("dfs",v1,v2,v,p,PACKAGE="gRapHD")
  rm(v1,v2)
  return(as.vector(result))
}
################################################################################

################################################################################
# Searches for the tree that minimizes the LR.
# In: dataset = matrix n by p.
#     homog = boolean - if the covariance is homogeneous
#     forbEdges = matrix k by 2, forbidden edges
#     stat = "LR" (default), "BIC", or "AIC"
# As defined in Lauritzen (1996, pg 168): saturated mixed model - lets X have
# an arbritary CG distribution. The homogeneous saturated mixed model restricts
# the distribution of X to a CG distribution which is homogeneous, i.e. the
# conditional covariance given the discrete variables does not depend on the
# conditional.
# Out: list = edges (matrix k by 2 with the edges, integer);
#             logL (numeric, log L of the selected model)
#             DF (integer, degrees of freedom) = p*(p-1)/2 - (p-1)
#             p (integer, number of variables)
#             stat.minForest(string, indicates which was used to build: LR,AIC,
#                              BIC)
#             end (integer, the last edge in the tree/forest)
#             statSeq (weight of each edges - LR,BIC,AIC)
#             numCat (vector, number of catagories in each discrete variable;
#                     if the variable is continuous numCat=0)
#             homog (boolean - if the covariance is homogeneous)
#             numP (vector, number of parameters in each edge)
################################################################################
minForest <- function(dataset,homog=TRUE,forbEdges=NULL,stat="BIC",cond=NULL,...)
{
  if (length(stat)==0)
    stop("No valid measure. Options: LR, BIC, AIC, or a user defined function.")
  if (mode(stat)=="function")
  {
    FUN <- match.fun(stat)
    stat <- 3
  }
  else
  {
    stat <- toupper(stat)
    if (is.element(stat,c("LR","AIC","BIC")))
      stat <- switch(stat,LR=0,BIC=1,AIC=2)
    else
      stop("No valid measure. Options: LR, BIC, AIC, or a user defined function.")
  }

  gClique <- function(v)
  {
    v <- sort(unique(v))
    x <- function(v,y)
    {
      if (v<y)
        return(cbind(v,(v+1):y))
      else
        return(NULL)
    }
    n <- length(v)
    aux <- sapply(1:n,x,y=n)
    edges <- NULL
    for (i in 1:length(aux))
      edges <- rbind(edges,matrix(v[aux[[i]]],,2))
    return(edges)
  }

  p <- ncol(dataset)
  n <- nrow(dataset)
  statS <- stat
  aux <- convData(dataset)
  dataset <- aux$ds
  numCat <- aux$numCat
  result <- aux$vertNames
  rm(aux)

  if (sum(is.na(dataset))>0)
    stop("Missing values not allowed.")

  if (is.null(forbEdges))
    forbEdges <- 0

###################################
  condEdges <- NULL
  compType <- varType <- as.integer(numCat!=0)
  if (is.null(cond))
    comp <- 1:p
  else
  {
    k <- length(cond) #number of conditional sets
    ed <- matrix(integer(0),,2)
    if (k > 1)
    {
      inter <- matrix(0,k,k) #matrix giving which sets are not disjoint
      for (i in 1:(k-1))
        for (j in (i+1):k)
          inter[i,j] <- length(intersect(cond[[i]],cond[[j]]))>0
      ed <- which(inter==1,arr.ind=T) #edges of a graph in which the nodes are de sets
    }
    if (nrow(ed) == 0)
      condJ <- cond
    else
    {
      compJ <- rep(0,k) #indicates which sets are in the same connected component in ed
      i <- 1
      while (min(compJ)==0)
      {
        aux <- DFS(edges=ed,v=i,p=k) #all nodes reachable from i
        if (compJ[i] == 0)
          compJ[c(i,aux)] <- max(compJ)+1
        i <- i + 1
      }
      condJ <- list() #the list of nodes in the original graph that are in the same connected component
      for (i in 1:max(compJ))
        condJ[[i]] <- sort(unique(unlist(cond[which(compJ==i)])))
      #  forbEdges <- rbind(forbEdges,gClique(unlist(cond[which(compJ==i)])))
    }

    comp <- rep(0,p)
    for (i in 1:length(condJ))
    {
      comp[condJ[[i]]] <- max(comp)+1
      x <- unique(varType[condJ[[i]]])
      compType[condJ[[i]]] <- ifelse(length(x)==1,x,2)
      forbEdges <- rbind(forbEdges,gClique(condJ[[i]]))
    }
    for (i in 1:length(cond))
      condEdges <- rbind(condEdges,gClique(cond[[i]]))
    condEdges <- unique(condEdges,MARGIN=1)
    ind <- which(comp==0)
    if (length(ind)>0)
    {
      comp[ind] <- max(comp) + (1:length(ind))
    }
  }
###################################
  if (forbEdges[1] != 0)
  {
    forbEdges <- matrix(forbEdges,,2)
    X <- function(x){sort(x)}
    forbEdges <- t(apply(forbEdges,1,X))
    x <- forbEdges[,1]
    y <- forbEdges[,2]
    forbEdges <- (x-1)*p-(x-1)*x/2 + y-x #simpler representation of the edges
  }
  values <- NULL
  if (stat == 3)
  {
    for (i in 1:(p-1))
      for (j in (i+1):p)
      {
        x <- (i-1)*p-(i-1)*i/2 + j-i
        if (!is.element(x,forbEdges))
          values <- rbind(values,FUN(c(i,j),numCat,dataset,...))
        else
          values <- rbind(values,NA)
      }
  }

  storage.mode(values) <- storage.mode(dataset) <- "double"
  storage.mode(numCat) <- storage.mode(homog) <- "integer"
  storage.mode(forbEdges) <- storage.mode(stat) <- "integer"
  storage.mode(comp) <- storage.mode(compType) <- "integer"
  aux <- .Call("minForest",dataset,numCat,homog,
               forbEdges,stat,values,comp,compType,PACKAGE="gRapHD")
  tree <- aux$tree
  storage.mode(homog) <- "logical"
  n <- NROW(tree)

  if (is.null(condEdges))
  {
    edges <- matrix(tree[,1:2],n,2)
    statSeq <- tree[,3]
    numP <- tree[,4]
  }
  else
  {
    edges <- rbind(matrix(tree[,1:2],n,2),condEdges)
    statSeq <- c(tree[,3],rep(NA,nrow(condEdges)))
    numP <- c(tree[,4],rep(NA,nrow(condEdges)))
  }
  result <- list(edges = edges,
                p = p,
                stat.minForest = switch(stat+1,"LR","BIC","AIC","User's function"),
                statSeq = statSeq,
                numCat = numCat,
                homog = homog,
                numP = numP,
                vertNames = result,
                minForest = c(1,nrow(edges)))
  class(result)<-"gRapHD"
  if (dim(aux$errors)[1]!=0)
  {
    result$error <- aux$errors
    warning("Check model$errors for edges with problems.", call. = FALSE)
  }
  return(result)
}
################################################################################

################################################################################
# Generates a random tree (no cycles) with "p" vertices.
# In: p = integer (number of vertices)
#     seed = for the random generator.
# Out: list = edges (matrix p-1 by 2, with edges)
#             seed
#             p
################################################################################
randTree <- function(p,seed=1)
{
  set.seed(seed,kind="Mersenne-Twister")
  V <- 1:p
  r <- sample(1:p,1)
  V[c(r,p)] <- V[c(p,r)]
  tree <- matrix(0,nrow=p-1,ncol=2)
  for (i in 1:(p-1))
  {
    x <- ifelse(i==p-1,1,sample(1:(p-i),1))
    y <- ifelse(i==1,p,sample((p-i+1):p,1))
    tree[i,] <- sort(c(V[x],V[y]))
    V[c(x,p-i)] <- V[c(p-i,x)]
  }
  result <- list(edges=tree,seed=seed,p=p)
  return(result)
}
################################################################################

################################################################################
# Calculates the change in the model's BIC resulting by the addition of all
# possible edges (edges that preserve the graph's triangulation).
# In: model = gRapHD object
#     dataset = matrix n by p.
#     previous = same output of chStat()
#     forbEdges = matrix k by 2, forbidden edges
# As defined in Lauritzen (1996, pg 168): the saturated mixed model lets X have
# an arbritary CG distribution. The homogeneous saturated mixed model restricts
# the distribution of X to a CG distribution which is homogeneous, i.e. the
# conditional covariance given the discrete variables does not depend on the
# conditional.
# Out: list:
#      edges.to.test = matrix (k by 4), column:
#                      1 - first vertex of the tested edge
#                      2 - second vertex (the two values are ordered)
#                      3 - index for the separator in S
#                      4 - change in the LR by adding this edge
#                      5 - number of parameters for that edge
#      S = list with the separators
################################################################################
chStat <- function(model,dataset,previous=NULL,forbEdges=NULL)
{
  if (sum(is.na(dataset))>0)
    stop("Missing values not allowed.")

  p <- ncol(dataset) # number of variables (vertices)
  n <- nrow(dataset) # number of observations
  edges.to.test <- previous$edges.to.test # (v1; v2; index of the minimal separator;
                             # a previous value for the LR change, or 0
                             # indicating that it has to be calculated;
                             # the number of parameters for that edge)
  SS <- previous$S # minimal separators
  rm(previous)

  num <- nrow(edges.to.test)

  numbPar <- function(y)
  {
    a <- sum(rowSums(y)!=0)-1
    b <- sum(colSums(y)!=0)-1
    return(a*b*(a>0))
  }
  whichSeq <- function(x,y){return(sum(abs(x-y))==0)}
  Gcliq <- function(v)
  {
    v <- sort(unique(v))
    i <- 1
    n <- length(v)
    ed <- sort(rep(i:(i+n-1),n))
    ed <- cbind(ed,rep(i:(i+n-1),n))
    ed <- matrix(ed[-unlist(lapply(1:n,function(x,n){1:x+(x-1)*n},n=n)),],,2)
    dimnames(ed) <- NULL
    ed[,1] <- v[ed[,1]]
    ed[,2] <- v[ed[,2]]
    return(ed)
  }

  if (num > 0)
    for (i in 1:num)
    {
      x <- edges.to.test[i,1]
      y <- edges.to.test[i,2]
      if (!is.na(edges.to.test[i,4]) && is.finite(edges.to.test[i,4]))
        if ((edges.to.test[i,4]==0) &
            (!is.element((x-1)*p-(x-1)*x/2+y-x,forbEdges)))
        { # the change has to be calculated
          # this is the clique that will result from the inclusion of this edge
          S <- SS[[edges.to.test[i,3]]]
          clique <- c(edges.to.test[i,1:2],S)

          # the 2 cliques thal will be merged
          clique1 <- c(edges.to.test[i,1],S)
          clique2 <- c(edges.to.test[i,2],S)

          if (sum(model$numCat[clique]!=0) == 0) # continuous
          {
            CM <- cov(dataset[,clique],use="pairwise.complete.obs")*(nrow(dataset)-1)
            if (sum(is.na(CM))) # meaning that the new edge does not have more
                                # than 1 observation (not NA)
            {
              edges.to.test[i,4] <- 0
              edges.to.test[i,5] <- 0
            }
            else
            {
              a <- match(edges.to.test[i,1:2],clique)
              # because of join in findEd
              d <- log(det(CM)) + ifelse(length(S)==0,0,log(det(matrix(CM[-a,-a],length(clique)-2))))
              ##
              d <- d - (log(det(matrix(CM[-a[1],-a[1]],length(clique)-1))) +
                        log(det(matrix(CM[-a[2],-a[2]],length(clique)-1))))
              edges.to.test[i,4] <- nrow(dataset)*d
              edges.to.test[i,5] <- 1
            }
          }
          else
            if (sum(model$numCat[clique]==0)==0) # discrete
            {
              t12 <- table(as.data.frame(dataset[,clique]))
              t1 <- margin.table(t12,match(clique1,clique))
              t2 <- margin.table(t12,match(clique2,clique))
              # because of join in findEd
              tS <- if(length(S)==0) 0 else margin.table(t12,match(S,clique))
              # to calculate the df: Whittaker, pg 223, Proposition 7.6.2 and
              # pg 224, Corolary 7.6.3
              # using MARGIN=(1:length(clique))[-(1:2)] uses the other margins but
              # not the first 2 (that are the 2 vertices being tested)

              # because of join in findEd
              if (length(S)==0)
                numP <- prod(model$numCat[clique]-1)#length(t12)-1
              else
                numP <- sum(apply(t12,MARGIN=(1:length(clique))[-(1:2)],numbPar))
              ##
              numObsUsed <- sum(t12)
              t12 <- t12[t12!=0]
              t1 <- t1[t1!=0]
              t2 <- t2[t2!=0]
              tS <- tS[tS!=0]

              if (numObsUsed>0)
              {
                # BIC(n+1) - BIC(n) = log(f(C)/(f(C1)*f(C2)/f(S)))
                v12 <- sum(t12*log(t12/numObsUsed))
                v1  <- sum(t1 *log(t1 /numObsUsed))
                v2  <- sum(t2 *log(t2 /numObsUsed))
                # because of join in findEd
                vS  <- ifelse(length(S)==0,0,sum(tS *log(tS /numObsUsed)))
                ##
                edges.to.test[i,4] <- -2*(v12 + vS - v1 - v2)
                edges.to.test[i,5] <- numP
              }
              else
              {
                edges.to.test[i,4] <- 0
                edges.to.test[i,5] <- 0
              }
            }
            else # mixed
            {
              discrete   <- intersect(clique,which(model$numCat!=0))
              continuous <- setdiff(clique,discrete)
              tab <- table(as.data.frame(dataset[,discrete]))
              ssd <- diag(0,length(continuous))
              if (sum(model$numCat[edges.to.test[i,1:2]])==0) # both are continuous
                if (model$homog)
                {
                  for (j in 1:length(tab))
                    if (tab[j]>0)
                    {
                      ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,model$numCat[discrete]))
                      if (sum(ind)>1) #otherwise there is only one observation, and the contribution to the variance is zero
                        ssd <- ssd + cov(dataset[ind,continuous],use="pairwise.complete.obs")*(tab[j]-1)
                    }
                  ind <- match(edges.to.test[i,1:2],continuous)
                  edges.to.test[i,4] <- n*(log(det(ssd)) +
                                            ifelse(length(continuous)>2,log(det(matrix(ssd[-ind,-ind],length(continuous)-2))),0) -
                                            log(det(matrix(ssd[-ind[1],-ind[1]],length(continuous)-1))) -
                                            log(det(matrix(ssd[-ind[2],-ind[2]],length(continuous)-1))))
                  edges.to.test[i,5] <- 1
                }
                else
                {
                  edges.to.test[i,4] <- 0
                  ind1 <- match(edges.to.test[i,1:2],continuous)
                  for (j in 1:length(tab))
                    if (tab[j]>0)
                    {
                      ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,model$numCat[discrete]))
                      if (sum(ind)>1) #otherwise there is only one observation, and the contribution to the variance is zero
                        ssd <- cov(dataset[ind,continuous],use="pairwise.complete.obs")*(tab[j]-1)
                      edges.to.test[i,4] <- edges.to.test[i,4] +
                                            +tab[j]*(log(det(ssd)) +
                                                    ifelse(length(continuous)>2,log(det(matrix(ssd[-ind1,-ind1],length(continuous)-2))),0) -
                                                    log(det(matrix(ssd[-ind1[1],-ind1[1]],length(continuous)-1))) -
                                                    log(det(matrix(ssd[-ind1[2],-ind1[2]],length(continuous)-1))))
                    }
                  edges.to.test[i,5] <- sum(tab>(length(continuous)+1))
                }
              else
              {
                vDiscr <- (edges.to.test[i,1:2])[model$numCat[edges.to.test[i,1:2]]!=0]
                vCont  <- setdiff(edges.to.test[i,1:2],vDiscr)
                discrMarg <- setdiff(discrete,vDiscr)
                tabMarg <- margin.table(tab,match(discrMarg,discrete))
                ssdMarg <- diag(0,length(continuous))
                if (model$homog)
                {
                  for (j in 1:length(tab))
                  {
                    if (tab[j]>0)
                    {
                      ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,model$numCat[discrete]))
                      if (sum(ind)>1) #otherwise there is only one observation, and the contribution to the variance is zero
                        ssd <- ssd + var(dataset[ind,continuous],use="pairwise.complete.obs")*(tab[j]-1)
                    }
                    if (j<=length(tabMarg))
                    {
                      if (tabMarg[j]>0)
                      {
                        if (length(discrMarg)>0)
                          ind <- apply(matrix(dataset[,discrMarg],,length(discrMarg)),1,whichSeq,y=seqLevels(j,model$numCat[discrMarg]))
                        else
                          ind <- rep(TRUE,tabMarg[j])
                        if (sum(ind)>1) #otherwise there is only one observation, and the contribution to the variance is zero
                          ssdMarg <- ssdMarg + var(dataset[ind,continuous],use="pairwise.complete.obs")*(tabMarg[j]-1)
                      }
                    }
                  }
                  ind <- match(vCont,continuous)
                  edges.to.test[i,4] <- n*(log(det(ssd)) +
                                            ifelse(length(continuous)>1,log(det(matrix(ssdMarg[-ind,-ind],length(continuous)-1))),0) -
                                            ifelse(length(continuous)>1,log(det(matrix(ssd[-ind,-ind],length(continuous)-1))),0) -
                                            log(det(ssdMarg)))
                  if (length(discrete)==1) #means that the separator is continuous
                    edges.to.test[i,5] <- model$numCat[vDiscr]-1
                  else
                  {
                    aux <- match(vDiscr,discrete) #position of the discrete variable in the new edge (in the list of discrete variables in the edge/separator)
                    # dimension 1 in aux1 is the discrete variable in the tested edge
                    aux1 <- table(as.data.frame(dataset[,c(vDiscr,discrete[-aux])])) #it's only the table with the discrete in the edge as the 1st dimension
                    aux.tab <- margin.table(aux1,2:length(discrete)) #the marginal table without the edge's discrete
                    aux.ind <- which(aux.tab>0,arr.ind=TRUE) #cells with at least one observation
                    numP <- 0
for (iii in 1:nrow(aux.ind))
  numP <- numP + (sum(aux1[cbind(1:model$numCat[vDiscr],aux.ind[iii,])]>0) - 1)
                    rm(aux,aux1,aux.tab,aux.ind)
                    edges.to.test[i,5] <- numP
                  }
                }
                else
                {
                  edges.to.test[i,4] <- 0
                  ind1 <- match(vCont,continuous)
                  for (j in 1:length(tab))
                  {
                    if (tab[j]>0)
                    {
                      ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,model$numCat[discrete]))
                      ssd <- var(dataset[ind,continuous],use="pairwise.complete.obs")*(tab[j]-1)
                      dim(ssd) <- rep(length(continuous),2)
                    }
                    if (j<=length(tabMarg))
                    {
                      if (tabMarg[j]>0)
                      {
                        if (length(discrMarg)>0)
                          ind <- apply(matrix(dataset[,discrMarg],,length(discrMarg)),1,whichSeq,y=seqLevels(j,model$numCat[discrMarg]))
                        else
                          ind <- rep(TRUE,tabMarg[j])
                        ssdMarg <- var(dataset[ind,continuous],use="pairwise.complete.obs")*(tabMarg[j]-1)
                        dim(ssdMarg) <- rep(length(continuous),2)
                      }
                    }
                    ind <- match(vCont,continuous)
                    edges.to.test[i,4] <- edges.to.test[i,4] +
                                          -tab[j]*(log(det(ssd))-ifelse(length(continuous)>1,log(det(matrix(ssd[-ind,-ind],length(continuous)-1))),0)) +
                                          +ifelse(j<=length(tabMarg),tabMarg[j],0)*(log(det(ssdMarg))-ifelse(length(continuous)>1,log(det(matrix(ssdMarg[-ind,-ind],length(continuous)-1))),0))
                  }
                  edges.to.test[i,4] <- edges.to.test[i,4] +
                                        sum(tab*log(tab),na.rm=TRUE) -
                                        sum(tabMarg*log(tabMarg),na.rm=TRUE)
                  edges.to.test[i,4] <- -edges.to.test[i,4]
                  if (length(discrete)==1) #means that the separator is continuous
                    edges.to.test[i,5] <- sum(tab>(length(continuous)+1))-1
                  else
                  {
                    aux <- match(vDiscr,discrete)
                    # dimension 1 in aux1 is the discrete variable in the tested edge
                    aux1 <- table(as.data.frame(dataset[,c(vDiscr,discrete[-aux])]))
                    aux.tab <- margin.table(aux1,2:length(discrete))
                    aux.ind <- which(aux.tab>(length(continuous)+1),arr.ind=TRUE)
                    numP <- 0
for (iii in 1:nrow(aux.ind))
  numP <- numP + (sum(aux1[cbind(1:model$numCat[vDiscr],aux.ind[iii,])]>0) - 1)
                    rm(aux,aux1,aux.tab,aux.ind)
                    edges.to.test[i,5] <- numP
                  }
                }
              }
            }
        }
    }
  return(list(edges.to.test=edges.to.test,S=SS))
}
################################################################################

################################################################################
# Does a stepwise forward, starting from model and using the change in BIC to
# update the model. Tests every possible edge and selects the one that reduces
# the BIC the most. Stops when no more redution is possible, or there is no
# more possible edge to test.
# In: model = gRapHD object
#     dataset = matrix n by p
#     stat = "BIC","AIC","LR"
#     saveCH = if not NULL save each iteration's "chStat", with saveCH as
#              part of the file name
#     forbEdges = matrix k by 2, forbidden edges
#     exact = TRUE or FALSE
#     initial = continue the algorithm from a previous point. The paramenter
#                must have the same structure as the result of chStat.
#     threshold = values greater than the threshold are not included in the
#                 model
# Out: model (the same updated object)
################################################################################
stepw <- function(model,dataset,stat="BIC",saveCH=NULL,forbEdges=NULL,
                  exact=FALSE,initial=NULL,threshold=0,join=FALSE)
{
  if (class(model)!="gRapHD")
    stop("Model is not of class gRapHD.")
  else
  {
    star <- NULL
    aux <- range(model$numCat)
    if ((aux[1]==0) & (aux[2]!=0)) #mixed
    {
      # there are at least 2 discrete variables, which means that forbidden paths
      # are possible to happen
      # star graph
      star <- cbind(which(model$numCat!=0),model$p+1)
      perfect <- perfSets(edges=rbind(model$edges,star),
                          p=model$p+1,varType=c(model$numCat,1),from=0)
      if (!is.list(perfect))
        stop("The model is not strongly decomposable.")
      rm(perfect)
    }
    else
      if (MCS(edges=model$edges,,p=model$p)[1]==0)
        stop("The model is not triangulated.")
    rm(aux)
  }
  if (length(stat)==0)
    stop("No valid measure. Options: LR, BIC, AIC, or a user defined function.")
  if (mode(stat)=="function")
  {
    FUN <- match.fun(stat)
    model$stat.stepw <- "User's function"
    stat <- "USER"
  }
  else
  {
    stat <- toupper(stat)
    if (!is.element(stat,c("LR","AIC","BIC")))
      stop("No valid measure. Options: LR, BIC, AIC, or a user defined function.")
    model$stat.stepw <- stat
  }

  edges.stepw <- nrow(model$edges)+1

  if (model$p != NCOL(dataset))
    stop("model$p doesn't agree with NCOL(dataset).")

  aux <- convData(dataset)
  dataset <- aux$ds
  numCat <- aux$numCat
  if (!all.equal(numCat,model$numCat))
    stop("model$numCat doesn't agree with data.")
  if (!identical(sort(aux$vertNames),sort(model$vertNames)))
    stop("model$vertNames doesn't agree with data.")
  dataset <- dataset[,match(model$vertNames,aux$vertNames)]
  rm(aux)
  if (sum(is.na(dataset))>0)
    stop("Missing values not allowed.")

  p <- ncol(dataset) # number of variables (vertices)
  n <- nrow(dataset) # number of observations
  stat <- toupper(stat)

  if (is.null(forbEdges))
    forbEdges <- 0
  else
  {
    forbEdges <- matrix(forbEdges,,2)
    X <- function(x){sort(x)}
    forbEdges <- t(apply(forbEdges,1,X))
    x <- forbEdges[,1]
    y <- forbEdges[,2]
    forbEdges <- (x-1)*p-(x-1)*x/2 + y-x; #simpler representation of the edges
    rm(X,x,y)
  }

  STOP <- FALSE # flag indicating when to stop
  iteration <- 0
  if (!is.null(model$minForest))
    iteration <- nrow(model$edges) - model$minForest[2]
  ch <- initial
  SS <- NULL

  while ((!STOP))
  {
    iteration <- iteration + 1
    ch <- findEd(model$edges,p,ch,model$numCat,0,exact,join)
    if (stat!="USER")
      ch <- chStat(model,dataset,ch,forbEdges)
    else
      ch <- FUN(model,dataset,ch,forbEdges)
    change <- ch$edges.to.test
    if (!is.null(saveCH))
    {
      fileN <- paste(saveCH,"_",
                     paste(rep("0",6-(as.integer(log10(iteration))+1)),sep="",
                           collapse=""),iteration,".RData",sep="")
      save(ch,model,file=fileN)
    }
    if (nrow(change) == 0) # there is no more edges possible to test
      STOP <- TRUE
    else
    {
      continue <- TRUE
      add <- FALSE
      # before adding one the edge, it has to test if it generates a cycle
      df <- change[,5]
      df[df<=0] <- NA
      statValues <- change[,4] + ((stat=="BIC")*df+(stat=="AIC")*2)*log(n)
      rm(df)
      ind <- order(statValues)
      i <- 0
      while (continue)
      {
        # the edge which reduces the most (the stat)
        i <- i+1
        aux <- ind[i]
        value <- statValues[aux]
        if (is.na(value))
          continue <- FALSE
        else
          if (is.finite(value))
          {
            if ((value < threshold) & (exact))
            {
              continue <- FALSE
              add <- TRUE
            }
            else
              if (value < threshold)
              {
                # Lauritzen (1996, pg 11): Proposition 2.6 (Leimer) - An undirected,
                # marked graph G is decomposable if and only if G* is triangulated.
                continue <- (MCS(edges=rbind(model$edges,change[aux,1:2],star),v=1,p=p+!is.null(star)))[1] == 0
                add <- !continue #statValues[aux] <- !continue*statValues[aux] + continue*threshold
              }
              else
                continue <- FALSE
          }
      }
      if (add) # it improves the model
      {
        model$edges <- rbind(model$edges,change[aux,1:2])
        model$numP <- c(model$numP,change[aux,5])
        model$statSeq <- c(model$statSeq,change[aux,4])
      }
      else
        STOP <- TRUE
    }
  }
  if (edges.stepw<=nrow(model$edges))
    model$stepw <- c(edges.stepw,nrow(model$edges))
  else
    model$stepw <- integer(0)
  return(model)
}
################################################################################

################################################################################
# Finds the first neighbours of a vertex.
# In: model = gRapHD object
#     edges = matrix (n by 2, numeric).
#     v
# Out: fn = vector (first neighbours numbers, sortered)
################################################################################
neighbours <- function(model=NULL,edges=NULL,v)
{
  if ((is.null(model) && is.null(edges)) || (!is.null(model) && !is.null(edges)))
    stop("Either model or edges must be provided.")
  if (is.null(edges))
  {
    if (class(model)!="gRapHD")
      stop("Model must be of class gRapHD.")
    edges <- model$edges
    p <- model$p    
  }
  else
    if (NROW(edges)==0 || NCOL(edges)!=2)
      stop("No edges.")
    else
      p <- max(edges)
  if (v > p)
    stop(paste("v must be in [1,",p,"].",sep=""))  
  fn <- integer(0)
  edges <- matrix(edges,ncol=2)
  ind <- as.matrix(which(rbind(edges,c(0,0))==v,arr.ind=TRUE)) # edges in "v"
  n <- nrow(ind) # number of edges found
  if (n > 0)
    for (i in 1:n)
      fn <- c(fn,edges[ind[i,"row"],(1:2)[-ind[i,"col"]]])
  rm(ind,n)
  return(sort(fn))
}
################################################################################

################################################################################
# Given a triangulated graph (edges), get the sets:
#               C1,...,Ck: clique ordering (from a perfect numbering)
#               H1,...,Hk: histories
#               R1,...,Rk: residues
#               S1,...,Sk: separators
# From these sets, finds the possible edges (but here, it finds some edges that
# don't preserve the triangulation, so, before including a edge, it has to be
# tested). It also gives the minimal separator, what makes easier (and much
# faster) the calculations.
# If some previous data is provided, the last column in the resulting matrix
# indicates if the value should be recalculated or not.
# In: edges = matrix k by 2
#     p = integer (number of variables in the dataset)
#     previous = list: edges.to.test = matrix w by 4 (vertex 1; vertex 2; index
#                        for the separator in S; value of the change in the LR,
#                        number of parameters for that edge)
#                S: list with the minimal separators
#     varType = vertices type
#     from = first vertex in the perfect sequence
#     exact = TRUE or FALSE (if it's 1 (or TRUE) the QR decomp is used, if 2,
#             the gauss elimination)
#     join = T or F (is disjoind components can be joined)
# Out: list: edges = matrix w1 by 4 (same structure as edges.to.test)
#            S: list with the minimal separators
################################################################################
findEd <- function(edges,p,previous=NULL,numCat,from=0,exact=FALSE,join=FALSE)
{
  if (from > p)
    stop("from > p!!")
  if (is.null(edges))
    edges <- matrix(integer(0),1,2)
  if (is.null(previous))
    prev <- 0
  else
  {
    x <- previous$edges.to.test[,1]
    y <- previous$edges.to.test[,2]
    prev <- (x-1)*p-(x-1)*x/2 + y-x #simpler representation of the edges
  }

  v1 <- edges[,1]
  v2 <- edges[,2]
  storage.mode(v1) <- storage.mode(v2) <- "integer"
  storage.mode(p) <- storage.mode(prev) <- "integer"
  storage.mode(numCat) <- "integer"
  storage.mode(exact) <- storage.mode(from) <- "integer"
  storage.mode(join) <- "integer"

  result <- .Call("findEd", v1, v2, p, prev, numCat,
                  from, exact, join,PACKAGE="gRapHD")

  i <- result$total
  if (i > 0)
  {
    result <- list(edges.to.test=result$edges[1:result$total,1:5],S=result$S)
    dim(result$edges.to.test) <- c(i,5)
    result$edges.to.test[,5] <- 0
    if (!is.null(previous))
    {
      result$edges.to.test[result$edges.to.test[,4]>0,5] <-
            previous$edges.to.test[result$edges.to.test[,4],5]
      result$edges.to.test[result$edges.to.test[,4]>0,4] <-
            previous$edges.to.test[result$edges.to.test[,4],4]
    }
  }
  else
  {
    result$edges.to.test <- integer(0)
    dim(result$edges.to.test) <- c(0,5)
  }

  return(result)
}
################################################################################

################################################################################
# Given a triangulated graph (edges), get the sets:
#               C1,...,Ck: clique ordering (from a perfect numbering)
#               H1,...,Hk: histories
#               R1,...,Rk: residues
#               S1,...,Sk: separators
# In: model = gRapHD object
#     edges = matrix (n by 2, numeric).
#     p = integer (number of variables in the dataset)
#     varType = vertices type
#     from = first vertex in the perfect sequence
# Out: list - C,S,H,R
################################################################################
perfSets <- function(model=NULL,edges=NULL,p=NULL,varType=0,from=0)
{
  if ((is.null(model) && is.null(edges)) || (!is.null(model) && !is.null(edges)))
    stop("Either model or edges must be provided.")
  if (is.null(edges))
  {
    if (class(model)!="gRapHD")
      stop("Model must be of class gRapHD.")
    edges <- model$edges
    p <- model$p
    varType <- model$numCat
    varType[varType!=0] <- 1
  }
  if ((NROW(edges)==0) && is.null(p))
    stop("p must not be NULL.")
  else
    if (is.null(p))
      p <- max(edges)

  from <- max(c(0,from))
  if (is.null(edges))
    edges <- matrix(integer(0),,2)

  if (!is.element(length(varType),c(1,p)))
    stop("varType must be 0, 1, or lenght(varType)==p.")

  if (length(varType)==1)
    varType <- rep(as.integer(varType!=0),p)
  else
    varType[varType!=0] <- 1

  if (nrow(edges) > 0)
    if (max(edges)>p)
      stop("There are more vertices in edges than p.")
  v1 <- edges[,1]
  v2 <- edges[,2]
  storage.mode(v1) <- storage.mode(v2) <- "integer"
  storage.mode(p) <- storage.mode(varType) <- "integer"
  storage.mode(from) <- "integer"
  result <- .Call("perfSets",v1,v2,p,varType,from,PACKAGE="gRapHD")
  
  if (storage.mode(result)=="list")
    names(result$cliques) <- names(result$histories) <-
    names(result$separators) <- names(result$residuals) <- 1:length(result$cliques)

  return(result)
}
################################################################################

################################################################################
# Like "rowSums", evaluate the products of each row.
# In: x = matrix k by p
#     na.rm = TRUE/FALSE
# Out: vector k
################################################################################
rowProds <- function(x,na.rm=TRUE)
{
  if (is.data.frame(x))
      x <- as.matrix(x)
  if (!is.array(x) || length(dn <- dim(x)) < 2)
      stop("'x' must be an array of at least two dimensions")
      
  k <- nrow(x)
  p <- ncol(x)
  storage.mode(x) <- "double"
  storage.mode(k) <- "integer"
  storage.mode(na.rm) <- storage.mode(p) <- "integer"
  result <- .Call("rowProds",x,k,p,na.rm,PACKAGE="gRapHD")
  return(result)
}
################################################################################

################################################################################
# Calculates statistics (LR for independence, AIC, BIC) for each pair of
# variables in the dataset.
# In: dataset (n by p, numeric).
#     homog = T or F
#     forbEdges = list of forbidden edges
#     stat = LR, BIC, AIC
# Out: matrix = p*(p-1)/2 by 3
#               edge - columns 1 and 2
#               stat - column (2*log(L))
################################################################################
calcStat <- function(dataset,homog=TRUE,forbEdges=NULL,stat="LR")
{
  if (length(stat)==0)
    stop("No valid measure. Options: LR, BIC, AIC, or a user defined function.")
  if (mode(stat)=="function")
  {
    FUN <- match.fun(stat)
    stat <- 3
  }
  else
  {
    stat <- toupper(stat)
    if (is.element(stat,c("LR","AIC","BIC")))
      stat <- switch(stat,LR=0,BIC=1,AIC=2)
    else
      stop("No valid measure. Options: LR, BIC, AIC, or a user defined function.")
  }

  aux <- convData(dataset)
  dataset <- aux$ds
  numCat <- aux$numCat
  rm(aux)
  if (sum(is.na(dataset))>0)
    stop("Missing values not allowed.")

  p <- ncol(dataset)
  n <- nrow(dataset)

  if (is.null(forbEdges))
    forbEdges <- 0
  else
  {
    forbEdges <- matrix(forbEdges,,2)
    X <- function(x){sort(x)}
    forbEdges <- t(apply(forbEdges,1,X))
    x <- forbEdges[,1]
    y <- forbEdges[,2]
    forbEdges <- (x-1)*p-(x-1)*x/2 + y-x; #simpler representation of the edges
  }

  values <- NULL
  if (stat == 3)
  {
    for (i in 1:(p-1))
      for (j in (i+1):p)
      {
        x <- (i-1)*p-(i-1)*i/2 + j-i
        if (!is.element(x,forbEdges))
          values <- rbind(values,FUN(c(i,j),numCat,dataset))
        else
          values <- rbind(values,NA)
      }
  }

  storage.mode(values) <- storage.mode(dataset) <- "double"
  storage.mode(numCat) <- storage.mode(homog) <- "integer"
  storage.mode(forbEdges) <- storage.mode(stat) <- "integer"
  stat <- .Call("calcStat",dataset,numCat,homog,
                forbEdges,stat,values,PACKAGE="gRapHD")
  if (dim(stat$errors)[1]!=0)
    warning("Check $errors for edges with problems.", call. = FALSE)
  return(stat)
}
################################################################################

################################################################################
# Plot the graph using Fruchterman-Reingold algorithm.
# Fruchterman, T. M. J., & Reingold, E. M. (1991). Graph Drawing by
# Force-Directed Placement. Software: Practice and Experience, 21(11).
# In: x = gRapHD object
#     vert = list of vertices to be plotted
#     numIter = number of iterations for Fruchterman-Reingold algorithm
#     main = main title
#     plotVert = T or F (add the vertices)
#     energy = T or F (use the minimum energy as initial values)
#     useWeights = T or F (use the model$statSeq as edges lenght weights)
#     vert.hl = list of vertices to highlight
#     col.hl = colour to be used in the highlighted vertices
#     vert.radii = radii of the edges (scalar or vector with length p)
#     coord = matrix (number of vertices by 2) with initial coordenates
#     col.ed = colour of the edges
#     lty.ed = type of line to be used for the edges
#     lwd.ed = width of the edges
#     lwd.vert = width of the vertices' border
#     border = colour of the vertices borders (length 1 or number of vertices)
#     symbol.vert = symbol to be used in the vertices (length 1 or number of
#          vertices) 1 is a circle, 2 a square, 3 or higher represents the
#          number of sides
#     cex.vert.label  = numeric character expansion factor for the labels;
#            multiplied by par yields the final character size. NULL and NA are
#            equivalent to 1.0.
#     vert.labels = labels to be used in the vertices.
#     asp = the aspect ratio y/x (see plot.window for more details).
#     disp = plot (TRUE) or not (FALSE)
#
# Out: matrix (2xp) with the coordinates of the vertices
################################################################################
plot.gRapHD <- function(x,vert=NULL,numIter=50,main="",
                        plotVert=TRUE,plotEdges=TRUE,energy=FALSE,
                        useWeights=FALSE,vert.hl=NULL,col.hl="red",
                        vert.radii=0.01,coord=NULL,col.ed="darkgray",lty.ed=1,
                        lwd.ed=1,lwd.vert=1,border=0,symbol.vert=1,
                        cex.vert.label=.40,vert.labels=TRUE,asp=NA,disp=TRUE,
                        font=par("font"),col.labels=NULL,add=FALSE,...)
{
  model <- x
  rm(x)
  edges <- model$edges
  varType <- model$numCat
  p <- model$p
  originalOrder <- 1:p
  if (is.null(vert))
    vertices <- 1:p
  else
  {
    if ((max(vert) > p) || (min(vert) < 1))
      stop("vert has vertices that are not in the model.")
    if (length(vert)!=length(unique(vert)))
      stop("There are repeated vertices in vert.")
    originalOrder <- order(vert)
    vertices <- sort(vert)
    edges <- SubGraph(edges=edges,v=vertices,p=max(edges))$edges
    edges[,1] <- vertices[edges[,1]]
    edges[,2] <- vertices[edges[,2]]
  }
  p <- length(vertices)

  # vertices radii
  vert.radii <- rep(vert.radii,length.out=p)
  # test if all vertices to be highlighted are in the graph
  if (sum(is.element(vert.hl,vertices))<length(vert.hl))
     stop("Only vertices plotted can be highlighted.")
  # colour of vertices' borders
  border <- rep(border,length.out=p)
  # shape of the vertices
  symbol.vert <- rep(symbol.vert,length.out=p)

  # initialise vertices' coordenates
  if (!is.null(coord)) # if it is given by the user
  {
    if (NROW(coord)!=length(vertices))
      stop("coord must have the same number of rows as vertices.")
    w <- p^2/4
    coord <- w*(2*coord-1) # the algorithm uses the graph centered in (0,0)
  }
  else
    if (!energy) # or using a lattice
    {
      w <- floor(sqrt(p))
      x <- seq(from=-p^2/4,to=p^2/4,by=(2*p^2/4)/w)
      y <- seq(from=p^2/4,to=-p^2/4,by=-(2*p^2/4)/w)
      x <- matrix(rep(x,w+1),nrow=w+1,byrow=T)
      y <- matrix(rep(y,w+1),nrow=w+1)
      x <- t(x)
      y <- t(y)
      coord <- cbind(as.vector(x),as.vector(y))
      coord <- coord[1:p,]
    }
    else # using minimum energy
    {
      if (!useWeights)
      {
        A <- adjMat(edges=edges, p=p)
        Q <- -A
        diag(Q) <- rowSums(A)
      }
      else
      {
        A <- matrix(0,p,nrow(edges))
        A[cbind(edges[,1],1:nrow(edges))] <- -1
        A[cbind(edges[,2],1:nrow(edges))] <- 1
        Q <- A%*%diag(1/abs(model$statSeq))%*%t(A)
      }
      tt <- eigen(Q, symmetric=T)
      v1 <- tt$vectors[ ,p]
      v2 <- tt$vectors[ ,p-1]
      v3 <- tt$vectors[ ,p-2]
      u1 <- v1
      u2 <- v2 - sum(v2*u1)*u1
      u3 <- v3 - sum(v3*u1)*u1 - sum(v3*u2)*u2
      coord <- cbind(u3,u2)+cbind(1:p,1:p)/(1000*p)
    }
  aux <- rep(0,max(vertices))
  aux[vertices] <- 1:length(vertices)
  edgesL <- cbind(aux[edges[,1]],aux[edges[,2]])
  storage.mode(p) <- "integer"
  storage.mode(edgesL) <- "integer"
  storage.mode(numIter) <- "integer"
  storage.mode(coord) <- "double"
  coordV <- .Call("frucRein",edgesL,p,numIter,coord,PACKAGE="gRapHD")
  # rescale coordenates resulting from frucRein
  if (p > 1)
  {
    coordV[,1] <- ((coordV[,1]-min(coordV[,1]))/max(coordV[,1]-min(coordV[,1])))
    coordV[,2] <- ((coordV[,2]-min(coordV[,2]))/max(coordV[,2]-min(coordV[,2])))
    coordV[is.na(coordV)] <- .5
  }
  else
    coordV[1,] <- c(.5,.5*.9)

  if (disp) # if it is to be plotted
  {
    if (is.null(main)) main <- ""
    save.mar <- par("mar")
    save.oma <- par("oma")
    if (!add)
    {
      plot.new()
      par(mar=c(0,0,1.5*(main!=""),0),oma=c(0,0,0,0))
      plot.window(c(0,1),c(0,1),asp=asp)
    }
    arg <- list(...)
    if (!is.null(arg$xpd)) par(xpd=arg$xpd)
    rm(arg)
    coordE <- matrix(0,nrow(edges),4)
    lty.ed <- rep(lty.ed,length.out=length(edgesL))
    lwd.ed <- rep(lwd.ed,length.out=length(edgesL))
    col.ed <- rep(col.ed,length.out=length(edgesL))
    if ((plotEdges) & (length(edgesL)>0))
      for (i in 1:nrow(edgesL))
      {
        coordE[i,1:2] <- c(coordV[edgesL[i,1],1],coordV[edgesL[i,2],1])
        coordE[i,3:4] <- c(coordV[edgesL[i,1],2],coordV[edgesL[i,2],2])
        lines(x=coordE[i,1:2],
              y=coordE[i,3:4],
              col=col.ed[i],lty=lty.ed[i],lwd=lwd.ed[i])#col="darkgray",lty=1,lwd=1)
      }

    if (plotVert)
    {
      xy <- list(x=cos(seq(0,2*pi,2*pi/25)),y=sin(seq(0,2*pi,2*pi/25)))
      symbol.vert <- symbol.vert[originalOrder]
      vert.radii <- vert.radii[originalOrder]
      border <- border[originalOrder]
      Fill <- rep("lightgray",length(varType))
      Fill[varType!=0] <- "black"
      Fill[vert.hl] <- col.hl
      Fill <- Fill[vertices]
      for (j in 1:nrow(coordV))
      {
        if (symbol.vert[j]==0)
          polygon(xy$x*vert.radii[j]*2+coordV[j,1],xy$y*vert.radii[j]+coordV[j,2],border=border[j],col=Fill[j],lwd=lwd.vert)
        else
          if (symbol.vert[j]==1)
            symbols(x=coordV[j,1],y=coordV[j,2],circles=vert.radii[j],inches=FALSE,add=TRUE,bg=Fill[j],
                    fg=border[j],lwd=lwd.vert)
          else
            if (symbol.vert[j]==2)
              symbols(x=coordV[j,1],y=coordV[j,2],squares=vert.radii[j],inches=FALSE,add=TRUE,bg=Fill[j],
                      fg=border[j],lwd=lwd.vert)
            else
              if (symbol.vert[j]>=3)
              {
                z1 <- rep(vert.radii[j],symbol.vert[j])
                dim(z1) <- c(1,symbol.vert[j])
                symbols(x=coordV[j,1],y=coordV[j,2],stars=z1,inches=FALSE,add=TRUE,bg=Fill[j],
                        fg=border[j],lwd=lwd.vert)
              }
      }
    }
    pLabels <- if (vert.labels[1]==F) FALSE else TRUE
    if (pLabels)
    {
      if (is.null(col.labels))
      {
        Fill <- rep("white",length(varType))
        Fill[varType==0] <- "black"
        Fill <- Fill[vertices]
      }
      else
        Fill <- col.labels

      if (!is.logical(vert.labels[1]))
        text(x = coordV[, 1], y = coordV[, 2], vert.labels,
             cex = cex.vert.label, col = Fill, font = font)
      else
        text(x = coordV[, 1], y = coordV[, 2], model$vertNames[vertices],
             cex = cex.vert.label, col = Fill, font = font)    
    }
    if (!is.null(main))
      title(main=main)
    par(mar=save.mar,oma=save.oma,xpd=FALSE)
  }

  invisible(coordV)
}
################################################################################

################################################################################
# Find the sub-graph with all vertices with less or equal distance from a
# central vertex.
# In: model = gRapHD object
#     edges = matrix (n by 2, numeric).
#     orig = central vertex
#     rad = distance
# Out: list - subEdges = matrix (ncol=2)
#             v = matrix (vertex and the distance to the orig)
################################################################################
neighbourhood <- function(model=NULL,edges=NULL,orig=NULL,rad=1)
{
  if (is.null(orig))
    stop("orig must be not NULL.")
    
  if ((is.null(model) && is.null(edges)) || (!is.null(model) && !is.null(edges)))
    stop("Either model or edges must be provided.")
  if (is.null(edges))
  {
    if (class(model)!="gRapHD")
      stop("Model must be of class gRapHD.")
    edges <- model$edges
    p <- model$p
  }
  p <- 0
  if (NROW(edges)==0)
    return(list(edges=matrix(integer(0),,2),v=matrix(c(orig,0),,2)))
  else
    p <- max(p,max(edges))
    
  if (orig > p)
    stop(paste("orig must be in [1,",p,"].",sep=""))
  if (rad < 1)
    vertices <- orig
  else
  {
    vertices <- list()
    vertices[[1]] <- neighbours(edges=edges,v=orig) # first neighbours
    if (rad > 1)
      for (j in 2:rad)
      {
        vertices[[j]] <- integer(0)
        for (i in 1:length(vertices[[j-1]]))
          vertices[[j]] <- c(vertices[[j]],neighbours(edges=edges,v=vertices[[j-1]][i]))
        vertices[[j]] <- unique(setdiff(vertices[[j]],c(orig,unlist(vertices[1:(j-1)]))))
        if (length(vertices[[j]])==0)
        {
          vertices[[j]] <- NULL
          break
        }
      }
  }
  v <- c(orig,0)
  for (i in 1:length(vertices))
    if (length(vertices[[i]])>0)
      v <- rbind(v,cbind(vertices[[i]],i))
  colnames(v) <- c("vertex","rad")
  rownames(v) <- NULL
  V <- sort(unique(v[,1]))
  subEdges <- SubGraph(edges=edges,v=V,p=p)$edges
  subEdges[,1] <- V[subEdges[,1]]
  subEdges[,2] <- V[subEdges[,2]]
  return(list(edges=subEdges,v=v))
}
################################################################################

################################################################################
# Return the adjacency matrix.
# In: model = gRapHD object
#     edges = dataset (n by 2, numeric).
#     p = number of vertices
# Out: matrix (p by p)
################################################################################
adjMat <- function(model=NULL,edges=NULL,p=NULL)
{
  if ((is.null(model) && is.null(edges)) || (!is.null(model) && !is.null(edges)))
    stop("Either model or edges must be provided.")
  if (is.null(edges))
  {
    if (class(model)!="gRapHD")
      stop("Model must be of class gRapHD.")
    edges <- model$edges
    p <- model$p
  }
  if ((NROW(edges)==0) && is.null(p))
    stop("p must not be NULL.")
  if (is.null(p))
    p <- max(edges)

  result <- matrix(0,p,p)
  result[edges] <- 1
  result[edges[,2:1]] <- 1
  return(result)
}
################################################################################

################################################################################
# Return the degree of a list of vertices (or all vertices).
# In: model = gRapHD object
#     edges = list of edges (n by 2, numeric).
#     v = vertices
# Out: vector (length(v))
################################################################################
Degree <- function(model=NULL,edges=NULL,v=NULL)
{
  if ((is.null(model) && is.null(edges)) || (!is.null(model) && !is.null(edges)))
    stop("Either model or edges must be provided.")
  if (is.null(edges))
    if (class(model) != "gRapHD")
      stop("model must be of class gRapHD.")
    else
    {
      p <- model$p
      edges <- model$edges
    }
  else
    p <- max(c(max(edges),v))

  if (is.null(v))
    v <- 1:p
  else
    if (max(v)>p)
      stop(paste("v must be in [1,",p,"].",sep=""))
    
  aux <- sort(as.vector(edges))
  aux <- table(aux)
  result <- rep(0,length(v))
  names(result) <- v
  result[as.character(v)] <- aux[as.character(v)]
  result[is.na(result)] <- 0
  return(result)
}
################################################################################

################################################################################
# Return the shortest paths between v and all other vertices.
# Dijkstra's algorithm
# In: model = gRapHD object
#     edges = dataset (n by 2, numeric).
#     v = vertex
#     p = number of vertices
# Out: vector (max(edges))
################################################################################
shortPath <- function(model=NULL, edges=NULL, v=NULL, p=NULL)
{
  if ((is.null(model) && is.null(edges)) || (!is.null(model) && !is.null(edges)))
    stop("Either model or edges must be provided.")
  if (is.null(edges))
    if (class(model) != "gRapHD")
      stop("model must be of class gRapHD.")
    else
    {
      p <- model$p
      edges <- model$edges
    }
  else
    if (is.null(p))
      p <- max(c(max(edges),v))

  v1 <- edges[,1]
  v2 <- edges[,2]
  storage.mode(v1) <- storage.mode(v2) <- storage.mode(p) <- "integer"
  if (!is.null(v))
  {
    if (max(v)>p)
      stop(paste("v must be in [1,",p,"].",sep=""))
    storage.mode(v) <- "integer"
    result <- .Call("shortPath",v1,v2,v,p,PACKAGE="gRapHD")
    names(result) <- 1:p
  }
  else
  {
    result <- matrix(0,p,p)
    for (v in 1:p)
    {
      storage.mode(v) <- "integer"
      result[v,] <- .Call("shortPath",v1,v2,v,p,PACKAGE="gRapHD")
    }
  }
  result[result==p] <- Inf
  return(result)
}
################################################################################

################################################################################
# Converts the dataset to numeric format
# In: dataset = matrix (numeric) or data frame.
# Out: list = ds (matrix)
#             vertNames (vertices' names)
#             numCat (number of levels in each variable/column)
################################################################################
convData <- function(dataset)
{
  vertNames <- colnames(dataset)
  if (mode(dataset)=="numeric")
  {
    dim(dataset) <- c(NROW(dataset),NCOL(dataset))
    numCat <- rep(0,ncol(dataset))
  }
  else
    if (is.data.frame(dataset))
    {
      numCat <- rep(0,ncol(dataset))
      dsC <- NULL
      for (i in 1:ncol(dataset))
      {
        if (is.numeric(dataset[,i]))
          dsC <- cbind(dsC,dataset[,i])
        else
        {
          if (is.factor(dataset[,i]))
          {
            dsC <- cbind(dsC,as.numeric(dataset[,i]))
            numCat[i] <- nlevels(dataset[,i])
          }
          else
            if (is.logical(dataset[,i]))
            {
              dsC <- cbind(dsC,as.numeric(dataset[,i])+1)
              numCat[i] <- 2
            }
            else
              stop(paste("Column ",i," is not double, integer, factor or logic.",sep=""))
        }
      }
      dataset <- dsC
    }
    else
      stop("Input not numeric nor data frame.")
  if (is.null(vertNames))
    vertNames <- 1:ncol(dataset)
  return(list(ds=dataset,numCat=numCat,vertNames=vertNames))
}
################################################################################

################################################################################
# Calculate a model's log-likelihood, AIC, and BIC
# In: model = gRapHD object
#     edges = dataset (n by 2, numeric).
#     dataset = matrix (nxp)
#     homog = T of F
# Out: number of parameters, -2*logL, BIC, AIC
################################################################################
fit <- function(model=NULL,edges=NULL,dataset,homog=NULL)
{
  if (is.null(edges) && is.null(model))
    stop("Edges or model must be provided.")
  if (!is.null(edges) && !is.null(model))
    stop("Only one (edges or model) must be provided.")
  if (is.null(edges) && (class(model)!="gRapHD"))
    stop("Model must be of class gRapHD.")

  if (is.null(edges))
  {
    edges <- model$edges
    homog <- ifelse(is.null(homog),model$homog,homog)
  }
  else
    homog <- ifelse(is.null(homog),TRUE,homog)

  aux <- convData(dataset)
  dataset <- aux$ds
  varType <- aux$numCat
  varType[varType!=0] <- 1
  numCat <- aux$numCat
  p <- length(numCat)
  n <- nrow(dataset)
  if (sum(is.na(dataset))>0)
    stop("Missing values not allowed.")
  if (!is.null(model))
  {
    if (!all.equal(numCat,model$numCat))
      stop("model$numCat doesn't agree with data.")
    if (!identical(sort(aux$vertNames),sort(model$vertNames)))
      stop("model$vertNames doesn't agree with data.")
    dataset <- dataset[,match(model$vertNames,aux$vertNames)]
  }

  if (length(unique(varType))==1) #all discrete or all continuous
    perfect <- perfSets(edges=edges,p=p,varType=varType,from=0)
  else
  {
    # star graph
    perfect <- perfSets(edges=rbind(edges,cbind(which(varType==1),p+1)),p=p+1,varType=c(varType,1),from=0)
    if (!is.list(perfect))
      stop("The model is not strongly decomposable.")
    perfect <- perfSets(edges=edges,p=p,varType=varType,from=0)
  }
  if (!is.list(perfect))
    stop("The model is not triangulated.")

  L <- 0
  numP <- 0

  whichSeq <- function(x,y)
  {
    yes <- TRUE
    i <- 0
    while ((yes) & (i<length(x)))
      yes <- x[i<-i+1]==y[i]
    return(yes)
  }

  for (i in 1:length(perfect$cliques))
  {
    x <- unlist(perfect$cliques[i])
    discrete <- x[which(varType[x]==1)]
    continuous <- setdiff(x,discrete)
    ds <- dataset[,continuous]
    dim(ds) <- c(NROW(ds),NCOL(ds))

    if ((length(continuous)>0) & (length(discrete)==0))
    {
      m <- colMeans(ds)
      covm <- var(ds)*(n-1)/n
      L <- L + sum(normDens(ds,m,covm,log=TRUE))
      numP <- numP + length(continuous)*(length(continuous)+3)/2
    }
    else
      if ((length(continuous)==0) & (length(discrete)>0))
      {
        tc <- table(as.data.frame(dataset[,discrete]))
        tc <- tc[tc>0]
        L <- L + sum(tc*log(tc/sum(tc)))
        numP <- numP + length(tc)
      }
      else
        if (homog)
        {
          tc <- table(as.data.frame(dataset[,discrete]))
          covm <- matrix(0,NCOL(ds),NCOL(ds))
          numUsedObs <- 0
          for (j in 1:length(tc))
            if (tc[j]>length(continuous))
            {
              ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,numCat[discrete]))
              covm <- covm + var(ds[ind,])*(tc[j]-1)
              numUsedObs <- numUsedObs + tc[j]
            }
          covm <- covm/sum(tc)
          for (j in 1:length(tc))
            if (tc[j]>0)
            {
              ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,numCat[discrete]))
              m <- colMeans(matrix(ds[ind,],,length(continuous)))
              L <- L + sum(normDens(ds[ind,],m,covm,log=TRUE)) + tc[j]*log(tc[j]/sum(tc))
              numP <- numP + 1 + length(continuous)
            }
          numP <- numP + length(continuous)*(length(continuous)+1)/2
        }
        else
        {
          tc <- table(as.data.frame(dataset[,discrete]))
          covm <- matrix(0,NCOL(ds),NCOL(ds))
          numUsedObs <- 0
          for (j in 1:length(tc))
            if (tc[j]>length(continuous))
            {
              ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,numCat[discrete]))
              covm <- var(ds[ind,])*(tc[j]-1)/tc[j]
              m <- colMeans(matrix(ds[ind,],,length(continuous)))
              L <- L + sum(normDens(ds[ind,],m,covm,log=TRUE)) + tc[j]*log(tc[j]/sum(tc))
              numP <- numP + length(continuous)*(length(continuous)+3)/2 + 1
            }
        }

    x <- unlist(perfect$separators[i])
    if (length(x)>0)
    {
      discrete <- x[which(varType[x]==1)]
      continuous <- setdiff(x,discrete)
      ds <- dataset[,continuous]
      dim(ds) <- c(NROW(ds),NCOL(ds))

      if ((length(continuous)>0) & (length(discrete)==0))
      {
        m <- colMeans(ds)
        covm <- var(ds)*(n-1)/n
        L <- L - sum(normDens(ds,m,covm,log=TRUE))
        numP <- numP - length(continuous)*(length(continuous)+3)/2
      }
      else
        if ((length(continuous)==0) & (length(discrete)>0))
        {
          tc <- table(as.data.frame(dataset[,discrete]))
          tc <- tc[tc>0]
          L <- L - sum(tc*log(tc/sum(tc)))
          numP <- numP - length(tc)
        }
        else
          if (homog)
          {
            tc <- table(as.data.frame(dataset[,discrete]))
            covm <- matrix(0,NCOL(ds),NCOL(ds))
            numUsedObs <- 0
            for (j in 1:length(tc))
              if (tc[j]>length(continuous))
              {
                ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,numCat[discrete]))
                covm <- covm + var(ds[ind,])*(tc[j]-1)
                numUsedObs <- numUsedObs + tc[j]
              }
            covm <- covm/sum(tc)
            for (j in 1:length(tc))
              if (tc[j]>0)
              {
                ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,numCat[discrete]))
                m <- colMeans(matrix(ds[ind,],,length(continuous)))
                L <- L - sum(normDens(ds[ind,],m,covm,log=TRUE)) - tc[j]*log(tc[j]/sum(tc))
                numP <- numP - 1 - length(continuous)
              }
            numP <- numP - length(continuous)*(length(continuous)+1)/2
          }
          else
          {
            tc <- table(as.data.frame(dataset[,discrete]))
            covm <- matrix(0,NCOL(ds),NCOL(ds))
            numUsedObs <- 0
            for (j in 1:length(tc))
              if (tc[j]>length(continuous))
              {
                ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,numCat[discrete]))
                covm <- var(ds[ind,])*(tc[j]-1)/tc[j]
                m <- colMeans(matrix(ds[ind,],,length(continuous)))
                L <- L - sum(normDens(ds[ind,],m,covm,log=TRUE)) - tc[j]*log(tc[j]/sum(tc))
                numP <- numP - length(continuous)*(length(continuous)+3)/2 - 1
              }
          }
    }
  }

  if (!is.element(sum(varType)/p,c(1,0)))
  {
    if (is.null(model))
    {
      model <- list()
      model$p <- p
      model$edges <- edges
      model$varType <- varType
      model$numCat <- numCat
      class(model) <- "gRapHD"
    }
    model$homog <- homog
    numP <- modelDim(model)
  }
  result <- c(numP,-2*L,-2*L+2*numP,-2*L+numP*log(nrow(dataset)))
  names(result) <- c("Number of parameters","-2*Log-likelihood","AIC","BIC")
  return(result)
}
################################################################################

################################################################################
# To be used as useMethod by the summary function.
# In: gRapHD object
################################################################################
summary.gRapHD <- function(object,...)
{
  if (!is.gRapHD(object))
    stop("Not a gRapHD object.")

  cat(paste("Number of edges       = ",nrow(object$edges),"\n",sep=""))
  cat(paste("Number of vertices    = ",object$p,"\n",sep=""))
  p <- max(c(1,object$p))
  aux <- ifelse(sum(object$numCat)==0,"continuous",ifelse(sum(object$numCat!=0)/p==1,"discrete","mixed"))
  aux1 <- ifelse(aux=="mixed",paste(" and ",ifelse(object$homog,"homogeneous","heterogeneous"),sep=""),"")
  cat(paste("Model                 = ",aux,aux1," \n",sep=""))
  if (!is.null(object$stat.minForest))
    cat(paste("Statistic (minForest) = ",object$stat.minForest,"\n",sep=""))
  if (!is.null(object$stat.stepw))
    cat(paste("Statistic (stepw)     = ",object$stat.stepw,"\n",sep=""))
  if (!is.null(object$minForest))
    cat(paste("Edges (minForest)     = ",object$minForest[1],"...",object$minForest[2],"\n",sep=""))
  if (length(object$stepw)>0)
    cat(paste("Edges (stepw)         = ",object$stepw[1],"...",object$stepw[2],"\n",sep=""))
  else
    if (!is.null(object$stepw))
      cat(paste("Edges from stepw      = ","no edge added","\n",sep=""))
  cat("\n")
}
################################################################################

################################################################################
# see pages (Steffen's): 202-203; 215-216
# In: model - gRapHD object
# Out: list - discrete: (d,\empty)
#             linear: (d,\gamma^2)
#             quadratic: (d,\gamma) and (d,{\gamma,\mu})
#             quadratic2: : (d,c^2)
################################################################################
modelFormula <- function(model)
{
  if (class(model)!="gRapHD")
    stop("Model must be of class gRapHD.")

  psSub <- function(edges,p,v,varType)
  {
    n <- length(v)
    edgesV <- SubGraph(edges=edges,v=v,p=model$p)$edges
    pSeq <- perfSets(edges=edgesV,p=n,varType=varType,from=0)
    for (i in 1:length(pSeq$cliques))
      pSeq$cliques[[i]] <- v[pSeq$cliques[[i]]]
    return(pSeq$cliques)
  }

  edges <- model$edges
  varType <- model$numCat
  varType[varType!=0] <- 1
  p <- model$p

  discr <- which(varType==1)
  nDiscr <- length(discr)
  discrete <- list()
  # find all d terms
  if (nDiscr>0)
    discrete <- psSub(edges,p,discr,varType[discr])

  ind <- cbind(varType[edges[,1]],varType[edges[,2]])
  # using this cont, only those continuous that have at least one edge with
  # a discrete are considered
  cont <- sort(unique(c(edges[which((ind[,1]==0) & (ind[,2]==1)),1],
                        edges[which((ind[,2]==0) & (ind[,1]==1)),2])))
  otherCont <- sort(setdiff(which(varType==0),cont))
  # and like this, consider all continuous
  #cont <- which(varType==0)

  # find all terms (d,gamma^2)
  linear <- list()
  if (length(cont)>0)
    for (i in 1:length(cont))
    {
      aux <- list()
      pSeq <- psSub(edges,model$p,c(discr,cont[i]),c(rep(1,nDiscr),0))
      for (j in 1:length(pSeq))
        if (is.element(cont[i],pSeq[[j]]))
          aux <- c(aux,list(sort(pSeq[[j]])))
      linear <- c(linear,aux) #linear <- c(linear,list(aux))
    }

  # reduction of discrete terms
  if (length(linear)>0)
    for (i in 1:length(linear))
    {
      j <- 1
      while (j<=length(discrete))
      {
        if (length(setdiff(discrete[[j]],linear[[i]]))==0)
          discrete <- discrete[-j]
        else
          j <- j + 1
      }
    }

  if (model$homog)
  {
    quadratic <- linear
    linearX <- list()
    rm(linear)
    x <- sort(c(otherCont,cont))
    pSeq <- psSub(edges,model$p,x,varType[x])
    quadratic2 <- pSeq
  }
  else
  {
    # find all terms (d,{gamma,mu})
    quadratic <- list()
    if (length(cont)>0)
      for (i in 1:length(cont))
      {
        for (j in (i+1):length(cont))
          if ((i!=j) & (j<=length(cont)))
          {
            aux <- list()
            pSeq <- psSub(edges,model$p,c(discr,cont[c(i,j)]),
                                         c(rep(1,nDiscr),0,0))
            for (k in 1:length(pSeq))
              if ((is.element(cont[i],pSeq[[k]]))&(is.element(cont[j],pSeq[[k]])))
                aux <- c(aux,list(sort(pSeq[[k]])))#list(sort(pSeq[[k]])))
            quadratic <- c(quadratic,aux) #quadratic <- c(quadratic,list(aux))
          }
      }

    # find the indexes in linear and quadratic that have the same d
    # so, in position i of dLinead are all indexes of linear that have the same
    # d; idem for dQuadratic
    dLinear <- list()
    dQuadratic <- list()
    visitedL <- rep(FALSE,length(linear))
    visitedQ <- rep(FALSE,length(quadratic))
    v <- 1
    storage.mode(v) <- "integer"
    k <- 1
    while (sum(visitedL==0)>0)
    {
      d <- intersect(discr,linear[[v]])
      visitedL[v] <- TRUE
      dLinear[[k]] <- v
      for (i in 1:length(linear))
      {
        if (!visitedL[i])
          if ((length(setdiff(d,intersect(discr,linear[[i]])))==0) &&
              (length(setdiff(intersect(discr,linear[[i]]),d))==0))
          {
            dLinear[[k]] <- c(dLinear[[k]],i)
            visitedL[i] <- TRUE
          }
      }
      if (length(quadratic)>0)
      {
        k1 <- NULL
        for (i in 1:length(quadratic))
        {
          if (!visitedQ[i])
            if ((length(setdiff(d,intersect(discr,quadratic[[i]])))==0) &&
                (length(setdiff(intersect(discr,quadratic[[i]]),d))==0))
            {
              k1 <- c(k1,i)
              visitedQ[i] <- TRUE
            }
        }
        if (length(k1)>0)
          dQuadratic[[length(dQuadratic)+1]] <- k1
        else
          dQuadratic[[length(dQuadratic)+1]] <- integer(0)
      }
      k <- k + 1
      v <- which(visitedL==0)[1]
    }

    # remove redundant generators - for those which have a common element d for a
    # range of pairs \gamma, \mu \in c \subset \Gamma, the corresponding list
    # of generators are replaced for (d,c^2); but the way it's here, goes a bit
    # further than what is described in Steffen's book, getting some like what
    # is described in example 6.29 ("Here we have taken advantage of the fact...")
    quadratic2 <- list()
    linearX <- list()
    if (length(quadratic)>0)
    {
      for (i in 1:length(dLinear))
      {
        x <- sort(unique(c(unlist(linear[dLinear[[i]]]),
                           unlist(quadratic[dQuadratic[[i]]]))))
        pSeq <- psSub(edges,model$p,x,varType[x])
        for (j in 1:length(pSeq))
          if (length(intersect(cont,pSeq[[j]]))==1)
            linearX <- c(linearX,pSeq[j])
          else
          {
            quadratic2 <- c(quadratic2,pSeq[j])
            for (k in 1:length(dQuadratic[[i]]))
              quadratic[[dQuadratic[[i]][k]]] <- integer(0)
          }
      }
      for (i in 1:length(quadratic))
        if (length(quadratic[[i]])>0)
          quadratic2 <- c(quadratic2,quadratic[i])
    }
    else
      linearX <- linear
    rm(linear)
    # remove all empty elements of quadratic
    if (length(quadratic) > 0)
    {
      aux <- quadratic
      quadratic <- list()
      for (i in 1:length(aux))
        if (length(aux[[i]])>0)
        quadratic <- c(quadratic,aux[i])
    }
    # get the generators with only continuous variables
    x <- sort(c(otherCont,cont))
    pSeq <- psSub(edges,model$p,x,varType[x])
    for (i in 1:length(pSeq))
      if (length(intersect(pSeq[[i]],otherCont))>0)
        if (length(pSeq[[i]])==1)
          linearX <- c(linearX,pSeq[i])
        else
          quadratic2 <- c(quadratic2,pSeq[i])
  }
  return(list(discrete=discrete,
              linear=linearX,
              quadratic=quadratic,
              quadratic2=quadratic2))
}
################################################################################

################################################################################
# Model's dimension
# In: model - gRapHD object
# Out: number of parameters
################################################################################
modelDim <- function(model)
{
  if (class(model)!="gRapHD")
    stop("Model must be of class gRapHD.")

  mf <- modelFormula(model)
  exf <- c(rep(0,length(mf$discrete)),
           rep(2,length(mf$linear)),
           rep(1,length(mf$quadratic)),
           rep(2,length(mf$quadratic2)))
  mf <- c(mf$discrete,mf$linear,mf$quadratic,mf$quadratic2)
  for (i in 1:length(mf))
    mf[i] <- list(as.integer(mf[[i]]))

  numCat <- as.integer(model$numCat)
  storage.mode(exf) <- "integer"
  result <- .Call("modelDim",mf,exf,numCat,PACKAGE="gRapHD")
  return(result)
}
################################################################################

################################################################################
# To be used as useMethod by the as function.
# In: at least a list of edges
# Out: gRapHD object
################################################################################
as.gRapHD <- function(object,...)
{
  arg <- list(...)
  result <- list()
  n <- 0 # number of edges
  p <- 0 # number of vertices
  if (is.null(object))
    object <- matrix(integer(0),,2) # empty list of edges
  else
    if (NCOL(object)==2) # matrix of edges
    {
      n <- NROW(object)
      p <- ifelse(n==0,0,max(object))
      if (n>0)
      {
        a <- object[,1]-object[,2] 
        object[a>0,] <- cbind(object[a>0,2],object[a>0,1]) # sort rows
        if (0%in%a)
        {
          warning("Self-loops removed.")
          x <- x[-which(a==0),]
        }
        z0 <- (object[,1]-1)*p-(object[,1]-1)*object[,1]/2+object[,2]-object[,1]
        z1 <- unique(z0)
        if (length(z1)<n) # there are multiple edges
        {
          warning("Multiple edges removed.")
          object <- object[match(z1,z0),]
        }
      }
      n <- min(c(n,1))
    }
    else
      stop("Object must be NULL, or with dimension (k,2).")
  result$edges <- object
  result$p <- ifelse(!is.null(arg$p),as.integer(arg$p),as.integer(p))
  result$stat.user <- ifelse(!is.null(arg$stat),as.character(arg$stat),"LR")
  result$statSeq <- if(!is.null(arg$statSeq)) as.numeric(arg$statSeq) else as.numeric(rep(NA,nrow(result$edges)))
  result$numCat <- if(!is.null(arg$numCat)) as.integer(arg$numCat) else integer(result$p)
  result$homog <- ifelse(!is.null(arg$homog),as.logical(arg$homog),TRUE)
  result$numP <- if(!is.null(arg$numP)) as.numeric(arg$numP) else as.numeric(rep(NA,nrow(result$edges)))
  if (!is.null(arg$vertNames))
    result$vertNames <- arg$vertNames
  else
    if (result$p == 0)
      result$vertNames <- NA
    else
      result$vertNames <- 1:result$p
  result$userDef <- c(n,nrow(result$edges))
  class(result)<-"gRapHD"
  return(result)
}
################################################################################

################################################################################
# To be used as useMethod by the print function.
# In: gRapHD object
################################################################################
print.gRapHD <- function(x,...)
{
  cat("gRapHD object\n")
  summary(x)
}
################################################################################

################################################################################
# To be used as useMethod by the is function.
# In: gRapHD object
# Out: T or F
################################################################################
is.gRapHD <- function(object)
{
  class(object) == "gRapHD"
}
################################################################################

################################################################################
# Finds a junction tree.
# In: gRapHD object
# Out: separators - list with unique minimal separators.
#      juncTree - edges in the tree (each vertex is a clique in the list below).
#      sepSubSetOfSep - list in which each element gives all the separators
#                       which contain this respective separator.
#      indSepOrig - index of the original separator (in the MCS result) in the
#                   list above.
#      cliques - list with cliques.
################################################################################
jTree <- function(model)
{
  if (!is.gRapHD(model))
    stop("Object is not of gRapHD class.")
    
  v1 <- model$edges[,1]
  v2 <- model$edges[,2]
  p <- model$p
  varType <- model$numCat
  varType[varType!=0] <- 1
  storage.mode(v1) <- storage.mode(v2) <- "integer"
  storage.mode(p) <- "integer"
  storage.mode(varType) <- "integer"

  result <- .Call("juncTree",v1,v2,p,varType,PACKAGE="gRapHD")

  return(result)
}
################################################################################

################################################################################
# Calculates the clustering coefficient for the given graph
# In: model = gRapHD object
#     edges = dataset (n by 2, numeric)
#     p = number of vertices
# Output: C
################################################################################
ccoeff <- function(model=NULL,edges=NULL,p=NULL)
{
  if ((is.null(model) && is.null(edges)) || (!is.null(model) && !is.null(edges)))
    stop("Either model or edges must be provided.")
  if (is.null(edges))
  {
    if (class(model)!="gRapHD")
      stop("Model must be of class gRapHD.")
    edges <- model$edges
    p <- model$p
  }
  if (NROW(edges)==0)
    return(0)
  if (is.null(p))
    p <- max(edges)
  adjM <- adjMat(edges=edges,p=p)
  # select all nodes with edge to node i (all rows and columns connected to i)
  cc <- rep(0,p)
  for (i in 1:p)
  {
    j <- which(adjM[i,]==1)
    if (length(j) > 1)
    {
      aux <- adjM[j,j]
      Ei <- sum(aux)/2
      k <- nrow(aux)
      cc[i] <- 2*Ei/(k*(k-1))
    }
  }
  return(cc)
}
################################################################################

################################################################################
# Calculates the deviance for the inclusion of a new edge.
# In: x = vertex 1
#     y = vertex 2
#     S = separator
#     dataset = dataframe
#     homog = homogeneous (T) or not (F)
# Out: deviance and number of parameters
################################################################################
CI.test <- function(x,y,S,dataset,homog=TRUE)
{
  edge <- sort(c(x,y))
  p <- ncol(dataset) # number of variables (vertices)
  n <- nrow(dataset) # number of observations
  result <- numP <- NA

  numbPar <- function(y)
  {
    a <- sum(rowSums(y)!=0)-1
    b <- sum(colSums(y)!=0)-1
    return(a*b*(a>0))
  }
  whichSeq <- function(x,y){return(sum(abs(x-y))==0)}
  Gcliq <- function(v)
  {
    v <- sort(unique(v))
    i <- 1
    n <- length(v)
    ed <- sort(rep(i:(i+n-1),n))
    ed <- cbind(ed,rep(i:(i+n-1),n))
    ed <- matrix(ed[-unlist(lapply(1:n,function(x,n){1:x+(x-1)*n},n=n)),],,2)
    dimnames(ed) <- NULL
    ed[,1] <- v[ed[,1]]
    ed[,2] <- v[ed[,2]]
    return(ed)
  }

  dataset <- convData(dataset)
  numCat <- dataset$numCat
  dataset <- dataset$ds
  if (sum(is.na(dataset))>0)
    stop("Missing values not allowed.")

  # this is the clique that will result from the inclusion of this edge
  clique <- c(edge,S)

  # the 2 cliques thal will be merged
  clique1 <- c(edge[1],S)
  clique2 <- c(edge[2],S)

  if (sum(numCat[clique]) == 0) # continuous
  {
    CM <- cov(dataset[,clique],use="pairwise.complete.obs")*(n-1)
    if (sum(is.na(CM))) # meaning that the new edge does not have more than 1 observation (not NA)
      result <- NA
    else
    {
      a <- match(edge,clique)
      d <- log(det(CM)) + ifelse(length(S)==0,0,log(det(matrix(CM[-a,-a],length(clique)-2))))
      d <- d - (log(det(matrix(CM[-a[1],-a[1]],length(clique)-1))) +
                log(det(matrix(CM[-a[2],-a[2]],length(clique)-1))))
      result <- nrow(dataset)*d
      numP <- 1
    }
  }
  else
    if (sum(numCat[clique]!=0)==length(clique)) # discrete
    {
      t12 <- table(as.data.frame(dataset[,clique]))
      t1 <- margin.table(t12,match(clique1,clique))
      t2 <- margin.table(t12,match(clique2,clique))
      tS <- if(length(S)==0) 0 else margin.table(t12,match(S,clique))
      if (length(S)==0)
        numP <- prod(numCat[clique]-1)#length(t12)-1
      else
        numP <- sum(apply(t12,MARGIN=(1:length(clique))[-(1:2)],numbPar))
      numObsUsed <- sum(t12)
      t12 <- t12[t12!=0]
      t1 <- t1[t1!=0]
      t2 <- t2[t2!=0]
      tS <- tS[tS!=0]

      if (numObsUsed>0)
      {
        v12 <- sum(t12*log(t12/numObsUsed))
        v1  <- sum(t1 *log(t1 /numObsUsed))
        v2  <- sum(t2 *log(t2 /numObsUsed))
        vS  <- ifelse(length(S)==0,0,sum(tS *log(tS /numObsUsed)))
        result <- -2*(v12 + vS - v1 - v2)
      }
    }
    else # mixed
    {
      discrete   <- clique[numCat[clique]!=0] #intersect(clique,which(numCat!=0))
      continuous <- clique[numCat[clique]==0] #setdiff(clique,discrete)
      tab <- table(as.data.frame(dataset[,discrete]))
      ssd <- diag(0,length(continuous))
      if (sum(numCat[edge]==0)==2) # both are continuous
        if (homog)
        {
          for (j in 1:length(tab))
            if (tab[j]>0)
            {
              ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,numCat[discrete]))
              if (sum(ind)>1) #otherwise there is only one observation, and the contribution to the variance is zero
                ssd <- ssd + cov(dataset[ind,continuous],use="pairwise.complete.obs")*(tab[j]-1)
            }
          ind <- match(edge,continuous)
          result <- n*(log(det(ssd)) +
                    ifelse(length(continuous)>2,log(det(matrix(ssd[-ind,-ind],length(continuous)-2))),0) -
                    log(det(matrix(ssd[-ind[1],-ind[1]],length(continuous)-1))) -
                    log(det(matrix(ssd[-ind[2],-ind[2]],length(continuous)-1))))
          numP <- 1
        }
        else
        {
          result <- 0
          ind1 <- match(edge,continuous)
          for (j in 1:length(tab))
            if (tab[j]>0)
            {
              ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,numCat[discrete]))
              if (sum(ind)>1) #otherwise there is only one observation, and the contribution to the variance is zero
                ssd <- cov(dataset[ind,continuous],use="pairwise.complete.obs")*(tab[j]-1)
              result <- result + tab[j]*(log(det(ssd)) +
                        ifelse(length(continuous)>2,log(det(matrix(ssd[-ind1,-ind1],length(continuous)-2))),0) -
                        log(det(matrix(ssd[-ind1[1],-ind1[1]],length(continuous)-1))) -
                        log(det(matrix(ssd[-ind1[2],-ind1[2]],length(continuous)-1))))
            }
          numP <- sum(tab>(length(continuous)+1))
        }
      else
      {
        vDiscr <- edge[numCat[edge]!=0]
        vCont  <- setdiff(edge,vDiscr)
        discrMarg <- setdiff(discrete,vDiscr)
        tabMarg <- margin.table(tab,match(discrMarg,discrete))
        ssdMarg <- diag(0,length(continuous))
        if (homog)
        {
          for (j in 1:length(tab))
          {
            if (tab[j]>0)
            {
              ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,numCat[discrete]))
              if (sum(ind)>1) #otherwise there is only one observation, and the contribution to the variance is zero
                ssd <- ssd + var(dataset[ind,continuous],use="pairwise.complete.obs")*(tab[j]-1)
            }
            if (j<=length(tabMarg))
            {
              if (tabMarg[j]>0)
              {
                if (length(discrMarg)>0)
                  ind <- apply(matrix(dataset[,discrMarg],,length(discrMarg)),1,whichSeq,y=seqLevels(j,numCat[discrMarg]))
                else
                  ind <- rep(TRUE,tabMarg[j])
                if (sum(ind)>1) #otherwise there is only one observation, and the contribution to the variance is zero
                  ssdMarg <- ssdMarg + var(dataset[ind,continuous],use="pairwise.complete.obs")*(tabMarg[j]-1)
              }
            }
          }
          ind <- match(vCont,continuous)
          result <- n*(log(det(ssd)) +
                    ifelse(length(continuous)>1,log(det(matrix(ssdMarg[-ind,-ind],length(continuous)-1))),0) -
                    ifelse(length(continuous)>1,log(det(matrix(ssd[-ind,-ind],length(continuous)-1))),0) -
                    log(det(ssdMarg)))
          if (length(discrete)==1) #means that the separator is continuous
            numP <- numCat[vDiscr]-1
          else
          {
            aux <- match(vDiscr,discrete)
            # dimension 1 in aux1 is the discrete variable in the tested edge
            aux1 <- table(as.data.frame(dataset[,c(vDiscr,discrete[-aux])]))
            aux.tab <- margin.table(aux1,2:length(discrete))
            aux.ind <- which(aux.tab>0,arr.ind=TRUE)
            numP <- 0
            for (i in 1:nrow(aux.ind))
            {
              aux2 <- (1:numCat[vDiscr]) + sum((aux.ind[i,]-1)*cumprod(numCat[discrete[-aux]]))
              numP <- numP + (sum(aux1[aux2]>0) - 1)
            }
            rm(aux,aux1,aux2,aux.tab,aux.ind)
          }
        }
        else
        {
          result <- 0
          ind1 <- match(vCont,continuous)
          for (j in 1:length(tab))
          {
            if (tab[j]>0)
            {
              ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,numCat[discrete]))
              ssd <- var(dataset[ind,continuous],use="pairwise.complete.obs")*(tab[j]-1)
              dim(ssd) <- rep(length(continuous),2)
            }
            if (j<=length(tabMarg))
            {
              if (tabMarg[j]>0)
              {
                if (length(discrMarg)>0)
                  ind <- apply(matrix(dataset[,discrMarg],,length(discrMarg)),1,whichSeq,y=seqLevels(j,numCat[discrMarg]))
                else
                  ind <- rep(TRUE,tabMarg[j])
                ssdMarg <- var(dataset[ind,continuous],use="pairwise.complete.obs")*(tabMarg[j]-1)
                dim(ssdMarg) <- rep(length(continuous),2)
              }
            }
            ind <- match(vCont,continuous)
            result <- result +
                      -tab[j]*(log(det(ssd))-ifelse(length(continuous)>1,log(det(matrix(ssd[-ind,-ind],length(continuous)-1))),0)) +
                      +ifelse(j<=length(tabMarg),tabMarg[j],0)*(log(det(ssdMarg))-ifelse(length(continuous)>1,log(det(matrix(ssdMarg[-ind,-ind],length(continuous)-1))),0))
          }
          result <- result +
                    sum(tab*log(tab),na.rm=TRUE) -
                    sum(tabMarg*log(tabMarg),na.rm=TRUE)
          result <- -result
          if (length(discrete)==1) #means that the separator is continuous
            numP <- sum(tab>(length(continuous)+1))-1
          else
          {
            aux <- match(vDiscr,discrete)
            # dimension 1 in aux1 is the discrete variable in the tested edge
            aux1 <- table(as.data.frame(dataset[,c(vDiscr,discrete[-aux])]))
            aux.tab <- margin.table(aux1,2:length(discrete))
            aux.ind <- which(aux.tab>(length(continuous)+1),arr.ind=TRUE)
            numP <- 0
            for (i in 1:nrow(aux.ind))
            {
              aux2 <- (1:numCat[vDiscr]) + sum((aux.ind[i,]-1)*cumprod(numCat[discrete[-aux]]))
              numP <- numP + (sum(aux1[aux2]>0) - 1)
            }
            rm(aux,aux1,aux.tab,aux.ind)
          }
        }
      }
    }

  return(list(deviance=-result,numP=numP))
}
################################################################################

################################################################################
# Given a "position" in a list of combination of variables levels, return the
# respective combination.
# E.g., if numCat=c(2,3,2)
#   comb = dim1 dim2 dim3
#             1    1    1
#             2    1    1
#             1    2    1
#             2    2    1
#             1    3    1
#             2    3    1
#             1    1    2
#             2    1    2
#             1    2    2
#             2    2    2
#             1    3    2
#             2    3    2
#  if k=7 then the result is c(1,1,2)
################################################################################
seqLevels <- function(k,numCat)
{
  n <- length(numCat)
  x <- rep(0,n)
  if (n>0)
    for (i in 1:n)
    {
      y <- k%%numCat[i]
      x[i] <- ifelse(y==0,numCat[i],y)
      k <- (k-y)/numCat[i] + ifelse(y>0,1,0)
    }
  return(x)
}
################################################################################

################################################################################
# Calculates the normal (uni or multivariate) density.
# In: x = vector or matrix of values
#     meanVec = mean vector (or scalar).
#     covMat = covariance matrix (or scalar)
#     logx = if to return the log or not
# Out: density (or log) values
################################################################################
normDens <- function (x, meanVec, covMat, logx = FALSE)
{
    if (NCOL(covMat)==1)
      result <- dnorm(x,mean=meanVec,sd=sqrt(covMat),log=logx)
    else
    {
      distM <- mahalanobis(x,center=meanVec,cov=covMat)
      logD <- sum(log(eigen(covMat,symmetric=TRUE,only.values=TRUE)$values))
      logVal <- -(ncol(x)*log(2*pi)+logD+distM)/2
      if (logx)
        result <- logVal
      else
        result <- exp(logVal)
    }
  return(result)
}
################################################################################

################################################################################
# Converts objects between gRapHD and graphNEL classes.
# In: object = gRapHD or graphNEL object
# Out: object with class oposite to the one in the input
################################################################################
convertClass <- function(object)
{
  if (as.character(class(object))=="graphNEL")
  {
    v <- object@nodes
    p <- length(object@nodes)
    edges <- NULL
    for (i in 1:p)
    {
      x <- object@edgeL[[v[i]]]$edges
      if (length(x)>0)
        edges <- rbind(edges,cbind(i,x))
    }
    edges <- t(apply(edges,1,sort))
    edges <- unique(edges,MARGIN=1)

    g1 <- as.gRapHD(edges,p=p)
  }
  else
    if (class(object)=="gRapHD")
    {
      require(graph)
      edgeL <- vector("list",length=object$p)
      names(edgeL) <- 1:object$p
      I <- adjMat(edges=object$edges,p=object$p)
      for (i in 1:object$p)
        edgeL[[i]] <- list(edges=(1:object$p)[I[i,]==1],weights=rep(1,sum(I[i,])))
      g1 <- new("graphNEL", nodes=as.character(1:object$p), edgeL=edgeL, edgemode = "undirected")
    }
    else
      stop("object must be of class graphNEL or gRapHD.")
  return(g1)
}
################################################################################