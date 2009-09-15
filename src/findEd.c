/******************************************************************************/
//  This file is part of gRapHD R package.
//
//  gRapHD R package
//  Copyright (C) 2009 Gabriel Coelho Goncalves de Abreu, Rodrigo Labouriau,
//  and David Edwards
//
//  gRapHD R package is free software: you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or any later version.
//
//  gRapHD R package program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
//  Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with this program.  If not, see <http://www.gnu.org/licenses/>.
//
/******************************************************************************/

#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdbool.h>
#include <R.h>
#include <Rdefines.h>
#include "comfunc.h"

/******************************************************************************/
// structure used in a linked list containing (per node) possible edges, the
// index of each edge in the previous edges (i.e., if it was already
// considered in a previous iteration), and the separator for these edges.
// There is a index to the next node.
/******************************************************************************/
struct nodeFE
{
  unsigned int *vert1; //first vertex in the edge
  unsigned int *vert2; //second vertex in the edge
  unsigned int *index; //index of the edge in a previous list
  SEXP S; //separator for all edges in this node
  struct nodeFE *next; //pointer to the next node
};

/******************************************************************************/
// Creates a new node, with the parameters.
// In: vert1 - pointer to int vector (vert1[0] = number of items in the vector)
//             containing the first vertex of the edge
//     vert2 - idem, for the second vertex
//     index - pointer to int vector (idem) with the index for that edge in the
//             list of previous edges (if is a new edge, index[i]=0)
//     S - pointer SEXP with the separator
// Out: pointer to the new node
/******************************************************************************/
struct nodeFE* newNodeFE(unsigned int *vert1, unsigned int *vert2, unsigned int *index, SEXP S)
{
  struct nodeFE* nodeFE = malloc(sizeof *nodeFE); //new node
  nodeFE->vert1 = vert1;
  nodeFE->vert2 = vert2;
  nodeFE->index = index;
  nodeFE->S = S;
  nodeFE->next = NULL; //there is no next
  return(nodeFE);
}

/******************************************************************************/
// Like R's "which.max", returns the index for the maximum value in a vector.
// In: d - pointer to int vector
// Out: int*, pointer to an array with the indexes for the maximum values in d.
/******************************************************************************/
unsigned int *whichMaxFE(unsigned int *d)
{
  unsigned int i;
  unsigned int m;
  unsigned int *v, *result;

  v = (unsigned int *)malloc((d[0]+1)*sizeof(unsigned int));
  m = d[1]; //set the maximum as the first element
  v[0] = 1; //set the number of elements in the result as 1
  v[1] = 1; //the index of the first "maximum"
  for(i=2; i<=d[0]; i++)
  {
    if (m < d[i])
    {
      m = d[i];
      v[0] = 1;
      v[1] = i;
    }
    else
      if (m == d[i])
      {
        v[0]++;
        v[v[0]] = i;
      }
  }
  
  result = (unsigned int *)malloc((v[0]+1)*sizeof(unsigned int));
  for (i=0; i<=v[0]; i++)
    result[i] = v[i];
    
  free(v);
    
  return(result);
}

/******************************************************************************/
// Maximum Cardinality Search (Tarjan & Yannakakis, in Jayson Rome, 2002).
// Searches for a perfect numbering.
// In: v1 - pointer SEXP (double, but considered ind), with the first vertex of
//          the edge
//     v2 - pointer SEXP (double, but considered ind), with the second vertex of
//          the edge
//     p  - int, total number of vertices
//     varType - pointer SEXP, 0 if the vertex is continuous or 1 if discrete
//     from - starting vertex
// Out: pointer to a int vector with the perfect numbering: the first element
//      [0] indicates how many items in the vector, if it is zero means that the
//      graph is not triangulated; the second element [1] is the vertex numbered
//      as 1; [2] the vertex numbered as 2, ...
/******************************************************************************/
unsigned int *mcsFE(SEXP v1, SEXP v2, unsigned int p, SEXP numCat, unsigned int from, unsigned int **neighbourhood)
{
  bool bTriangulated, found;
  unsigned int i, j, n_edges, v, first;
  unsigned int *ind2, *ind3, *ind4, *num_neigh, *numbered, *perfN;
  unsigned int n_numb = 0;
//  unsigned int **neighbourhood;


  if (ISNA(INTEGER(v1)[0]) | ISNA(INTEGER(v2)[0]))
  {
    perfN = (unsigned int *)calloc(p+1,sizeof(unsigned int));
    perfN[0] = p;
    for (i=1;i<=p;i++) perfN[i] = i;
    return(perfN);
  }

  bTriangulated = true;
  n_edges = length(v1);

  if (from == 0)
  {
    // start from a discrete vertex
    i = v = 1;
    found = false;
//    while ((i<=p) & !found)
//      if (INTEGER(numCat)[i-1] != 0) found = true;
//      else i++;
    while ((i<=p) & !found)
      (INTEGER(numCat)[i-1] != 0) ? found = true : i++;
    if (found) v = i;
  }
  else
    v = from;

  numbered = (unsigned int *)calloc(p+1,sizeof(unsigned int)); // initialise with zeros
  perfN = (unsigned int *)calloc(p+1,sizeof(unsigned int)); // initialise with zeros
  numbered[0] = p; //number of elements in the vector
  perfN[0] = p; //idem
  numbered[v] = 1; //v is numbered as 1
  perfN[1] = v; //and is the first element in the result
  n_numb++; //number of already numbered vertices
  num_neigh = (unsigned int *)malloc((p+1)*sizeof(unsigned int));
  num_neigh[0] = p;  //number of elements in the vector

  //while not numbered all vertices and the graph is triangulated
  while ((n_numb<p) & (bTriangulated))
  {
    first = 0;
    for (i=1; i<=p; i++)
    {
      num_neigh[i] = 0;
      if (numbered[i]==0)
      {
        if (first==0) first = i;
        ind2 = neighbourhood[i];
        for (j=1; j<=ind2[0]; j++)
          if (numbered[ind2[j]] > 0) num_neigh[i]++;
      }
    }
    if (max(num_neigh)==0) //if there is no numbered neighbour
    {
      n_numb++;
      // give preference to discrete vertices
      i = 1;
      found = false;
//      while ((i<=p) & !found)
//        if ((INTEGER(numCat)[i-1] != 0) && (numbered[i]==0)) found = true;
//        else i++;
      while ((i<=p) & !found)
        ((INTEGER(numCat)[i-1] != 0) && (numbered[i]==0)) ? found = true : i++;
      if (!found) i = first;
      numbered[i] = n_numb;
      perfN[n_numb] = i;
    }
    else
    {
      ind2 = whichMaxFE(num_neigh); //the one(s) with more numbered neighbours
      // give preference to discrete vertices
      i = 1;
      found = false;
//      while ((i<=ind2[0]) & !found)
//        if (INTEGER(numCat)[ind2[i]-1] != 0) found = true;
//        else i++;
      while ((i<=ind2[0]) & !found)
        (INTEGER(numCat)[ind2[i]-1] != 0) ? found = true : i++;
      if (!found) i = 1;
      v = ind2[i];
      free(ind2);
      n_numb++;
      numbered[v] = n_numb;
      perfN[n_numb] = v;
      ind2 = neighbourhood[v];
      ind3 = which(numbered,0,false);
      ind4 = intersect(ind2,ind3);
      ind4[0]++;
      ind4[ind4[0]] = v;
      free(ind3);

      if (!isComplete(v1,v2,ind4,p)) //tests if it complete
        bTriangulated = false;
      free(ind4);
    }
  }
  if (!bTriangulated) perfN[0] = 0;

  free(num_neigh);
  free(numbered);
  return(perfN);
}

/******************************************************************************/
// Like R's "c". Binds 2 vectors of int - c(A,B)
// In: A - pointer to a int vector
//     B - pointer to a int vector
// Out: pointer to a int vector with A[0]+B[0]+1 elements.
/******************************************************************************/
unsigned int *bindFE(unsigned int *A, unsigned int *B)
{
  unsigned int i, *result;
  
  result = (unsigned int *)malloc((A[0]+B[0]+1)*sizeof(unsigned int));
  for (i=1;i<=A[0];i++)
    result[i] = A[i];
  for (i=1;i<=B[0];i++)
    result[A[0]+i] = B[i];
  result[0] = A[0] + B[0];
  return(result);
}

/******************************************************************************/
// Like R's "unique". Given a vector of int, returns a vector with the same
// elements, but with no repeats.
// In: A - pointer to a int vector
// Out: pointer to a int vector.
/******************************************************************************/
unsigned int *uniqueFE(unsigned int *A)
{
  unsigned int *B, *R, *aux;
  unsigned int i;

  B = (unsigned int *)malloc((A[0]+1)*sizeof(unsigned int));
  B[0] = 0;
  for (i=1; i<=A[0]; i++)
  {
    aux = which(B,A[i],true);
    if (aux[0] == 0)
    {
      B[0]++;
      B[B[0]] = A[i];
    }
    free(aux);
  }
  R = (unsigned int *)malloc((B[0]+1)*sizeof(unsigned int));
  for (i=0;i<=B[0];i++) R[i] = B[i];
  free(B);
  return(R);
}

/******************************************************************************/
// Given a vector of pointers to vectors of int, the number of vectors, and a
// particular vector, search all vectors to find out if the particular one is
// among the possible ones.
// In: S - vector of pointers (each [i] points to a vector of int)
//     numC - number of pointers in S
//     set - pointer to a vector of int
// Out: boolean (true if the set is in S; false otherwise).
/******************************************************************************/
bool findSetFE(unsigned int *S[], unsigned int numC, unsigned int *set)
{
  bool found = false;
  unsigned int k = 1;

  while ((!found) & (k<=numC))
  {
    found = setEqual(S[k-1],set);
    k++;
  }
  return(found);
}

/******************************************************************************/
// Sometimes it is easier to duplicate a vector, if 2 pointers make reference
// to the same vector. In this way, there is no chance of releasing a vector
// that is still being used, or trying to release more than once. This function
// just makes a copy of a vector of int.
// In: A - vector of int
// Out: pointer to a vector of int.
/******************************************************************************/
unsigned int *duplFE(unsigned int *A)
{
  unsigned int *result = (unsigned int *)malloc((A[0]+1)*sizeof(unsigned int));
  unsigned int i;
  
  for (i=0; i<=A[0]; i++) result[i] = A[i];
  return(result);
}

/******************************************************************************/
// Test if A is a subset of B (assuming that they are never equal).
// In: A - pointer to int vector
//     B - pointer to int vector
// Out: boolean.
/******************************************************************************/
bool subSetFE(unsigned int *A, unsigned int *B)
{
  bool isSS = true;
  unsigned int k = 1;
  unsigned int *ind;

  if ((A==NULL) | (B==NULL)) //if at least one is NULL
    isSS = false;
  else
    //no repetitions, so if A is bigger, cannot be a subset
    //assume that they are never equal
    if ((A[0] > B[0]) | (A[0] == B[0]))
      isSS = false;
    else
    {
      while (isSS & (k<=A[0])) //for all elements in one edge, seaches the other
      {
        ind = which(B,A[k],true);
        if (ind[0] != 1) isSS = false; // no repetitions allowed
        free(ind);
        k++;
      }
    }

  return(isSS);
}

/******************************************************************************/
/* verifies if x is in A                                                      */
/******************************************************************************/
bool isEl(unsigned int x, unsigned int *A)
{
  bool found=false;
  unsigned int i=1;
  if (A[0]!=0)
    while ((i<=A[0]) && (!found)) found = (x==A[i++]);
  return(found);
}

/******************************************************************************/
/* tests if S is a true separator of C1 and C2 in the graph giben by (v1,v2)  */
/******************************************************************************/
bool isTrueSep(SEXP v1, SEXP v2, unsigned int *C1, unsigned int *C2,
               unsigned int *S, unsigned int p, unsigned int **neighbourhood)
{
  unsigned int *L; //list of vertices that have to be visited; always the last
                   //element is the next to be visited
  unsigned int w, v;
  unsigned int *arr, *aux, *C;
  char *visited;
  bool found = false;

  C1 = setDiff(C1,S);
  C2 = setDiff(C2,S);
  if (C1[0] > C2[0])
  {
    C = C1;
    v = C2[1];
  }
  else
  {
    C = C2;
    v = C1[1];
  }

  visited = (char*)calloc(p+1,sizeof(char));
  visited[v] = 1; //the reference vertex starts as visited

  L = setDiff(neighbourhood[v],S);
  while ((L[0] > 0) && !found)
  {
    w = L[L[0]--]; //removes the next vertex to be visited from the list
    if (visited[w]==0) //hasn't already been visited
    {
      found = isEl(w,C);
      if (!found)
      {
        visited[w] = 1;
        aux = setDiff(neighbourhood[w],S);//the neighbours of w have to be visited
        arr = bindFE(aux,L);
        free(L);
        free(aux);
        L = arr;
      }
    }
  }
  free(visited);
  free(L);
  return(!found);
}


/******************************************************************************/
// The main function. Given a list of edges, seaches for possible edges that
// could be added. Before adding one edge is needed to test if it doesn't
// create a cycle. Returns the list of possible edges, the index in a previous
// list of edges (if it apperead in a previous iteration), the index to the
// minimal separator in the graph, and the list of minimal separators. It also
// returns the number of edges that must be considered in the list.
// In: v1 - pointer SEXP (double, but considered int), with the first vertex of
//          the edge
//     v2 - pointer SEXP (double, but considered int), with the second vertex of
//          the edge
//     pp - pointer SEXP to the number of vertices
//     previous - pointer SEXP, edges that accured in a previous iteration (the
//                indexes returned are the indexes to the elements in this
//                vector)
//     varType - pointer SEXP, 0 if the vertex is continuous or 1 if discrete
//     FROM
//     EXACT - "boolean"
// Out: R list
//        edges - matrix k by 5:
//                  column 1: first vertex in the edge
//                  column 2: second vertex in the edge
//                  column 3: index to the separator
//                  column 4: index in the previous list
//                  column 5: length of the separator (internal use only)
//        S - list with kk elements (separators)
//        total - number of lines that must be used from the matrix edges (<=k)
/******************************************************************************/
SEXP findEd(SEXP v1, SEXP v2, SEXP pp, SEXP previous, SEXP numCat, SEXP FROM,
            SEXP EXACT, SEXP JOIN)
{
  unsigned int i, j, k, x, y, z, ii, jj, iii, kk, t, numB = 0;
  unsigned int numNewEdges, numPr = 0;
  unsigned int *perfN, *perfNaux, *aux, *aux1;
  unsigned int *newEdges1, *newEdges2, *newEdges3;
  unsigned int numC = 0; //number of sets
  bool disjoint, include, join, trueSep, found;
  char exact;
  unsigned int from, p, *vert;

  exact = INTEGER(EXACT)[0];
  join = INTEGER(JOIN)[0];
  from = INTEGER(FROM)[0];
  p = INTEGER(pp)[0];
  vert = (unsigned int*)malloc(p*sizeof(unsigned int));
  unsigned int *C[p-1]; //perfect sequence
  struct nodeFE *newEdges = NULL;
  struct nodeFE *currNewEdges = NULL;
  unsigned int *H1 = NULL;
  unsigned int *R1 = NULL;
  unsigned int *Hdisj = NULL;
  unsigned int **neighbourhood;
  SEXP S1;

  neighbourhood = (unsigned int **)malloc((p+1)*sizeof(unsigned int*));
  for (i=1;i<=p;i++) neighbourhood[i] = findNeigh(v1,v2,i,p);
/******************************************************************************/
// find a perfect sequence
  perfNaux = (unsigned int *)calloc(p+1,sizeof(unsigned int));
  perfNaux[0] = 0;
  perfN = mcsFE(v1,v2,p,numCat,from,neighbourhood); //the perfect numbering
  for (i=1;i<=p;i++)
  {
    perfNaux[0]++;
    perfNaux[i] = perfN[i];
    aux = intersect(neighbourhood[perfN[i]],perfNaux);
    aux[++aux[0]] = perfN[i];
    C[numC++] = aux;
  }
  free(perfN);

/******************************************************************************/
// to remove what is not a clique (subset of a clique)
  t = 0;
  numB = numC;
  while (t < numC)
  {
    if (C[t]!=NULL)
    {
      i = 0;
      found = false;
      while ((i<numC) & (!found))
      {
        if (t != i) found = subSetFE(C[t],C[i]);
        i++;
      }
      i--;
      if (found)
      {
        numB--;
        free(C[t]);
        if (i<t)
          C[t] = NULL;
        else
        {
          C[t] = C[i];
          C[i] = NULL;
          t--; //as C[i] is in t-th position, this one has to be tested again
        }
      }
    }
    t++;
  }
  for (i=0; i<numC; i++)
  {
    if (C[i]==NULL)
    {
      j = i;
      while ((C[j]==NULL) & (j<numC-1)) j++;
      C[i] = C[j];
      C[j] = NULL;
    }
  }
  numC = numB;
  unsigned int *H[numC]; //list of histories
  unsigned int *R[numC]; //list of residuals
  unsigned int *S[numC]; //list of separators
/******************************************************************************/
// find histories (H), residuals (R), ans separators (S)
  S[0] = NULL;
  H[0] = duplFE(C[0]);
  R[0] = duplFE(C[0]);
  for (i=2; i<=numC; i++)
  {
    aux  = bindFE(H[i-2],C[i-1]); //remember that the index starts at 0
    aux1 = uniqueFE(aux);
    free(aux);
    H[i-1] = aux1;
    aux  = setDiff(C[i-1],H[i-2]);
    aux1 = uniqueFE(aux);
    free(aux);
    R[i-1] = aux1;
    S[i-1] = intersect(H[i-2],C[i-1]);
  }
/******************************************************************************/
// find add-eligible edges
  numNewEdges = 0;

  for (j=2; j<=numC; j++)
  {
    i = j;
    disjoint = false;
    while ((!disjoint) & (i>1))
    {
      if (S[i-1] == NULL) disjoint = true;
      else
        if (S[i-1][0] == 0) disjoint = true;
      i--;
      include = false;
      if (disjoint && join)
      {
        include = true;
        R1 = duplFE(R[j-1]);
        H1 = duplFE(H[i-1]);
      }
      else
      {
        Hdisj = intersect(C[i-1],C[j-1]);
        if (Hdisj[0] != 0) include = findSetFE(S,j,Hdisj);
        if (include)
        {
          trueSep = true;
          if (exact > 0) trueSep = isTrueSep(v1,v2,C[i-1],C[j-1],Hdisj,p,neighbourhood);
          include = trueSep;
          if (include)
          {
            R1 = duplFE(R[j-1]);
            H1 = setDiff(C[i-1],C[j-1]);
          }
        }
      }
      if (include)
      {
        numNewEdges += H1[0]*R1[0];
        newEdges1 = (unsigned int *)calloc((H1[0]*R1[0]+1),sizeof(unsigned int)); //first vertex
        newEdges2 = (unsigned int *)calloc((H1[0]*R1[0]+1),sizeof(unsigned int)); //second vertex
        newEdges3 = (unsigned int *)calloc((H1[0]*R1[0]+1),sizeof(unsigned int)); //index in previous
        for (jj=1; jj<=H1[0]; jj++)
          for (ii=1; ii<=R1[0]; ii++)
          {
            x = imin2(R1[ii],H1[jj]);
            y = imax2(R1[ii],H1[jj]);
            z = (x-1)*p-(x-1)*x/2 + y-x;
            newEdges1[0]++;
            newEdges2[0]++;
            newEdges3[0]++;
            newEdges1[newEdges1[0]] = x;
            newEdges2[newEdges2[0]] = y;
            aux1 = whichS(previous,z,true);
            if ((aux1[0] > 0) & (x!=INTEGER(v1)[length(v1)-1]) &
                (y!=INTEGER(v2)[length(v2)-1]) & (y!=INTEGER(v1)[length(v1)-1]) &
                (x!=INTEGER(v2)[length(v2)-1]))
              newEdges3[newEdges3[0]] = aux1[aux1[0]]+1;
              //the index in R starts at 1; it was initialized with 0s
            free(aux1);
          }
        if (newEdges1[0] > 0)
        {
          if (!disjoint)
          {
            PROTECT(S1 = allocVector(INTSXP, Hdisj[0]));
            for (iii=1; iii<=Hdisj[0]; iii++) INTEGER(S1)[iii-1] = Hdisj[iii];
            free(Hdisj);
          }
          else
            PROTECT(S1 = allocVector(INTSXP,0));
          numPr++;
          if (newEdges == NULL)
          {
            newEdges = newNodeFE(newEdges1,newEdges2,newEdges3,S1);
            currNewEdges = newEdges;
          }
          else
          {
            currNewEdges->next = newNodeFE(newEdges1,newEdges2,newEdges3,S1);
            currNewEdges = currNewEdges->next;
          }
        }
        else
        {
          free(newEdges1);
          free(newEdges2);
          free(newEdges3);
        }
        free(H1);
        free(R1);
      }
      else
        if (!disjoint) free(Hdisj);
    }
  }
  free(vert);

  for (i=1; i<=numC; i++)
  {
    free(C[i-1]);
    free(H[i-1]);
    free(R[i-1]);
    free(S[i-1]);
  }

  // note that the first position in neighbourhood was never used
  for (i=1;i<=p;i++) free(neighbourhood[i]);
  free(neighbourhood);

/******************************************************************************/
// assemble the results in R structure
  SEXP subList; //the list with the separators
  PROTECT(subList = allocVector(VECSXP,numPr));
  numPr++;

  SEXP edges; //matrix with all edges information
  PROTECT(edges = allocMatrix(INTSXP, numNewEdges, 5));
  numPr++;

  SEXP listNames; //the names (in R) of the result components
  char *names[3] = {"edges", "S", "total"};
  PROTECT(listNames = allocVector(STRSXP,3));
  numPr++;
  SET_STRING_ELT(listNames,0,mkChar(names[0]));
  SET_STRING_ELT(listNames,1,mkChar(names[1]));
  SET_STRING_ELT(listNames,2,mkChar(names[2]));

  SEXP list; //the final result
  PROTECT(list = allocVector(VECSXP,3));
  numPr++;

  SEXP total; //number of lines in "edges" that must be considered
  PROTECT(total = allocVector(INTSXP,1));
  numPr++;
  SET_VECTOR_ELT(list, 0, edges);
  SET_VECTOR_ELT(list, 1, subList);
  SET_VECTOR_ELT(list, 2, total);
  setAttrib(list, R_NamesSymbol, listNames);
  
  k = 0; //index in edges
  kk = 0; //index in S (sublist)
  unsigned int *edgesAdded = (unsigned int *)malloc((numNewEdges+1)*sizeof(unsigned int));
  edgesAdded[0] = 0;
  unsigned int *indexEdgesAdded = (unsigned int *)malloc((numNewEdges+1)*sizeof(unsigned int));
  indexEdgesAdded[0] = 0;
  while (newEdges != NULL)
  {
    for (i=1; i<=newEdges->vert1[0]; i++)
    {
      x = newEdges->vert1[i];
      y = newEdges->vert2[i];
      z = (x-1)*p-(x-1)*x/2 + y-x;
      aux = which(edgesAdded,z,true); //see if the edge has already been added
      if (aux[0] == 0)
      {
        INTEGER(edges)[k] = newEdges->vert1[i]; //first vertex
        INTEGER(edges)[1*numNewEdges+k] = newEdges->vert2[i]; //second vertex
        INTEGER(edges)[2*numNewEdges+k] = kk + 1; //index for the separator in S
        INTEGER(edges)[3*numNewEdges+k] = newEdges->index[i]; //index previous
        //separator's length: only cliques should be considered in the perfect
        //sequence, but it makes the algorithm slower. So, here some subsets of
        //the cliques are considered. The problem is that the algorithm finds
        //some repeated edges, and with a incomplete separator. To fix this,
        //the separator must be the largest set among the separators found for
        //such edge.
        INTEGER(edges)[4*numNewEdges+k] = length(newEdges->S);
        //these below are just to make easier to find an added edge among the
        //added ones.
        edgesAdded[0]++;
        edgesAdded[edgesAdded[0]] = z;
        indexEdgesAdded[0]++;
        indexEdgesAdded[indexEdgesAdded[0]] = k;
        k++;
      }
      else
      {
        if (length(newEdges->S) >
            INTEGER(edges)[4*numNewEdges+indexEdgesAdded[aux[1]]])
        {
          INTEGER(edges)[2*numNewEdges+indexEdgesAdded[aux[1]]] = kk + 1;
          INTEGER(edges)[4*numNewEdges+indexEdgesAdded[aux[1]]] = length(newEdges->S);
        }
      }
      free(aux);
    }
    SET_VECTOR_ELT(subList, kk, newEdges->S);
    kk++;
    currNewEdges = newEdges;
    newEdges = newEdges->next;
    free(currNewEdges->vert1);
    free(currNewEdges->vert2);
    free(currNewEdges->index);
    free(currNewEdges);
  }
  INTEGER(total)[0] = k;
  free(edgesAdded);
  free(indexEdgesAdded);

  UNPROTECT(numPr);
  return(list);
}


/******************************************************************************/
// Given a list of edges, gets the cliques, separators, histories and residuals.
// In: v1 - pointer SEXP (double, but considered int), with the first vertex of
//          the edge
//     v2 - pointer SEXP (double, but considered int), with the second vertex of
//          the edge
//     pp - pointer SEXP to the number of vertices
//     varType - pointer SEXP, 0 if the vertex is continuous or 1 if discrete
//     from - inicial vertex for MCS
// Out: R list
//        cliques - list
//        histories - list
//        separators - list
//        residuals - list
/******************************************************************************/
SEXP perfSets(SEXP v1, SEXP v2, SEXP pp, SEXP numCat, SEXP FROM)
{
  unsigned int i, j, p;
  unsigned int numPr=0;
  unsigned int *perfN, *perfNaux, *aux, *aux1;
  unsigned int numC = 0; //number of sets
  unsigned int from, **neighbourhood;
  SEXP list; //the final result
  from = INTEGER(FROM)[0];
  p = INTEGER(pp)[0];
  unsigned int *C[p-1]; //perfect sequence
  unsigned int *H[p-1]; //list of histories
  unsigned int *R[p-1]; //list of residuals
  unsigned int *S[p-1]; //list of separators

  neighbourhood = (unsigned int **)malloc((p+1)*sizeof(unsigned int*));
  for (i=1;i<=p;i++)
    neighbourhood[i] = findNeigh(v1,v2,i,p);
  //////////////////////////////////////////////////////////////////////////////
  // FIND THE PERFECT SEQUENCE
  perfNaux = (unsigned int *)calloc(p+1,sizeof(unsigned int));
  perfNaux[0] = 0;
  perfN = mcsFE(v1,v2,p,numCat,from,neighbourhood); //the perfect numbering
  if (perfN[0]!=0)
  {
    for (i=1;i<=p;i++)
    {
      perfNaux[0]++;
      perfNaux[i] = perfN[i];
      aux = intersect(neighbourhood[perfN[i]],perfNaux);
      aux[0]++;
      aux[aux[0]] = perfN[i];
      C[numC] = aux;
      numC++;
    }
    for (i=1;i<=p;i++) // note that the first position was never used
      free(neighbourhood[i]);
    free(neighbourhood);
    free(perfNaux);
    free(perfN);
  /******************************************************************************/
  // to remove what is not a clique (subset of a clique)
    unsigned int numB = 0;
    unsigned int t;
    bool found;

    t = 0;
    numB = numC;
    while (t < numC)
    {
      if (C[t]!=NULL)
      {
        i = 0;
        found = false;
        while ((i<numC) & (!found))
        {
          if (t != i) found = subSetFE(C[t],C[i]);
          i++;
        }
        i--;
        if (found)
        {
          numB--;
          free(C[t]);
          if (i<t)
            C[t] = NULL;
          else
          {
            C[t] = C[i];
            C[i] = NULL;
            t--; //as C[i] is in t-th position, this one has to be tested again
          }
        }
      }
      t++;
    }
    for (i=0; i<numC; i++)
    {
      if (C[i]==NULL)
      {
        j = i;
        while ((C[j]==NULL) & (j<numC-1)) j++;
        C[i] = C[j];
        C[j] = NULL;
      }
    }
    numC = numB;
/******************************************************************************/
    //////////////////////////////////////////////////////////////////////////////
    // FIND THE HISTORIES (H), RESIDUALS (R), AND SEPARATORS (S)
    S[0] = NULL;
    H[0] = duplFE(C[0]);
    R[0] = duplFE(C[0]);
    for (i=2; i<=numC; i++)
    {
      aux  = bindFE(H[i-2],C[i-1]); //remember that the index starts at 0
      aux1 = uniqueFE(aux);
      free(aux);
      H[i-1] = aux1;
      aux  = setDiff(C[i-1],H[i-2]);
      aux1 = uniqueFE(aux);
      free(aux);
      R[i-1] = aux1;
      S[i-1] = intersect(H[i-2],C[i-1]);
    }
    //////////////////////////////////////////////////////////////////////////////
    // ASSEMBLE THE RESULTS IN R STRUCTURE
    numPr = 0;
    SEXP cliques;
    PROTECT(cliques = allocVector(VECSXP,numC));
    numPr++;
    SEXP histories;
    PROTECT(histories = allocVector(VECSXP,numC));
    numPr++;
    SEXP separators;
    PROTECT(separators = allocVector(VECSXP,numC));
    numPr++;
    SEXP residuals;
    PROTECT(residuals = allocVector(VECSXP,numC));
    numPr++;
    SEXP X;
    for (i=0; i<numC; i++)
    {
      PROTECT(X = allocVector(INTSXP, C[i][0]));
      numPr++;
      for (j=1; j<=C[i][0]; j++)
        INTEGER(X)[j-1] = C[i][j];
      free(C[i]);
      SET_VECTOR_ELT(cliques, i, X);

      PROTECT(X = allocVector(INTSXP, H[i][0]));
      numPr++;
      for (j=1; j<=H[i][0]; j++)
        INTEGER(X)[j-1] = H[i][j];
      free(H[i]);
      SET_VECTOR_ELT(histories, i, X);

      if (i==0)
      {
        PROTECT(X = allocVector(REALSXP,0));
        numPr++;
      }
      else
      {
        PROTECT(X = allocVector(INTSXP, S[i][0]));
        numPr++;
        for (j=1; j<=S[i][0]; j++)
          INTEGER(X)[j-1] = S[i][j];
        free(S[i]);
        SET_VECTOR_ELT(separators, i, X);
      }

      PROTECT(X = allocVector(INTSXP, R[i][0]));
      numPr++;
      for (j=1; j<=R[i][0]; j++)
        INTEGER(X)[j-1] = R[i][j];
      free(R[i]);
      SET_VECTOR_ELT(residuals, i, X);
    }
    SEXP listNames; //the names (in R) of the result components
    char *names[4] = {"cliques", "histories", "separators", "residuals"};
    PROTECT(listNames = allocVector(STRSXP,4));
    numPr++;
    SET_STRING_ELT(listNames,0,mkChar(names[0]));
    SET_STRING_ELT(listNames,1,mkChar(names[1]));
    SET_STRING_ELT(listNames,2,mkChar(names[2]));
    SET_STRING_ELT(listNames,3,mkChar(names[3]));
    PROTECT(list = allocVector(VECSXP,4));
    numPr++;
    SET_VECTOR_ELT(list, 0, cliques);
    SET_VECTOR_ELT(list, 1, histories);
    SET_VECTOR_ELT(list, 2, separators);
    SET_VECTOR_ELT(list, 3, residuals);
    setAttrib(list, R_NamesSymbol, listNames);
  }
  else
  {
    for (i=1;i<=p;i++) // note that the first position was never used
      free(neighbourhood[i]);
    free(neighbourhood);
    free(perfNaux);
    free(perfN);
    PROTECT(list = allocVector(INTSXP,1));
    numPr++;
    INTEGER(list)[0] = 0;
  }

  UNPROTECT(numPr);
  return(list);
}

/******************************************************************************/
// In: C - array with pointers to the cliques
//     S - array with pointers to the separators
//     numC - number of cliques
//     Sep - array with integers indicating which separators are to be kept
//           after removing copies
//     jT1,jT2 - are the lists of edges in the junction tree (cliques)
//     subSep - array with pointers indicating of which others separators in
//              Sep (i.e. after deleting copies) each separator is a subset
//     indSep - array with integers indicating that the original separator is
//              in position x in Sep
//     numS - pointer to the number os separators in Sep
// Out:
/******************************************************************************/
void jTree(unsigned int *C[],unsigned int *S[],unsigned int numC,
           unsigned int Sep[], unsigned int *jT1, unsigned int *jT2,
           unsigned int *subSep[], unsigned int *indSep, unsigned int *numS)
{
  unsigned int i, j;
  unsigned int *aux;
  bool STOP;

  char *del = (char*)calloc(numC,sizeof(char));
  if (numC==1) //if numC==1 the graph is complete
  {
    jT1[0] = 0;
    jT2[0] = 0;
    numS[0] = 0;
    indSep[0] = 0;
  }
  else
  {
    Sep[0] = 0;
    numS[0] = 1;
    indSep[0] = 0;

    unsigned int *cumsum = (unsigned int *)calloc(numC,sizeof(unsigned int));

    for (i=0;i<numC-1;i++)
      if (del[i]==0)
      {
        if (i>0)
          cumsum[i] = cumsum[i-1];
        for (j=i+1;j<numC;j++)
          if (del[j]==0)
          {
            if (setEqual(S[i],S[j]) | (S[j][0]==0))
            {
              del[j] = 1;
              indSep[j] = i;
            }
            else
              indSep[j] = j;
          }
      }
      else
        cumsum[i] = cumsum[i-1]+1;
    cumsum[numC-1] = cumsum[numC-2] + del[numC-1];
    for (i=1;i<numC;i++)
    {
      jT2[i] = i+1;
      if (del[i]==0) Sep[numS[0]++] = i;
      indSep[i] = indSep[i] - cumsum[indSep[i]];
    }
    free(cumsum);
  }
  free(del);

  if (numS[0]>1)
  {
    subSep[0] = (unsigned int *)malloc((numS[0]-1+1)*sizeof(unsigned int));
    subSep[0][0] = numS[0]-1;
    for (i=1;i<=numS[0]-1;i++)
    {
      subSep[0][i] = i;
      subSep[i] = (unsigned int *)malloc((numS[0]-2+1)*sizeof(unsigned int));
      subSep[i][0] = 0;
    }
    if (numS[0]>2)
      for (i=1;i<numS[0]-1;i++)
      {
        if (S[Sep[i]][0]==0) subSep[i][++subSep[i][0]] = 0;
        for (j=i+1;j<numS[0];j++)
        {
          aux = intersect(S[Sep[i]],S[Sep[j]]);
          if (aux[0]==S[Sep[i]][0]) subSep[i][++subSep[i][0]] = j;//no problems with empty
          else //never equal
            if (aux[0]==S[Sep[j]][0]) subSep[j][++subSep[j][0]] = i;
          free(aux);
        }
      }
  }

  //find a junction tree
  for (i=1;i<numC;i++)
  {
    STOP = false;
    j = i-1;
    if (S[i][0]!=0)
      while((!STOP) && (j>=0))
      {
        aux = intersect(S[i],C[j]);
        STOP = aux[0]==S[i][0];
        if (STOP)
        {
          jT1[0]++;
          jT2[0]++;
          jT1[jT1[0]] = j;
          jT2[jT2[0]] = i;
        }
        j--;
        free(aux);
      }
  }
}

/******************************************************************************/
// In: v1,v2 - list os edges
//     pp - number of vertices
//     varType - type os each vertex
// Out: list wth -
//      separators - list of unique separators
//      juncTree - matrix (k by 2) with the edges in the junction tree
//      sepSubSetOfSep - sepSubSetOfSep[i] holds the list of separators that
//                       include separators[i]
//      indSepOrig - separator[i] in the result of MCS is in position
//                   indSepOrig[i] in separators
//      cliques - cliques from MCS
/******************************************************************************/
SEXP juncTree(SEXP v1, SEXP v2, SEXP pp, SEXP numCat)
{
  unsigned int i, j;
  unsigned int numPr = 0;
  unsigned int *perfN, *perfNaux, *aux, *aux1;
  unsigned int numC = 0; //number of sets
  unsigned int from = 0, **neighbourhood;
  short unsigned p = INTEGER(pp)[0];
  unsigned int *C[p-1]; //perfect sequence
  unsigned int *H[p-1]; //list of histories
  unsigned int *R[p-1]; //list of residuals
  unsigned int *S[p-1]; //list of separators

  neighbourhood = (unsigned int **)malloc((p+1)*sizeof(unsigned int*));
  for (i=1;i<=p;i++)
    neighbourhood[i] = findNeigh(v1,v2,i,p);
  //////////////////////////////////////////////////////////////////////////////
  // FIND THE PERFECT SEQUENCE
  perfNaux = (unsigned int *)calloc(p+1,sizeof(unsigned int));
  perfNaux[0] = 0;
  perfN = mcsFE(v1,v2,p,numCat,from,neighbourhood); //the perfect numbering
  for (i=1;i<=p;i++)
  {
    perfNaux[0]++;
    perfNaux[i] = perfN[i];
    aux = intersect(neighbourhood[perfN[i]],perfNaux);
    aux[0]++;
    aux[aux[0]] = perfN[i];
    C[numC] = aux;
    numC++;
  }
  free(perfNaux);
  free(perfN);
/******************************************************************************/
// to remove what is not a clique (subset of a clique)
  unsigned int numB = 0;
  unsigned int t;
  bool found;
  t = 0;
  numB = numC;
  while (t < numC)
  {
    if (C[t]!=NULL)
    {
      i = 0;
      found = false;
      while ((i<numC) & (!found))
      {
        if (t != i) found = subSetFE(C[t],C[i]);
        i++;
      }
      i--;
      if (found)
      {
        numB--;
        free(C[t]);
        if (i<t)
          C[t] = NULL;
        else
        {
          C[t] = C[i];
          C[i] = NULL;
          t--; //as C[i] is in t-th position, this one has to be tested again
        }
      }
    }
    t++;
  }
  for (i=0; i<numC; i++)
  {
    if (C[i]==NULL)
    {
      j = i;
      while ((C[j]==NULL) & (j<numC-1))
        j++;
      C[i] = C[j];
      C[j] = NULL;
    }
  }
  numC = numB;
/******************************************************************************/
  for (i=1;i<=p;i++) // note that the first position was never used
    free(neighbourhood[i]);
  free(neighbourhood);
  //////////////////////////////////////////////////////////////////////////////
  // FIND THE HISTORIES (H), RESIDUALS (R), AND SEPARATORS (S)
  S[0] = NULL;
  H[0] = duplFE(C[0]);
  R[0] = duplFE(C[0]);
  for (i=2; i<=numC; i++)
  {
    aux  = bindFE(H[i-2],C[i-1]); //remember that the index starts at 0
    aux1 = uniqueFE(aux);
    free(aux);
    H[i-1] = aux1;
    aux  = setDiff(C[i-1],H[i-2]);
    aux1 = uniqueFE(aux);
    free(aux);
    R[i-1] = aux1;
    S[i-1] = intersect(H[i-2],C[i-1]);
  }

  unsigned int Sep[numC];
  unsigned int *jT1 = (unsigned int *)calloc(numC+1,sizeof(unsigned int));
  unsigned int *jT2 = (unsigned int *)calloc(numC+1,sizeof(unsigned int));
  unsigned int *subSep[numC];
  unsigned int *indSep = (unsigned int*)malloc(numC*sizeof(unsigned int));
  unsigned int *numS = (unsigned int*)calloc(1,sizeof(unsigned int));

  jTree(C,S,numC,Sep,jT1,jT2,subSep,indSep,numS);

  SEXP SEP,JT,SUBSEP,INDSEP,list,listNames,AUX,CLIQUES;
  char *names[5] = {"separators","juncTree","sepSubSetOfSep","indSepOrig","cliques"};
  PROTECT(listNames = allocVector(STRSXP,5));
  numPr++;
  SET_STRING_ELT(listNames,0,mkChar(names[0]));
  SET_STRING_ELT(listNames,1,mkChar(names[1]));
  SET_STRING_ELT(listNames,2,mkChar(names[2]));
  SET_STRING_ELT(listNames,3,mkChar(names[3]));
  SET_STRING_ELT(listNames,4,mkChar(names[4]));

  PROTECT(SEP = allocVector(VECSXP,numS[0]));
  numPr++;
  PROTECT(AUX = allocVector(INTSXP,0));
  numPr++;
  SET_VECTOR_ELT(SEP, 0, AUX);
  for (i=1;i<numS[0];i++)
  {
    PROTECT(AUX = allocVector(INTSXP,S[Sep[i]][0]));
    numPr++;
    for (j=1;j<=S[Sep[i]][0];j++)
      INTEGER(AUX)[j-1] = S[Sep[i]][j];
    SET_VECTOR_ELT(SEP, i, AUX);
  }

  PROTECT(JT = allocMatrix(INTSXP,jT1[0],2));
  numPr++;
  for (i=0;i<jT1[0];i++)
  {
    INTEGER(JT)[i] = jT1[i+1] + 1; //+1 because the index here starts at 0
    INTEGER(JT)[jT1[0]+i] = jT2[i+1] + 1; //and in R starts at 1
  }
  free(jT1);
  free(jT2);

  PROTECT(SUBSEP = allocVector(VECSXP,numS[0]));
  numPr++;
  for (i=0;i<numS[0];i++)
  {
    PROTECT(AUX = allocVector(INTSXP,subSep[i][0]));
    numPr++;
    for (j=1;j<=subSep[i][0];j++)
      INTEGER(AUX)[j-1] = subSep[i][j]+1;//+1 as above
    SET_VECTOR_ELT(SUBSEP, i, AUX);
    free(subSep[i]);
  }

  PROTECT(INDSEP = allocVector(INTSXP,numC));
  numPr++;
  for (i=0;i<numC;i++)
    INTEGER(INDSEP)[i] = indSep[i]+1;//+1 as above
  free(indSep);
  free(numS);

  PROTECT(CLIQUES = allocVector(VECSXP,numC));
  numPr++;
  for (i=0; i<numC; i++)
  {
    PROTECT(AUX = allocVector(INTSXP,C[i][0]));
    numPr++;
    for (j=1;j<=C[i][0];j++)
      INTEGER(AUX)[j-1] = C[i][j];
    SET_VECTOR_ELT(CLIQUES, i, AUX);
    free(C[i]);
    free(R[i]);
    free(H[i]);
    if (i>0)
      free(S[i]);
  }

  PROTECT(list = allocVector(VECSXP,5));
  numPr++;
  SET_VECTOR_ELT(list, 0, SEP);
  SET_VECTOR_ELT(list, 1, JT);
  SET_VECTOR_ELT(list, 2, SUBSEP);
  SET_VECTOR_ELT(list, 3, INDSEP);
  SET_VECTOR_ELT(list, 4, CLIQUES);
  setAttrib(list, R_NamesSymbol, listNames);

  UNPROTECT(numPr);
  return(list);
}
