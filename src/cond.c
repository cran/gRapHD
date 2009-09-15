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
#include "comfunc.h"


/******************************************************************************/
// Find the neighbours of a vertex. Similar to findNeigh in comfunc.c, but
// receives unsigned int instead of SEXP
// In: v1 - unsigned int, with the first vertex of the edge
//     v2 - unsigned int, with the second vertex of the edge
//     v  - unsigned int, the vertex which neighbours are required
//     p  - unsigned int, total number of vertices
// Out: pointer to unsigned int vector with the neighbours (vec[0]=number of
//      elements).
/******************************************************************************/
static unsigned int *findNeighC(unsigned int *v1, unsigned int *v2, 
                                unsigned int v, unsigned int p)
{
  unsigned int i, *ind1, *ind2, *result;

  ind1 = which(v1,v,true); //edges in which v is the lowest vertex
  ind2 = which(v2,v,true); //edges in which v is the highest vertex
  result = (unsigned int *)malloc((ind1[0]+ind2[0]+1)*sizeof(unsigned int));
  for (i=1;i<=ind1[0];i++)
    result[i] = v2[ind1[i]]; //gets the highest when v is the lowest
  for (i=1;i<=ind2[0];i++)
    result[ind1[0]+i] = v1[ind2[i]]; //opposite
  result[0] = ind1[0] + ind2[0];
  free(ind1);
  free(ind2);
  return result;
}
/******************************************************************************/

/******************************************************************************/
// Depth-first search, searches for all vertices reacheable from one specifc
// vertex (considering that there is no cycle). Similar to dfs.c, but the other
// returns R objects.
// In: v1 - unsigned int, with the first vertex of the edge
//     v2 - unsigned int, with the second vertex of the edge
//     v  - int, the reference vertex
//     p  - int, total number of vertices
// Out: vector with all vertices reacheable from v
/******************************************************************************/
static unsigned int *dfsS(unsigned int *v1, unsigned int *v2, unsigned int v, 
                          unsigned int p)
{
  unsigned int *L; //list of vertices that have to be visited; always the last
                   //element is the next to be visited
  unsigned int *K; //list of vertices reacheable from the reference vertex
  unsigned int i,w;
  unsigned int *arr, *aux;

  char visited[p+1]; //list of visited vertices (1 if visited, 0 otherwise)
  for (i=0;i<=p;i++)
    visited[i] = 0;
  visited[v] = 1; //the reference vertex starts as visited
  L = findNeighC(v1,v2,v,p);
  K = (unsigned int *)malloc((p+1)*sizeof(unsigned int));
  K[0] = 0;
  while (L[0] > 0)
  {
    w = L[L[0]--]; //removes the next vertex to be visited from the list
    if (visited[w] == 0) //hasn't already been visited
    {
      K[++K[0]] = w; //w is reacheable
      visited[w] = 1;
      aux = findNeighC(v1,v2,w,p); //the neighbours of w have to be visited
      arr = (unsigned int *)malloc((L[0]+aux[0]+1)*sizeof(unsigned int));
      arr[0] = 0;
      for (i=1;i<=aux[0];i++)
        arr[++arr[0]] = aux[i];
      for (i=1;i<=L[0];i++)
        arr[++arr[0]] = L[i];
      free(L);
      free(aux);
      L = arr;
    }
  }
  free(L);
  return(K);
}
/******************************************************************************/

/******************************************************************************/
// Finds out if 2 sets are disjoint.
// In: A - SEXP
//     B - SEXP
// Out: boolean (T if disjoint).
/******************************************************************************/
static bool isDisjoint(SEXP A, SEXP B)
{
  unsigned int nR, i, j;
  bool inc;
  bool disjoint = true;
  SEXP C, D;

  if (length(A) > length(B)) //selects the shortest vector to run through
  {
    nR = length(B);
    C = B;
    D = A;
  }
  else
  {
    nR = length(A);
    C = A;
    D = B;
  }

  i = 0;
  while ((i<nR) && disjoint)
  {
    j = 0;
    inc = false;
    while (!inc && (j<length(D)))
      inc = INTEGER(C)[i]==INTEGER(D)[j++];
    if (inc)
      disjoint = false;
    i++;
  }

  return(disjoint);
}

/******************************************************************************/
// The minimum value in a vector.
/******************************************************************************/
static unsigned int minimum(unsigned int *A)
{
  unsigned int i,m;
  m = A[1];
  for (i=1;i<=A[0];i++)
    if (A[i]<m)
      m = A[i];
  return(m);
}
/******************************************************************************/

/******************************************************************************/
// The maximum value in a vector.
/******************************************************************************/
static unsigned int maximum(unsigned int *A)
{
  unsigned int i,m;
  m = A[1];
  for (i=1;i<=A[0];i++)
    if (A[i]>m)
      m = A[i];
  return(m);
}
/******************************************************************************/

/******************************************************************************/
// Determines the components in a graph.
// In: cond - pointer SEXP to a list in which each element is a clique
// Out: int vector with the number of the componet that vertex is in.
/******************************************************************************/
SEXP components(SEXP cond)
{
  unsigned int i,j,n,M;
  unsigned int *v1,*v2,*comp,*aux;
  n = length(cond); //number of cliques
  SEXP result;
  // v1 and v2 are hyper-edges in a graph where teh cliques are the vertices
  // and 2 cliques are connected by an edge if their intersection is not empty
  v1 = (unsigned int*)malloc((n*(n-1)/2+1)*sizeof(unsigned int));
  v2 = (unsigned int*)malloc((n*(n-1)/2+1)*sizeof(unsigned int));
  v1[0] = v2[0] = 0;
  // finds the hyper-edges
  for (i=0; i<n-1; i++)
    for (j=i+1; j<n; j++)
      if (!isDisjoint(VECTOR_ELT(cond,i),VECTOR_ELT(cond,j)))
      {
        v1[++v1[0]] = i+1;        
        v2[++v2[0]] = j+1;        
      }
  
  comp = (unsigned int*)calloc((n+1),sizeof(unsigned int));
  comp[0] = n;
  i = 0;
  M = 0;
  while (minimum(comp)==0)
    if (comp[++i]==0)
    {
      aux = dfsS(v1,v2,i,n);
      M++;//x = maximum(comp);
      comp[i] = M;//x+1;
      for (j=1;j<=aux[0];j++)
        comp[aux[j]] = M;//x+1;
      free(aux);
    }
  free(v1);
  free(v2); 
  PROTECT(result=allocVector(INTSXP,n));
  for (i=0;i<n;i++)
    INTEGER(result)[i] = comp[i+1];
  UNPROTECT(1);
  return(result);
}
