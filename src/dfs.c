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
// The main function.
// Depth-first search, searches for all vertices reacheable from one specifc
// vertex (considering that there is no cycle).
// In: v1 - pointer SEXP (double, but considered int), with the first vertex of
//          the edge
//     v2 - pointer SEXP (double, but considered int), with the second vertex of
//          the edge
//     v  - int, the reference vertex
//     p  - int, total number of vertices
// Out: vector with all vertices reacheable from v
/******************************************************************************/
SEXP dfs(SEXP v1, SEXP v2, SEXP vv, SEXP pp)
{
  unsigned int *L; //list of vertices that have to be visited; always the last
                   //element is the next to be visited
  unsigned int *K; //list of vertices reacheable from the reference vertex
  unsigned int i, p, w, v;
  unsigned int *arr, *aux;
  SEXP result;

  v = INTEGER(vv)[0]; //reference vertex
  p = INTEGER(pp)[0]; //total number of vertices
  char visited[p+1]; //list of visited vertices (1 if visited, 0 otherwise)
  for (i=0;i<=p;i++)
    visited[i] = 0;
  visited[v] = 1; //the reference vertex starts as visited
  L = findNeigh(v1,v2,v,p);
  K = (unsigned int *)malloc((p+1)*sizeof(unsigned int));
  K[0] = 0;
  while (L[0] > 0)
  {
    w = L[L[0]--]; //removes the next vertex to be visited from the list
    if (visited[w] == 0) //hasn't already been visited
    {
      K[++K[0]] = w; //w is reacheable
      visited[w] = 1;
      aux = findNeigh(v1,v2,w,p); //the neighbours of w have to be visited
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
  PROTECT(result = allocVector(INTSXP, K[0]));
  for (i=1; i<=K[0]; i++)
    INTEGER(result)[i-1] = K[i]; // because R indexing starts in 1
  free(K);
  free(L);
  UNPROTECT(1);
  return(result);
}
