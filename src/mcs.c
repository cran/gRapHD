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
// Like R's "which.max", returns the index for the maximum value in a vector.
// In: d - pointer to int vector
// Out: int, the index for the maximum value in d.
/******************************************************************************/
static unsigned int whichMaxMCS(unsigned int *d)
{
  unsigned int i, j, m;

  m = d[1];
  j = 1;
  for(i=2; i<=d[0]; i++)
    if (m < d[i])
    {
      m = d[i];
      j = i;
    }
  return(j);
}


struct nodeMCS
{
  unsigned int neighbour;
  struct nodeMCS *next;
};

/******************************************************************************/
// Maximum Cardinality Search (Tarjan & Yannakakis, in Jayson Rome, 2002).
// Searches for a perfect numbering.
// In: v1 - pointer SEXP (double, but considered int), with the first vertex of
//          the edge
//     v2 - pointer SEXP (double, but considered int), with the second vertex of
//          the edge
//     v  - pointer SEXP, initial vertex
//     p  - pointer SEXP, total number of vertices
// Out: int vector with the perfect numbering: if the first element is zero
//      means that the graph is not triangulated; otherwise, each position
//      gives the order of that vertex in the perfect numbering.
/******************************************************************************/
SEXP mcs(SEXP v1, SEXP v2, SEXP vv, SEXP pp)
{
  unsigned int i, j, p, v, n_edges, n_numb = 0, first;
  unsigned int *ind2, *ind3, *ind4, *num_neigh, *numbered;
  unsigned int **neighbourhood;
  bool bTriangulated = true;
  char numPr = 0;
  SEXP result;

  v = INTEGER(vv)[0];
  p = INTEGER(pp)[0];
  n_edges = length(v1);

  numbered = (unsigned int *)calloc(p+1,sizeof(unsigned int)); // initialise with zeros
  numbered[0] = p; //number of elements in the vector
  numbered[v] = 1; //v is numbered as 1
  n_numb++; //number of already numbered vertices
  num_neigh = (unsigned int *)calloc(p+1,sizeof(unsigned int));
  num_neigh[0] = p;
  
  neighbourhood = (unsigned int **)malloc((p+1)*sizeof(unsigned int*));
  for (i=1;i<=p;i++)
    neighbourhood[i] = findNeigh(v1,v2,i,p);

  //while not numbered all vertices and the graph is triangulated
  while ((n_numb<p) & (bTriangulated))
  {
    first = 0;
    for (i=1; i<=p; i++)
    {
      num_neigh[i] = 0;
      if (numbered[i]==0)
      {
        if (first==0)
          first = i;
        num_neigh[i] = 0;
        ind2 = neighbourhood[i];//findNeigh(v1,v2,i,p);

        //number of numbered neighbours for each vertex neighbour of i
        for (j=1; j<=ind2[0]; j++)
          if (numbered[ind2[j]] > 0)
            num_neigh[i]++;
      }
    }
    if (max(num_neigh)==0) //if there is no numbered neighbour
    {
      n_numb++;
      numbered[first] = n_numb; //just number the first on in the vector
    }
    else
    {
      v = whichMaxMCS(num_neigh); //the one with more numbered neighbours
      numbered[v] = max(numbered) + 1;
      n_numb++;
      
      ind2 = neighbourhood[v];//findNeigh(v1,v2,v,p);
      ind3 = which(numbered,0,false);
      ind4 = intersect(ind2,ind3);
      ind4[0]++;
      ind4[ind4[0]] = v;
      //free(ind2);
      free(ind3);
      if (!isComplete(v1,v2,ind4,p)) //tests if it complete
        bTriangulated = false;
      free(ind4);
    }
  }
  
  free(num_neigh);
  for (i=1;i<=p;i++) // note that the fisrt position was never used
    free(neighbourhood[i]);
  free(neighbourhood);

  if (!bTriangulated)
  {
    PROTECT(result = allocVector(INTSXP, 1));
    numPr++;
    INTEGER(result)[0] = 0;
  }
  else
  {
    PROTECT(result = allocVector(INTSXP, p));
    numPr++;
    for (i=1;i<=p;i++)
      INTEGER(result)[i-1] = numbered[i];
  }
  free(numbered);

  UNPROTECT(numPr);
  return(result);
}
