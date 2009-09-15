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
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "comfunc.h"

unsigned int whichMinSP(unsigned int *d,unsigned int *q)
{
  unsigned int i;
  unsigned int m;
  unsigned int result=0;

  i = 1;
  while((i<=q[0]) && (q[i]==0))
    i++;
  m = d[q[i+1]]; //set the minimum as the first element
  result = q[i+1]; //the index of the first "minimum"
  for(i=1; i<=d[0]; i++)
  {
    if (q[i] > 0)
      if (m >= d[q[i]])
      {
        m = d[i];
        result = q[i];
      }
  }

  return(result);
}

SEXP shortPath(SEXP v1, SEXP v2, SEXP V, SEXP P)
{
  unsigned int *distV, *Q, *ind;
  unsigned int v, p, nQ, u, j, alt, i, x;
  unsigned int **neighbourhood;
  
  p = INTEGER(P)[0];
  v = INTEGER(V)[0];
  distV = (unsigned int *)malloc((p+1)*sizeof(unsigned int)); //distance from v to i (dist[i])
  Q = (unsigned int *)malloc((p+1)*sizeof(unsigned int));
  for (i=1;i<=p;i++)
  {
    distV[i] = p; //miximum distance + 1
    Q[i] = i; //set of vertices
  }
  distV[0] = p; //number of elements in the array
  Q[0] = p; //idem
  distV[v] = 0; //distance from v to v

  neighbourhood = (unsigned int **)malloc((p+1)*sizeof(unsigned int*));
  for (i=1;i<=p;i++) neighbourhood[i] = findNeigh(v1,v2,i,p);
  
  nQ = p; //number of elements in Q not yet visited
  while (nQ > 0)
  {
    u = whichMinSP(distV,Q); //the vertex with minimum distance from v
    x = Q[u];
    Q[u] = 0;
    nQ--;
    ind = neighbourhood[x];//findNeigh(v1,v2,x,p);
    for (j=1;j<=ind[0];j++)
    {
      alt = distV[x] + 1;
      if (alt < distV[ind[j]])
        distV[ind[j]] = alt;
    }
    //free(ind);
  }
  for (i=1;i<=p;i++) free(neighbourhood[i]);
  free(neighbourhood);

  SEXP result;
  PROTECT(result = allocVector(INTSXP,p));
  for (i=0;i<p;i++)
    INTEGER(result)[i] = distV[i+1];
  free(distV);
  free(Q);

  UNPROTECT(1);
  return(result);
}

SEXP spAll(SEXP v1, SEXP v2, SEXP P)
{
  unsigned int *distV, *Q, *ind;
  unsigned int v, p, nQ, u, j, alt, i, x;
  unsigned int **neighbourhood;

  p = INTEGER(P)[0];

  SEXP result;
  PROTECT(result = allocMatrix(INTSXP,p,p));
  distV = (unsigned int *)malloc((p+1)*sizeof(unsigned int)); //distance from v to i (dist[i])
  Q = (unsigned int *)malloc((p+1)*sizeof(unsigned int));
  neighbourhood = (unsigned int **)malloc((p+1)*sizeof(unsigned int*));
  for (i=1;i<=p;i++) neighbourhood[i] = findNeigh(v1,v2,i,p);

  for (v=1;v<=p;v++)
  {
    for (i=1;i<=p;i++)
    {
      distV[i] = p; //miximum distance + 1
      Q[i] = i; //set of vertices
    }
    distV[0] = Q[0] = p; //number of elements in the array
    distV[v] = 0; //distance from v to v

    nQ = p; //number of elements in Q not yet visited
    while (nQ > 0)
    {
      u = whichMinSP(distV,Q); //the vertex with minimum distance from v
      x = Q[u];
      Q[u] = 0;
      nQ--;
      ind = neighbourhood[x];//findNeigh(v1,v2,x,p);
      for (j=1;j<=ind[0];j++)
      {
        alt = distV[x] + 1;
        if (alt < distV[ind[j]])
          distV[ind[j]] = alt;
      }
    }
    for (i=0;i<p;i++)
      INTEGER(result)[i*p+v-1] = distV[i+1];
  }

  for (i=1;i<=p;i++) free(neighbourhood[i]);
  free(neighbourhood);

  free(distV);
  free(Q);

  UNPROTECT(1);
  return(result);
}

