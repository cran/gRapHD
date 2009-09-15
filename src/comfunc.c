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
#include <R_ext/Utils.h>
#include "comfunc.h"

/******************************************************************************/
// called from: findEd.c; mcs.c; buildF.c
//
// Like R's "which", returns the indexes in d for elements equal (s=1) or not
// equal (s=0) k.
// In: d - pointer to unsigned int vector
//     k - unsigned int, the element of interest
//     equal - bool, if equal (true) or not equal (false)
// Out: pointer to an array of int with the indexes in d. If none was found,
//      then arr[0]=0, otherwise arr[0] has the number of elements in arr.
/******************************************************************************/
unsigned int *which(unsigned int *d, unsigned int k, bool equal)
{
  unsigned int i, j, *arr, *result;

  //the maximum possible result
  arr = (unsigned int *)malloc((d[0]+1)*sizeof(unsigned int));
  j = 1; //number of elements added to arr

  if (equal)
    for (i=1; i<=d[0]; i++) //tests all elements in d
    {
      if (d[i]==k) arr[j++] = i; //looking for something equal
    }
  else
    for (i=1; i<=d[0]; i++) //tests all elements in d
      if (d[i]!=k) arr[j++] = i; //looking for something equal
  arr[0] = j-1; //as j starts in 1
  result = (unsigned int *)malloc(j*sizeof(unsigned int)); //precise lenght
  for (i=0;i<=arr[0];i++) result[i] = arr[i];
  free(arr);

  return(result);
}

/******************************************************************************/
// called from: findEd.c; mcs.c; dfs.c
//
// Like R's "which", returns the indexes in d for elements equal (s=1) or not
// equal (s=0) k.
// In: d - pointer SEXP (numeric)
//     k - unsigned int, the element of interest
//     equal - bool, if equal (true) or not equal (false)
// Out: pointer to an array of int with the indexes in d. If none was found,
//      then arr[0]=0, otherwise arr[0] has the number of elements in arr.
/******************************************************************************/
unsigned int *whichS(SEXP d, unsigned int k, bool equal)
{
  unsigned int i, p, j, *arr, *result;

  p = length(d);
  //the maximum possible result
  arr = (unsigned int *)malloc((p+1)*sizeof(unsigned int));
  j = 1; //number of elements added to arr

  if (equal)
    for (i=0; i<p; i++) //tests all elements in d
    {
      if (INTEGER(d)[i]==k) arr[j++] = i;//looking for something equal
    }
  else
    for (i=0; i<p; i++) //tests all elements in d
      if (INTEGER(d)[i]!=k) arr[j++] = i;//looking for something equal

  arr[0] = j-1; //as starts in 1
  result = (unsigned int *)malloc(j*sizeof(unsigned int)); //precise lenght
  for (i=0;i<=arr[0];i++) result[i] = arr[i];
  free(arr);

  return(result);
}

/******************************************************************************/
// called from: mcs.c; findEd.c; dfs.c
//
// Find the neighbours of a vertex.
// In: v1 - pointer SEXP (int), with the first vertex of the edge
//     v2 - pointer SEXP (int), with the second vertex of the edge
//     v  - unsigned int, the vertex which neighbours are required
//     p  - unsigned int, total number of vertices
// Out: pointer to unsigned int vector with the neighbours (vec[0]=number of
//      elements).
/******************************************************************************/
unsigned int *findNeigh(SEXP v1, SEXP v2, unsigned int v, unsigned int p)
{
  unsigned int i, *ind1, *ind2, *result;

  //this if is in case a empty edges list is used (included because of the
  //option of joining desconnected components
  if (ISNA(INTEGER(v1)[0]) | ISNA(INTEGER(v2)[0]))
  {
    result = (unsigned int *)malloc(1*sizeof(unsigned int));
    result[0] = 0;
    return(result);
  }

  ind1 = whichS(v1,v,true); //edges in which v is the lowest vertex
  ind2 = whichS(v2,v,true); //edges in which v is the highest vertex
  result = (unsigned int *)malloc((ind1[0]+ind2[0]+1)*sizeof(unsigned int));
  for (i=1;i<=ind1[0];i++)
    result[i] = INTEGER(v2)[ind1[i]]; //gets the highest when v is the lowest
  for (i=1;i<=ind2[0];i++)
    result[ind1[0]+i] = INTEGER(v1)[ind2[i]]; //opposite
  result[0] = ind1[0] + ind2[0];
  free(ind1);
  free(ind2);
  return result;
}

/******************************************************************************/
// called from: mcs.c; findEd.c
//
// Like R's "max", returns the maximum in a vector.
// In: d - pointer to unsigned int vector
// Out: unsigned int, the maximum in d.
/******************************************************************************/
unsigned int max(unsigned int *d)
{
  unsigned int i, m;

  m = d[1];
  for(i=2; i<=d[0]; i++)
    if (m < d[i]) m = d[i];
  return m;
}

/******************************************************************************/
// called from: mcs.c; findEd.c
//
// Like R's "setdiff", returns the vector A\B.
// In: A - pointer to unsigned int vector
//     B - pointer to unsigned int vector
// Out: pointer to int vector, with A\B (vec[0]=number of elements in vec).
/******************************************************************************/
unsigned int *setDiff(unsigned int *A, unsigned int *B)
{
  unsigned int *R, *result, i, j, k = 1;
  bool del;

  //the maximum possible result
  R = (unsigned int *)malloc((A[0]+1)*sizeof(unsigned int));
  k = 1;

  for (i=1; i<=A[0]; i++)
  {
    del = false;
    j = 1;
    while (!del && (j<=B[0]))
    {
      del = A[i]==B[j];
      j++;
    }
    if (!del)
    {
      R[k] = A[i];
      k++;
    }
  }

  R[0] = k-1;
  result = (unsigned int *)malloc((k+1+1)*sizeof(unsigned int)); //precise length
  for (i=0;i<=R[0];i++) result[i] = R[i];
  free(R);

  return(result);
}

/******************************************************************************/
// called from: mcs.c; findEd.c
//
// Like R's "setequal", returns TRUE if A=B, FALSE otherwise.
// In: A - pointer to unsigned int vector
//     B - pointer to unsigned int vector
// Out: boolean.
/******************************************************************************/
bool setEqual(unsigned int *A, unsigned int *B)
{
  bool stop, isE = true;
  unsigned int j, k = 1;

  if ((A==NULL) | (B==NULL)) //if at least one is NULL
    isE = false;
  else
    if (A[0] != B[0]) //if the number of elements is different
      isE = 0;
    else
    {
      while (isE & (k<=A[0])) //for all elements in one edge, seaches the other
      {
        stop = false;
        j = 1;
        while (!stop && (j<=B[0]))
        {
          stop = A[k]==B[j]; // no repetition allowed
          j++;
        }
        isE = stop; //if stoped because found A[k]==B[j]
        k++;
      }
    }

  return(isE);
}

/******************************************************************************/
// called from: mcs.c; findEd.c
//
// Like R's "intersect", returns the vector with the intersection of A and B.
// In: A - pointer to unsigned int vector
//     B - pointer to unsigned int vector
// Out: pointer to unsigned int vector (vec[0]=number of elements in vec).
/******************************************************************************/
unsigned int *intersect(unsigned int *A, unsigned int *B)
{
  unsigned int *R, nR, *C, *D, i, k, *result, j;
  unsigned int nA = A[0], nB = B[0];
  bool inc;

  if (nA > nB) //selects the shortest vector to run through
  {
    nR = nB;
    C = B;
    D = A;
  }
  else
  {
    nR = nA;
    C = A;
    D = B;
  }
  //maximum possible lenght
  R = (unsigned int *)malloc((nR+1)*sizeof(unsigned int));
  k = 1;
  for (i=1; i<=nR; i++)
  {
    inc = false;
    j = 1;
    while (!inc && (j<=D[0]))
    {
      inc = C[i]==D[j];
      j++;
    }
    if (inc)
    {
      R[k] = C[i];
      k++;
    }
  }

  R[0] = k-1;
  //precise lenght + 1
  result = (unsigned int *)malloc((k+1+1)*sizeof(unsigned int));
  for (i=0;i<=R[0];i++) result[i] = R[i];
  free(R); //cannot free C and D because they are just references to the real
           //vectors!!!

  return(result);
}

/******************************************************************************/
// called from: mcs.c; findEd.c
//
// Tests if a subgraph is complete.
// In: v1 - pointer SEXP (double, but considered ind), with the first vertex of
//          the edge
//     v2 - pointer SEXP (double, but considered ind), with the second vertex of
//          the edge
//     vs - vector unsigned int, subset of vertices
//     p  - unsigned int, total number of vertices
// Out: boolean (true=complete; false=not complete).
/******************************************************************************/
bool isComplete(SEXP v1, SEXP v2, unsigned int *vs, unsigned int p)
{
  bool isC = true;
  unsigned int k = 1;
  unsigned int *fn, *aux, *aux1;

  while ((isC==1) & (k<=vs[0])) //for each vertex in the subgraph
  {
    fn = findNeigh(v1,v2,vs[k],p); //find all neighbours
    aux1 = setDiff(fn,vs); //all neighbours that are not in the subgraph
    aux = setDiff(fn,aux1); //all neighbours that are in the subgraph
    aux[0]++;
    aux[aux[0]] = vs[k]; //to the neighbours in the subgraph, add the vertex
    isC = setEqual(aux,vs); //if this and the subgraph are equal, isC=true
    free(aux);
    free(aux1);
    free(fn);
    k++;
  }
  return(isC);
}
