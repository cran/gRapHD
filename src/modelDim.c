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
struct node
{
  unsigned int *a; //generator
  char ex; //exponent
  struct node *next; //pointer to the next node
};

/******************************************************************************/
unsigned int dd(struct node* A, unsigned int p, char* numCat)
{
  unsigned int x, nc;
  unsigned int i, *J, *aux, *aux1;
  char ex;
  struct node *B, *currA, *currB, *prevB;
  bool found, empty, firstISempty;
  
  if (A->next == NULL)
  {
    if (A->a[0] == 0)
    {
      free(A->a);
      free(A);
      return(1);
    }
    else
    {
      x = 1;
      nc = 0;
      for (i=1;i<=A->a[0];i++)
        if (numCat[A->a[i]] == 0)
          nc++;
        else
          x *= numCat[A->a[i]];
      if (nc == 0)
        ex = 1;
      else
        if (nc == 1)
          ex = 1 + A->ex;
        else
          if ((nc == 2) && (A->ex == 1))
            ex = 4;
          else
            ex = (nc+1)*(nc+2)/2;
      free(A->a);
      free(A);
      return(x*ex);
    }
  }
  else
  {
    J = A->a;
    ex = A->ex;
    currA = A;
    A = A->next;
    free(currA);
    
    

    B = NULL;
    currA = A;
    firstISempty = false;
    while (currA != NULL)
    {
      aux = intersect(J,currA->a);
      empty = (aux[0]==0);
      currB = B;
      prevB = B;
      found = false;
      while (!found && (currB != NULL))
      {
        //found = identical(aux,currB->a);
        aux1 = intersect(aux,currB->a);
        //considering that there is no repetitions
        found = ((aux1[0]==aux[0]) && (aux1[0]==currB->a[0]));
        free(aux1);
        prevB = currB;
        currB = currB->next;
      }
      if (!found)
      {
        if (((B == NULL) && empty) | (!empty))
        {
          currB = malloc(sizeof *currB);
          currB->next = NULL;
          currB->a = aux;
          if (ex < currA->ex)
            currB->ex = ex;
          else
            currB->ex = currA->ex;
          if (B == NULL)
          {
            firstISempty = empty;
            B = currB;
          }
          else
            prevB->next = currB;
        }
      }
      else
        free(aux);
      currA = currA->next;
    }
    if ((firstISempty) && (B->next!=NULL))
    {
      currB = B;
      B = B->next;
      free(currB->a);
      free(currB);
    }


    x = 1;
    nc = 0;
    for (i=1;i<=J[0];i++)
      if (numCat[J[i]] == 0)
        nc++;
      else
        x *= numCat[J[i]];
    if (nc == 0)
      ex = 1;
    else
      if (nc == 1)
        ex = 1 + ex;
      else
        if ((nc == 2) && (ex == 1))
          ex = 4;
        else
          ex = (nc+1)*(nc+2)/2;
    free(J);
    return(x*ex + dd(A,p,numCat) - dd(B,p,numCat));
  }
}

/******************************************************************************/
SEXP modelDim(SEXP list, SEXP expon, SEXP NUMCAT)
{
  unsigned int N, i, j, p;
  struct node *A, *curr, *prev;
  SEXP aux, result;
  char *numCat;

  N = length(expon);
  A = curr = prev = NULL;

  for (i=0;i<N;i++)
  {
    aux = VECTOR_ELT(list, i);
    curr = malloc(sizeof *curr);
    curr->a = (unsigned int*)malloc((length(aux)+1)*sizeof(unsigned int));
    curr->a[0] = length(aux);
    for (j=0;j<curr->a[0];j++)
      curr->a[j+1] = INTEGER(aux)[j];
    curr->ex = INTEGER(expon)[i];
    curr->next = NULL;
    if (A==NULL)
      A = prev = curr;
    else
    {
      prev->next = curr;
      prev = curr;
    }
  }
  
  p = length(NUMCAT);
  numCat = (char*)malloc((p+1)*sizeof(char));
  for (i=0;i<p;i++)
    numCat[i+1] = INTEGER(NUMCAT)[i];
  
  PROTECT(result = allocVector(INTSXP,1));
  INTEGER(result)[0] = dd(A,p,numCat);
  free(numCat);

  UNPROTECT(1);
  return(result);
}
