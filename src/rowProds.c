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
SEXP rowProds(SEXP matrix, SEXP NROW, SEXP NCOL, SEXP NARM)
{
  int nrow, ncol, i, j;
  int numPr = 0;
  double prod;
  bool naRm = INTEGER(NARM)[0];
  SEXP result;
  
  nrow = INTEGER(NROW)[0];
  ncol = INTEGER(NCOL)[0];
  PROTECT(result = allocVector(REALSXP, nrow));
  numPr++;

  for (i=0; i<nrow; i++)
  {
    prod = 1;
    j = 0;
    while ((j<ncol) & (prod!=0))
    {
      if ((!ISNA(REAL(matrix)[j*nrow+i]) | (ISNA(REAL(matrix)[j*nrow+i]) & !naRm)))
        prod = prod*REAL(matrix)[j*nrow+i];
      j++;
    }
    REAL(result)[i] = prod;
  }

  UNPROTECT(numPr);
  return(result);
}
