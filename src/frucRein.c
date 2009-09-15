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

/******************************************************************************/
// Find the coordenates of the vertices.
// Fruchterman, T. M. J., & Reingold, E. M. (1991). Graph Drawing by
// Force-Directed Placement. Software: Practice and Experience, 21(11).
//
// Reference: Csardi G, Nepusz T: The igraph software package for complex
// network research, InterJournal, Complex Systems 1695. 2006.
// http://igraph.sf.net
//
// In: edges - pointer SEXP (double), a matrix numVert by 2
//     p - pointer SEXP (double), number of vertices
//     NUMITER - pointer SEXP (double), number of iterations
//     POS - pointer SEXP (double), inicial coordenates
// Out: matrix numVert by 2 with the updated coordenates
/******************************************************************************/
SEXP frucRein(SEXP edges, SEXP p, SEXP NUMITER, SEXP POS)
{
  unsigned int i, j, k, numVert, numEd, numIter, v1, v2;
  double *disp;
  double *pos;
  double area, K, temp, Dx, Dy, delta, repF, atrF, x, y, z;
  short unsigned numPr = 0;

  numIter = INTEGER(NUMITER)[0];
  numEd = length(edges)/2;
  numVert = INTEGER(p)[0];
  disp = (double *)calloc(2*numVert,sizeof(double));
  pos  = (double *)calloc(2*numVert,sizeof(double));
  area = numVert*numVert;
  for (i=0; i<numVert; i++)
  {
    pos[i] = REAL(POS)[i];//runif(-area/4,area/4);
    pos[i+numVert] = REAL(POS)[i+numVert];//runif(-area/4,area/4);
  }
  K = sqrt(area/(double)numVert);
  
  for (i=numIter; i>0; i--)
  {
    x = numVert;
    y = i;
    z = numIter;
    temp = x*pow(y/z,1.5);
    //temp = numVert*pow(i/(double)numIter,1.5);
    for (j=0; j<numVert; j++)
    {
      disp[j] = 0;
      disp[j+numVert] = 0;
    }
    for (j=0; j<numVert; j++)
      for (k=j+1; k<numVert; k++)
      {
        Dx = pos[j] - pos[k];
        Dy = pos[j+numVert] - pos[k+numVert];
        delta = sqrt(Dx*Dx+Dy*Dy);
        Dx /= delta;
        Dy /= delta;
        repF = K*K * (1/delta - delta*delta/(area*numVert));
        disp[j] += Dx*repF;
        disp[k] -= Dx*repF;
        disp[j+numVert] += Dy*repF;
        disp[k+numVert] -= Dy*repF;
      }
    for (j=0; j<numEd; j++)
    {
      v1 = INTEGER(edges)[j]-1;
      v2 = INTEGER(edges)[j+numEd]-1;
      Dx = pos[v1] - pos[v2];
      Dy = pos[v1+numVert] - pos[v2+numVert];
      delta = sqrt(Dx*Dx+Dy*Dy);
      if (delta != 0)
      {
        Dx /= delta;
        Dy /= delta;
      }
      atrF = delta*delta/K;
      disp[v1] -= Dx*atrF;
      disp[v2] += Dx*atrF;
      disp[v1+numVert] -= Dy*atrF;
      disp[v2+numVert] += Dy*atrF;
    }
    for (j=0; j<numVert; j++)
    {
      delta = sqrt(disp[j]*disp[j] + disp[j+numVert]*disp[j+numVert]);
      if (delta > temp)
      {
        delta = temp/delta;
        disp[j] *= delta;
        disp[j+numVert] *= delta;
      }
      pos[j] += disp[j];
      pos[j+numVert] += disp[j+numVert];
    }
  }
  free(disp);

  SEXP result;
  
  PROTECT(result = allocMatrix(REALSXP, numVert, 2));
  numPr++;
  for (i=0; i<numVert; i++)
  {
    REAL(result)[0*numVert+i] = pos[i];
    REAL(result)[1*numVert+i] = pos[i+numVert];
  }

  UNPROTECT(numPr);
  return(result);
}
