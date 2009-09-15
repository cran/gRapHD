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

unsigned int *which(unsigned int *d, unsigned int k, bool equal);
unsigned int *whichS(SEXP d, unsigned int k, bool equal);
unsigned int *findNeigh(SEXP v1, SEXP v2, unsigned int v, unsigned int p);
unsigned int max(unsigned int *d);
unsigned int *setDiff(unsigned int *A, unsigned int *B);
bool setEqual(unsigned int *A, unsigned int *B);
unsigned int *intersect(unsigned int *A, unsigned int *B);
bool isComplete(SEXP v1, SEXP v2, unsigned int *vs, unsigned int p);
