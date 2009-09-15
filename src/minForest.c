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

/******************************************************************************/
// Structure used in the binary tree constructed to order the p*(p-1)/2 possible
// edges according to its LR.
/******************************************************************************/
struct nodeBT
{
  double value; //LR for this edge
  int edge; //edge as the index in Combination(p,2)
            //e.g., Combination(4,2): index | vert1 | vert2
            //                            1 |     1 |     2
            //                            2 |     1 |     3
            //                            3 |     1 |     4
            //                            4 |     2 |     3
            //                            5 |     2 |     4
            //                            6 |     3 |     4
  unsigned short numP;
               //number of parameters for f12/(f1*f2)
               //1. if both are continuous: 5-(2+2)
               //2. if both are discrete: (k1-1)*(k2-1)
               //3. if continuous (f1) and discrete (f2):
               //   if homogeneous: (k2*1+1)+(k2-1)-(2+k2-1)
               //   if heterogeneous: (k2*2)+(k2-1)-(2+k2-1)
               //4. and the other way around
               //in 2, the marginals with zero are not counted:
               //   (k1-#zero.margins-1)*(k2-#zero.margins-1)
  struct nodeBT *greater; //left
  struct nodeBT *less; //right
};

/******************************************************************************/
// Like R's "max", but returns the maximum and the '2nd' maximum in a vector.
// In: d - SEXP
// Out: pointer to a vector of unsigned int (length 2).
/******************************************************************************/
unsigned int *maxS(SEXP d)
{
  unsigned int i,n;
  unsigned int *M;
  
  n = length(d);

  M = (unsigned int *)malloc(2*sizeof(unsigned int));
  M[0] = INTEGER(d)[0];
  if (n>1)
    M[1] = INTEGER(d)[1];
  else
    M[1] = M[0];
  for(i=1; i<n; i++)
    if (M[0] < INTEGER(d)[i])
    {
      M[1] = M[0];
      M[0] = INTEGER(d)[i];
    }
    else
      if (M[1] < INTEGER(d)[i])
        M[1] = INTEGER(d)[i];
  return M;
}

/******************************************************************************/
// Rounds a number.
// In: x - double, number to be rounded
//     places - unsigned char, number of decimal places
// Out: rounded number (double)
/******************************************************************************/
static double Round(double x, unsigned char places)
{
  double const shift = powf(10.0f, places);

  x *= shift;
  x = floorf( x + 0.5f );
  x /= shift;

  return(x);
}

/******************************************************************************/
// Creates a new node, with the parameters.
// In: v - double, the LR for the edge
//     edge - int, index of the edge
// Out: pointer to the new node
/******************************************************************************/
struct nodeBT* NewNodeBT(double v, unsigned int edge, unsigned int numP)
{
  struct nodeBT* nodeBT = malloc(sizeof *nodeBT);
  nodeBT->value = v;
  nodeBT->edge = edge;
  nodeBT->numP = numP;
  nodeBT->greater = NULL;
  nodeBT->less = NULL;
  return(nodeBT);
}

/******************************************************************************/
// Inserts the new node in the right position in the binary tree. It is a
// recursive function.
// In: node - pointer to nodeBT, previous checked node (root for a subtree)
//     v - double, LR for the edge
//     edge - int, index of the edge
// Out: pointer to nodeBT (if found the place, the new node, otherwise, the
//      next nde tha has to be tested).
/******************************************************************************/
static struct nodeBT* insertBT(struct nodeBT* node, double v, unsigned int edge,
                               unsigned int numP)
{
  if (node == NULL)
  {
    return(NewNodeBT(v,edge,numP)); //if the subtree is empty, return a new node
  }
  else //otherwise, recur down the tree, in the proper subtree
  {
    if (v <= node->value)
      node->less = insertBT(node->less,v,edge,numP);
    else
      node->greater = insertBT(node->greater,v,edge,numP);
    return(node); //return the (unchanged) node pointer
  }
}

/******************************************************************************/
// Given a index in a Combination(n,2), returns the 2 elements.
// In: x - int, the index
//     n - int, number of elements
// Out: pointer to int vector of 2 ([0] = vertex 1; [1] = vertex 2)
/******************************************************************************/
static unsigned int *getVertBT(unsigned int x, unsigned int n)
{
  unsigned int i, j, v1=0, v2=0, *a, *result;

  result = (unsigned int *)malloc(2*sizeof(unsigned int));
  i = n*(n-1)/2; //number of combinations
  if (x == i) //if the required is the last one
  {
    v1 = n-1;
    v2 = n;
  }
  else
  {
    if (x < i)
    {
      a = (unsigned int *)malloc((n-1+1)*sizeof(unsigned int));
      a[0] = n-1;
      for (i=1; i<=(n-1); i++)
        a[i] = i*n - i*(i+1)/2;

      j = 1;
      while ((j <= (n-1)) & (x >= (a[j]-(n-j-1))))
        j++;

      if (j > (n-1))
        v1 = n - 1;
      else
        v1 = j - 1;

      v2 = v1 + 1 + x - (a[v1]-(n-v1-1));
      free(a);
    }
  }
  result[0] = v1;
  result[1] = v2;
  return(result);
}

/******************************************************************************/
// Like R's "which", returns the indexes in d for elements equal (s=1) or not
// equal (s=0) k. But what is received is a SEXP to the dataset, and the values
// are searched in one column of the dataset, and is returned the row indexes.
// In: d - SEXP
//     k - int, the element of interest
//     s - int, if equal (1) or not equal (0)
//     c - int, column index (as a matrix, and not as a vector in the memory)
//         from 1 to ncol
//     ncol - int, number of columns in the matrix
// Out: pointer to an array of int with the indexes in d. If none was found,
//      then arr[0]=0, otherwise arr[0] has the number of elements in arr.
/******************************************************************************/
static unsigned int *whichSDsBT(SEXP d, unsigned int k, bool equal,
                                unsigned int c, unsigned int ncol)
{
  unsigned int i, j, nrow, *arr, *result;

  nrow = length(d)/ncol;
  //the maximum length result
  arr = (unsigned int *)malloc((nrow+1)*sizeof(unsigned int));
  j = 1; //number of elements added to arr
  if (equal)
    for (i=1; i<=nrow; i++) //tests all elements in d
    {
      if (INTEGER(d)[c*nrow+(i-1)]==INTEGER(d)[c*nrow+(i-1)]) //not a NA
        if (INTEGER(d)[c*nrow+(i-1)]==k) //looking for something equal
          arr[j++] = i-1;
    }
  else
    for (i=1; i<=nrow; i++) //tests all elements in d
      if (INTEGER(d)[c*nrow+(i-1)]==INTEGER(d)[c*nrow+(i-1)]) //not a NA
        if (INTEGER(d)[c*nrow+(i-1)]!=k) //looking for something equal
          arr[j++] = i-1;

  arr[0] = j-1; //as j starts in 1
  result = (unsigned int *)malloc(j*sizeof(unsigned int)); //precise lenght
  result[0] = arr[0];
  for (i=0;i<=arr[0];i++)
    result[i] = arr[i];
  free(arr);

  return(result);
}

/******************************************************************************/
// Like R's "c". Binds 2 vectors of int - c(A,B)
// In: A - pointer to a int vector
//     B - pointer to a int vector
// Out: pointer to a int vector with A[0]+B[0]+1 elements.
/******************************************************************************/
static unsigned int *bindBT(unsigned int *A, unsigned int *B)
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
// Vectorise the binary tree of LR values, releasing the nodes. Recursive.
// In: curr - pointer to nodeBT
//     stat - pointer to double vector, will contain the sorted (decreasing) LR
//     edges - pointer to int, will contain the respective egdes indexes
// Out: none (it updates the two vectors in "In")
/******************************************************************************/
static void vecBT(struct nodeBT *curr, double *stat, int *edges,
                  unsigned short *numP, unsigned int *errors)
{
  if (curr == NULL) return; //reached the end of a branch
  
  vecBT(curr->greater,stat,edges,numP,errors); //recur to a subtree

  edges[0]++;
  edges[edges[0]] = curr->edge; //add the current edge to the vector
  stat[0]++;
  stat[edges[0]] = curr->value;
  numP[edges[0]] = curr->numP;
  if (isnan(curr->value) || !isfinite(curr->value))
  {
    errors[0]++;
    errors[errors[0]] = curr->edge;
  }

  vecBT(curr->less,stat,edges,numP,errors); //recur to a subtree

  free(curr); //release the node (it has already been added)
}

/******************************************************************************/
// If there is any error in the marginal tests for the edges, a list with all
// edges with problems is build. Also free the memory.
// Recursive.
// In: curr - pointer to nodeBT
//     vert1 - pointer to int, will contain the first vertex of the edge
//     vert2 - pointer to int, will contain the second vertex of the edge
//     p - int, number of vertices
// Out: none (it updates the three vectors in "In")
/******************************************************************************/
static void errorEdges(struct nodeBT *curr, unsigned int *vert1,
                       unsigned int *vert2, unsigned int p)
{
  unsigned int *w;

  if (curr == NULL) return; //reached the end of a branch

  //recur to a subtree
  errorEdges(curr->greater,vert1,vert2,p);

  if ((curr->value != curr->value) || (!finite(curr->value))) // NaN or Inf
  {
    w = getVertBT(curr->edge,p); //get the vertices relative to that index
    vert1[0]++;
    vert2[0]++;
    vert1[vert1[0]] = w[0];
    vert2[vert2[0]] = w[1];
    free(w);
  }
  //recur to a subtree
  errorEdges(curr->less,vert1,vert2,p);
  free(curr); //release the node (it has already been visited)
}

/******************************************************************************/
// To clean the memory (the tree) when needed.
// Recursive.
// In: root - pointer to nodeBT
// Out: -
/******************************************************************************/
static void cleanMem(struct nodeBT *root)
{
  if (root == NULL) return; //reached the end of a branch
  cleanMem(root->greater); //recur to a subtree
  cleanMem(root->less); //recur to a subtree
  free(root); //release the node (it has already been visited)
}

/******************************************************************************/
// Build the real tree, based on the LR calculated, and release the nodes.
// Recursive.
// In: curr - pointer to nodeBT
//     vert1 - pointer to int, will contain the first vertex of the edge
//     vert2 - pointer to int, will contain the second vertex of the edge
//     LR - pointer to double vector, will contain the respective LR
//     p - int, number of vertices
//     comp - pointer int, indicates in which component the vertex is
//     type - pointer int, if the component is continuous (0), discrete (1), or
//            mixed (2)
//     varType - pointer SEXP, type of each variable
//     numP - pointer unsigned short, number of parameters in the edge
// Out: none (it updates the three vectors in "In")
/******************************************************************************/
static void bTree(struct nodeBT *curr, unsigned int *vert1, unsigned int *vert2,
                  double *LR, unsigned int p, unsigned int *comp,
//                  unsigned char *type, SEXP varType, unsigned short *numP,
                  unsigned char *type, SEXP numCat, unsigned short *numP,
                  unsigned int *errors)
{
  unsigned int *w, *ind1, *ind2, *ind3, i, min;
  unsigned char mix;
  bool mix1, mix2, forbiddenPath;
  
  if (curr == NULL) return; //reached the end of a branch

  //recur to a subtree
//  bTree(curr->greater,vert1,vert2,LR,p,comp,type,varType,numP,errors);
  bTree(curr->greater,vert1,vert2,LR,p,comp,type,numCat,numP,errors);
  if (isnan(curr->value) || !isfinite(curr->value))
  {
    errors[0]++;
    errors[errors[0]] = curr->edge;
  }
  else
    if (vert1[0]<(p-1)) //if the tree is not complete (still a forest)
    {
      if (curr->value > 0) //what is used here is -BIC, so, it has to be maximized
      {
        w = getVertBT(curr->edge,p); //get the vertices relative to that index
        // finds out if vertex x[2] is reacheable from x[1] in the current graph
        // if it is true, then a cycle would be generated by adding this edge

        if (comp[w[0]] != comp[w[1]])
        {
//          mix1 = (INTEGER(varType)[w[0]-1]==0) & (type[w[0]]==2);
          mix1 = (INTEGER(numCat)[w[0]-1]==0) & (type[w[0]]==2);
//          mix2 = (INTEGER(varType)[w[1]-1]==0) & (type[w[1]]==2);
          mix2 = (INTEGER(numCat)[w[1]-1]==0) & (type[w[1]]==2);
//          forbiddenPath = ((INTEGER(varType)[w[0]-1]==1) & mix2) |
//                          ((INTEGER(varType)[w[1]-1]==1) & mix1) |
//                          (mix1 & mix2);
          forbiddenPath = ((INTEGER(numCat)[w[0]-1]!=0) & mix2) |
                          ((INTEGER(numCat)[w[1]-1]!=0) & mix1) |
                          (mix1 & mix2);
          if (!forbiddenPath)
          {
            vert1[0]++;
            vert2[0]++;
            LR[0]++;
            vert1[vert1[0]] = w[0];
            vert2[vert2[0]] = w[1];
            numP[0]++;
            numP[vert1[0]] = curr->numP;
            LR[vert1[0]] = curr->value;
            ind1 = which(comp,comp[w[0]],1);
            ind2 = which(comp,comp[w[1]],1);
            ind3 = bindBT(ind1,ind2);
            min  = imin2(comp[w[0]],comp[w[1]]);
            //the new component has vertex of only one type
            if (type[w[0]]==type[w[1]])
              mix = type[w[0]];
            else
              mix = 2; //the new component is mixed
            for (i=1; i<=ind3[0]; i++)
            {
              comp[ind3[i]] = min;
              if (mix==0) type[ind3[i]] = 0;
              else if (mix==1) type[ind3[i]] = 1;
              else type[ind3[i]] = 2;
            }
            free(ind1);
            free(ind2);
            free(ind3);
          }
        }
        free(w);
      }
    }
  //recur to a subtree
//  bTree(curr->less,vert1,vert2,LR,p,comp,type,varType,numP,errors);
  bTree(curr->less,vert1,vert2,LR,p,comp,type,numCat,numP,errors);

  free(curr); //release the node (it has already been visited)
}

/******************************************************************************/
// Calculates the pairwise statistics.
// In: dataset - pointer SEXP (double), a matrix n by p (n observations in p
//               variables)
//     mean - pointer SEXP (double), vector with the p sample means
//     cov - pointer SEXP (double), p by p (covariance matrix)
//     varType - pointer SEXP (double), p; 0=variable is continuous; 1=discrete
//     numCat - pointer SEXP, (double), p; number of catagories in each discrete
//              variable; if the variable is continuous numCat=0;
//              If numCat=k, then the categories are: 1,2,...,k.
//     homog - if the model is homog (TRUE) or heterog (FALSE)
//     forbEdges - pointer SEXP (double), k; list of forbidden edges
//     stat - 0 (LR), 1 (BIC), or 2 (AIC)
// Out: pointer to a nodeBT
/******************************************************************************/
//struct nodeBT* calc(SEXP dataset, SEXP varType, SEXP numCat, bool homog,
//                    SEXP forbEdges, unsigned short stat, unsigned int *errors)
struct nodeBT* calc(SEXP dataset, SEXP numCat, bool homog,
                    SEXP forbEdges, unsigned short stat, unsigned int *errors)

{
  unsigned int i, j, p, n, ii, k, N, iii, jjj, kkk, identical, idX, idY;
  unsigned int continuous, discrete, *ind, *maxNumCat;
  double M1, M2, V1, V2, C12, x, X, Y;
  double v12, constK, Xi, Xi2;
  double value = 0;
  struct nodeBT *root = NULL;
  double *n1, *n2, *n12, aux;
  unsigned short numP = 0;
  double *auxXi, *auxXi2, *auxN;
  bool addEdge; //indicates if the edge is eligible, meaning that its not NA,...
  unsigned int numUsedObs; //to count the number of obs not NA in the mixed case

  p = length(numCat); //number of variables
  n = length(dataset)/p; //number of observations
  N = p*(p-1)/2; //number of possible edges

  maxNumCat = maxS(numCat);
  n1  = (double *)calloc((maxNumCat[0]+1),sizeof(double));
  n2  = (double *)calloc((maxNumCat[0]+1),sizeof(double));
  n12 = (double *)calloc((maxNumCat[0]*maxNumCat[1]+1),sizeof(double));
  auxXi = (double *)calloc((maxNumCat[0]),sizeof(double));
  auxXi2 = (double *)calloc((maxNumCat[0]),sizeof(double));
  auxN = (double *)calloc((maxNumCat[0]),sizeof(double));
  free(maxNumCat);

  k = 0;
  for (i=0;i<(p-1);i++)
  {
    for (j=(i+1);j<p;j++)
    {
      ind = whichSDsBT(forbEdges,(i-1+1)*p-(i-1+1)*(i+1)/2+(j+1)-(i+1),1,0,1);
      if (ind[0]==0)
      {
        free(ind);
        addEdge = true;
//        if ((INTEGER(varType)[i]==0) & (INTEGER(varType)[j]==0)) //both continuous
        if ((INTEGER(numCat)[i]==0) & (INTEGER(numCat)[j]==0)) //both continuous
        {
          numUsedObs = V1 = V2 = C12 = M1 = M2 = 0;
          identical = idX = idY = 0;
          X = REAL(dataset)[i*n+0];
          Y = REAL(dataset)[j*n+0];
          for (ii=0;ii<n;ii++)
          {
            if(!ISNA(REAL(dataset)[i*n+ii]) && !ISNA(REAL(dataset)[j*n+ii]))
            {
              numUsedObs++;
              V1  += REAL(dataset)[i*n+ii]*REAL(dataset)[i*n+ii];
              V2  += REAL(dataset)[j*n+ii]*REAL(dataset)[j*n+ii];
              M1  += REAL(dataset)[i*n+ii];
              M2  += REAL(dataset)[j*n+ii];
              C12 += REAL(dataset)[i*n+ii]*REAL(dataset)[j*n+ii];
              identical += REAL(dataset)[i*n+ii]==REAL(dataset)[j*n+ii];
              idX += X==REAL(dataset)[i*n+ii];
              idY += Y==REAL(dataset)[j*n+ii];
            }
          }
          V1 = sqrt((V1-M1*M1/numUsedObs)/(numUsedObs));
          V2 = sqrt((V2-M2*M2/numUsedObs)/(numUsedObs));
          C12 = (C12-M1*M2/numUsedObs)/(numUsedObs);
          C12 /= (V1*V2);
          M1 /= numUsedObs;
          M2 /= numUsedObs;
          v12 = 0;
          if ((identical==numUsedObs) || (idX==numUsedObs) || (idY==numUsedObs))
            value = NAN;
          else
            value = n*log(1-pow(C12,2))*(-1);
          numP = 5-(2+2);
        }
        else
//          if ((INTEGER(varType)[i]==1) & (INTEGER(varType)[j]==1)) //both discrete
          if ((INTEGER(numCat)[i]!=0) & (INTEGER(numCat)[j]!=0)) //both discrete
          {
            if ((INTEGER(numCat)[i]==1) || (INTEGER(numCat)[j]==1))
            {
              value = 0;
              numP = 0;
            }
            else
              if (INTEGER(numCat)[i]*INTEGER(numCat)[j] > 0)
              {
                v12 = 0;
                n12[0] = INTEGER(numCat)[i]*INTEGER(numCat)[j];
                n1[0]  = INTEGER(numCat)[i];
                n2[0]  = INTEGER(numCat)[j];
                for (ii=0; ii<n; ii++)
                {
                  //the index cannot be double
                  //the index in n12 is relative to the number of categories
                  iii = REAL(dataset)[i*n+ii];
                  jjj = REAL(dataset)[j*n+ii];
                  kkk = (REAL(dataset)[i*n+ii]-1)*INTEGER(numCat)[j]+REAL(dataset)[j*n+ii];
                  n12[kkk]++;
                  n1[iii]++;
                  n2[jjj]++;
                }
                for (ii=1; ii<=n12[0]; ii++)
                {
                  if (n12[ii]>0)
                  {
                    aux = n12[ii];
                    v12 = v12 + aux*log(aux/n);
                    n12[ii] = 0;
                  }
                  if (ii<=n1[0])
                  {
                    aux = n1[ii];
                    v12 = v12 - aux*log(aux/n);
                    n1[ii] = 0;
                  }
                  if (ii<=n2[0])
                  {
                    aux = n2[ii];
                    v12 = v12 - aux*log(aux/n);
                    n2[ii] = 0;
                  }
                }
                value = 2*v12;//2*(v12-(V[i+1]+V[j+1]));
                numP = (INTEGER(numCat)[i]-1)*(INTEGER(numCat)[j]-1);
              }
              else
                addEdge = false;
          }
          else //mixed
          {
            continuous = j;
            discrete = i;
//            if (INTEGER(varType)[i]==0)
            if (INTEGER(numCat)[i]==0)
            {
              continuous = i;
              discrete = j;
            }
            if (INTEGER(numCat)[discrete]==1)
            {
              value = 0;
              numP = 0;
            }
            else
            {
              v12 = 0;
              V1 = V2 = 0;
              numUsedObs = 0;
              Xi = Xi2 = 0;
              idX = 0;
              X = REAL(dataset)[continuous*n+0];
              for (ii=0; ii<n; ii++)
                if (!ISNA(REAL(dataset)[continuous*n+ii]))
                {
                  numUsedObs++;
                  iii = REAL(dataset)[discrete*n+ii] - 1;
                  x = REAL(dataset)[continuous*n+ii];
                  auxXi[iii] += x;
                  auxXi2[iii] += x*x;
                  Xi += x;
                  Xi2 += x*x;
                  auxN[iii]++;
                  idX += X==REAL(dataset)[continuous*n+ii];
                }
              if (idX==numUsedObs)
                value = NAN;
              else
              {
                V1 = (Xi2 - Xi*Xi/numUsedObs)/numUsedObs;
                if (homog)
                {
                  for (ii=0; ii<INTEGER(numCat)[discrete]; ii++)
                  {
                    V2 += auxXi2[ii] - auxXi[ii]*auxXi[ii]/auxN[ii];
                    auxXi[ii] = 0;
                    auxXi2[ii] = 0;
                    auxN[ii] = 0;
                  }
                 V2 /= numUsedObs;
                 v12 = n*log(V1/V2);
                }
                else
                {
                  v12 = n*log(V1);
                  for (ii=0; ii<INTEGER(numCat)[discrete]; ii++)
                  {
                    v12 -= auxN[ii]*log((auxXi2[ii] - auxXi[ii]*auxXi[ii]/auxN[ii])/auxN[ii]);
                    auxXi[ii] = 0;
                    auxXi2[ii] = 0;
                    auxN[ii] = 0;
                  }
                }
              }
              value = v12;
              numP = (INTEGER(numCat)[discrete]*(1+!homog)+homog+INTEGER(numCat)[discrete]-1) - (INTEGER(numCat)[discrete]+1);
            }
          }
        k = (i-1+1)*p-(i-1+1)*(i+1)/2+(j+1)-(i+1);
        switch(stat)
        {
          case 0:  constK = 0; break;
          case 1:  constK = -numP*log(n); break;
          case 2:  constK = -numP*2; break;
          default: constK = 0;
        }
        if (!addEdge)
          value = NAN;
        if (isnan(value) || !isfinite(value))
          errors[0]++;
        root = insertBT(root,value+constK,k,numP);
      }
      else
        free(ind);
    }
  }
  free(n1);
  free(n2);
  free(n12);
  free(auxXi);
  free(auxXi2);
  free(auxN);

  return(root);
}

/******************************************************************************/
// The main function. Just receives all the information from R, and calculates
// the LR.
// In: dataset - pointer SEXP (double), a matrix n by p (n observations in p
//               variables)
//     mean - pointer SEXP (double), vector with the p sample means
//     cov - pointer SEXP (double), p by p (covariance matrix)
//     varType - pointer SEXP (double), p; 0=variable is continuous; 1=discrete
//     numCat - pointer SEXP, (double), p; number of catagories in each discrete
//              variable; if the variable is continuous numCat=0;
//              If numCat=k, then the categories are: 1,2,...,k.
//     homog - if the model is homog (TRUE) or heterog (FALSE)
//     forbEdges - pointer SEXP (double), k; list of forbidden edges
//     stat - pointer SEXP, 0 (LR), 1 (BIC), 2 (AIC), or 3 (user values)
//     values - pointer SEXP, if stat>3
//     comp - pointer SEXP, to which component each node is assigned (if there
//            is no initial edge, then it's 1:p
//     compType - pointer SEXP, component type, 0 (cont), 1 (discr), 2 (mixed)
// Out: matrix (p-1) by 3: [,1] = first vertices, [,2] = second vertices;
//                         [,3] = LR.
/******************************************************************************/
//SEXP minForest(SEXP dataset, SEXP varType, SEXP numCat, SEXP HOMOG,
//               SEXP forbEdges, SEXP STAT,SEXP VALUES)
SEXP minForest(SEXP dataset, SEXP numCat, SEXP HOMOG,
               SEXP forbEdges, SEXP STAT,SEXP VALUES,
               SEXP COMP, SEXP COMPTYPE)
{
  unsigned int i, p, n, k, N, Nv, *errors, *w;
  unsigned int numPr = 0;
  struct nodeBT *root = NULL;
  unsigned short stat;
  bool homog = INTEGER(HOMOG)[0];

  stat = INTEGER(STAT)[0];
  errors = (unsigned int *)calloc(1,sizeof(unsigned int));
  p = length(numCat); //number of variables
  n = length(dataset)/p; //number of observations
  N = p*(p-1)/2; //number of possible edges

  if (stat<3) //means LR, AIC, or BIC
//    root = calc(dataset,varType,numCat,homog,forbEdges,stat,errors);
    root = calc(dataset,numCat,homog,forbEdges,stat,errors);
  else
  {
    Nv = length(VALUES)/2;
    for (i=0; i<Nv; i++)
      if (!ISNA(REAL(VALUES)[i]))
        root = insertBT(root,REAL(VALUES)[i],i+1,REAL(VALUES)[Nv+i]);
  }

  /////////////////////////////////////////////
  // "load" the binary tree with the LRs,
  // build a ordered array and free the tree
  ////////////////////////////////////////////
  unsigned int *vert1, *vert2;
  unsigned int *comp; //to keep track of which connected component each vertex
                      //belongs to
  double *LR;
  unsigned char *type; //to identify if the component is continuous (0),
                       //discrete (1) or mixed (2)
  unsigned short *numParam; //number of parameters
  SEXP tree; //the resulting tree or list of edges with problem in the marginal
  SEXP errorList;
  
  k = 0;
  //the "longest" possible tree (p-1 edges)
  vert1 = (unsigned int *)malloc(p*sizeof(unsigned int));
  vert1[0] = 0;
  vert2 = (unsigned int *)malloc(p*sizeof(unsigned int));
  vert2[0] = 0;
  numParam = (unsigned short *)malloc(p*sizeof(unsigned short));
  numParam[0] = 0;
  LR = malloc(p*sizeof *LR);
  LR[0] = 0;
  comp = (unsigned int *)malloc((p+1)*sizeof(unsigned int));
  comp[0] = p;
  type = (unsigned char *)malloc((p+1)*sizeof(unsigned char));
  type[0] = p;
  for (i=1; i<=p; i++)
  {
    //comp[i] = i;
    comp[i] = INTEGER(COMP)[i-1];
    //type[i] = 1*(INTEGER(numCat)[i-1]!=0);
    type[i] = INTEGER(COMPTYPE)[i-1];
  }
  i = errors[0];
  free(errors);
  errors = (unsigned int *)calloc(i+1,sizeof(unsigned int));
//  bTree(root,vert1,vert2,LR,p,comp,type,varType,numParam,errors);
  bTree(root,vert1,vert2,LR,p,comp,type,numCat,numParam,errors);

  free(comp);
  PROTECT(tree = allocMatrix(REALSXP, vert1[0], 4));
  numPr++;
  for (i=1;i<=vert1[0];i++)
  {
    REAL(tree)[0*vert1[0]+i-1] = vert1[i];
    REAL(tree)[1*vert1[0]+i-1] = vert2[i];
    REAL(tree)[2*vert1[0]+i-1] = LR[i];
    REAL(tree)[3*vert1[0]+i-1] = numParam[i];
  }
  free(vert1),
  free(vert2);
  free(LR);
  free(numParam);
  free(type);
  
  PROTECT(errorList = allocMatrix(INTSXP, errors[0], 2));
  numPr++;
  for (i=1;i<=errors[0];i++)
  {
    w = getVertBT(errors[i],p);
    INTEGER(errorList)[0*errors[0]+i-1] = w[0];
    INTEGER(errorList)[1*errors[0]+i-1] = w[1];
    free(w);
  }
  free(errors);

  SEXP listNames; //the names (in R) of the result components
  char *names[2] = {"tree", "errors"};
  PROTECT(listNames = allocVector(STRSXP,2));
  numPr++;
  SET_STRING_ELT(listNames,0,mkChar(names[0]));
  SET_STRING_ELT(listNames,1,mkChar(names[1]));

  SEXP list;
  PROTECT(list = allocVector(VECSXP,2));
  numPr++;

  SET_VECTOR_ELT(list, 0, tree);
  SET_VECTOR_ELT(list, 1, errorList);
  setAttrib(list, R_NamesSymbol, listNames);

  UNPROTECT(numPr);
  return(list);
}


/******************************************************************************/
// Just a test function, returning the calculated pairwise statistics.
// In: dataset - pointer SEXP (double), a matrix n by p (n observations in p
//               variables)
//     mean - pointer SEXP (double), vector with the p sample means
//     cov - pointer SEXP (double), p by p (covariance matrix)
//     varType - pointer SEXP (double), p; 0=variable is continuous; 1=discrete
//     numCat - pointer SEXP, (double), p; number of catagories in each discrete
//              variable; if the variable is continuous numCat=0;
//              If numCat=k, then the categories are: 1,2,...,k.
//     homog - if the model is homog (TRUE) or heterog (FALSE)
// Out: matrix (p-1) by 3: [,1] = first vertices, [,2] = second vertices;
//                         [,3] = LR.
/******************************************************************************/
//SEXP calcStat(SEXP dataset, SEXP varType, SEXP numCat, SEXP HOMOG,
//              SEXP forbEdges, SEXP STAT,SEXP VALUES)
SEXP calcStat(SEXP dataset, SEXP numCat, SEXP HOMOG,
              SEXP forbEdges, SEXP STAT,SEXP VALUES)
{
  unsigned int i, p, n, k, N, Nv, *w;
  unsigned int numPr = 0;
  struct nodeBT *root = NULL;
  unsigned short measure;
  bool homog = INTEGER(HOMOG)[0];
  unsigned int *ind, *errors;

  measure = INTEGER(STAT)[0];

  p = length(numCat); //number of variables
  n = length(dataset)/p; //number of observations
  N = p*(p-1)/2; //number of possible edges
  errors = (unsigned int *)calloc(1,sizeof(unsigned int));

  if (measure<3) //means LR, AIC, or BIC
//    root = calc(dataset,varType,numCat,homog,forbEdges,measure,errors);
    root = calc(dataset,numCat,homog,forbEdges,measure,errors);
  else
  {
    Nv = length(VALUES)/2;
    for (i=0; i<Nv; i++)
      if (!ISNA(REAL(VALUES)[i]))
        root = insertBT(root,REAL(VALUES)[i],i+1,REAL(VALUES)[Nv+i]);
  }

  i = errors[0];
  free(errors);
  errors = (unsigned int *)calloc(i+1,sizeof(unsigned int));

  if (INTEGER(forbEdges)[0]!=0) //there is at least one forbidden
    N = N - length(forbEdges);
  SEXP stat; //the resulting tree
  SEXP errorList;
  
  PROTECT(stat = allocMatrix(REALSXP, N, 4));
  numPr++;
  double *v;
  int *edges;
  unsigned short *numParam;
  k = 1;
  v = (double *)malloc((N+1)*sizeof(double));
  edges = (int *)malloc((N+1)*sizeof(int));
  numParam = (unsigned short *)malloc((N+1)*sizeof(unsigned short));
  v[0] = 0;
  edges[0] = 0;
  numParam[0] = 0;

  vecBT(root,v,edges,numParam,errors);

  ind = (unsigned int *)malloc((2)*sizeof(unsigned int));
  for (i=1;i<=edges[0];i++)
  {
    ind = getVertBT(edges[i],p);
    REAL(stat)[0*N+i-1] = ind[0];
    REAL(stat)[1*N+i-1] = ind[1];
    REAL(stat)[2*N+i-1] = v[i];
    REAL(stat)[3*N+i-1] = numParam[i];
    free(ind);
  }
  free(v);
  free(edges);
  free(numParam);

  PROTECT(errorList = allocMatrix(INTSXP, errors[0], 2));
  numPr++;
  for (i=1;i<=errors[0];i++)
  {
    w = getVertBT(errors[i],p);
    INTEGER(errorList)[0*errors[0]+i-1] = w[0];
    INTEGER(errorList)[1*errors[0]+i-1] = w[1];
    free(w);
  }
  free(errors);

  SEXP listNames; //the names (in R) of the result components
  char *names[2] = {"stat", "errors"};
  PROTECT(listNames = allocVector(STRSXP,2));
  numPr++;
  SET_STRING_ELT(listNames,0,mkChar(names[0]));
  SET_STRING_ELT(listNames,1,mkChar(names[1]));

  SEXP list;
  PROTECT(list = allocVector(VECSXP,2));
  numPr++;

  SET_VECTOR_ELT(list, 0, stat);
  SET_VECTOR_ELT(list, 1, errorList);
  setAttrib(list, R_NamesSymbol, listNames);

  UNPROTECT(numPr);
  return(list);
}
