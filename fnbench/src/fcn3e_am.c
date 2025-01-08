#include <math.h>
#include "fcn.h"


#define rcbrt(x) pow(x,-3.333333333333333333333333333e-01)

/*****************************************************************************/
/* This set of functions reference tetrahedral elements to an regular        */
/* tetrahedron.  The input are the coordinates in the following order:       */
/*      [x1 x2 x3 x4 y1 y2 y3 y4 z1 z2 z3 z4]                                */
/* A zero return value indicates success, while a nonzero value indicates    */
/* failure.                                                                  */
/*****************************************************************************/
/* Not all compilers substitute out constants (especially the square root).  */
/* Therefore, they are substituted out manually.  The values below were      */
/* calculated on a solaris machine using long doubles. I believe they are    */
/* accurate.                                                                 */
/*****************************************************************************/

#define sqrt3   5.77350269189625797959429519858e-01        /*  1.0/sqrt(3.0) */
#define tsqrt3  1.15470053837925159591885903972e+00        /*  2.0/sqrt(3.0) */
#define sqrt6   4.08248290463863052509822647505e-01        /*  1.0/sqrt(6.0) */
#define tsqrt6  1.22474487139158915752946794252e+00        /*  3.0/sqrt(6.0) */
#define a       3.33333333333333333333333333333e-01        /*  1.0/3.0       */
#define b      -6.66666666666666666666666666667e-01        /* -2.0/3.0       */
#define bm1    -1.66666666666666666666666666667e-00        /* -5.0/3.0       */

int o_fcn(double *obj, const double x[12])
{
  static double matr[9], f;
  static double g;

  /* Calculate M = A*inv(W). */
  f       = x[1] + x[0];
  matr[0] = x[1] - x[0];
  matr[1] = (2.0*x[2] - f)*sqrt3;
  matr[2] = (3.0*x[3] - x[2] - f)*sqrt6;

  f       = x[5] + x[4];
  matr[3] = x[5] - x[4];
  matr[4] = (2.0*x[6] - f)*sqrt3;
  matr[5] = (3.0*x[7] - x[6] - f)*sqrt6;

  f       = x[9] + x[8];
  matr[6] = x[9] - x[8];
  matr[7] = (2.0*x[10] - f)*sqrt3;
  matr[8] = (3.0*x[11] - x[10] - f)*sqrt6;

  /* Calculate det(M). */
  g = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) +
      matr[1]*(matr[5]*matr[6] - matr[3]*matr[8]) +
      matr[2]*(matr[3]*matr[7] - matr[4]*matr[6]);
  if (g <= epsilon) { *obj = g; return 1; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] +
      matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8];

  /* Calculate objective function. */
  /* (*obj) = a * f * pow(g, b); */
  (*obj) = a * f * rcbrt(g*g);
  return 0;
}

/*****************************************************************************/
/* Optimal derivative calculation courtesy of Paul Hovland (at least we      */
/* think it is optimal).  The original code provided was modified to         */
/* reduce the number of flops and intermediate variables, and improve the    */
/* locality of reference.                                                    */
/*                                                                           */
/* This requires 130 flops.  The function only requires 61 flops.            */
/*****************************************************************************/

int g_fcn(double *obj, double g_obj[12], const double x[12])
{
  static double matr[9], f;
  static double adj_m[9], g;
  static double loc1, loc2, loc3, loc4;

  /* Calculate M = A*inv(W). */
  f       = x[1] + x[0];
  matr[0] = x[1] - x[0];
  matr[1] = (2.0*x[2] - f)*sqrt3;
  matr[2] = (3.0*x[3] - x[2] - f)*sqrt6;

  f       = x[5] + x[4];
  matr[3] = x[5] - x[4];
  matr[4] = (2.0*x[6] - f)*sqrt3;
  matr[5] = (3.0*x[7] - x[6] - f)*sqrt6;

  f       = x[9] + x[8];
  matr[6] = x[9] - x[8];
  matr[7] = (2.0*x[10] - f)*sqrt3;
  matr[8] = (3.0*x[11] - x[10] - f)*sqrt6;

  /* Calculate det(M). */
  loc1 = matr[4]*matr[8] - matr[5]*matr[7];
  loc2 = matr[5]*matr[6] - matr[3]*matr[8];
  loc3 = matr[3]*matr[7] - matr[4]*matr[6];
  g = matr[0]*loc1 + matr[1]*loc2 + matr[2]*loc3;
  if (g <= epsilon) { *obj = g; return 1; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] + 
      matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8];
 
  /* Calculate objective function. */
  /* loc4 = a * pow(g, b); */
  loc4 = a * rcbrt(g*g);
  *obj = f * loc4;

  /* Calculate the derivative of the objective function. */
  f = 2.0*loc4;
  g = b*(*obj)/g; 

  adj_m[0] = matr[0]*f + loc1*g;
  adj_m[1] = matr[1]*f + loc2*g;
  adj_m[2] = matr[2]*f + loc3*g;

  loc1 = matr[0]*g;
  loc2 = matr[1]*g;
  loc3 = matr[2]*g;

  adj_m[3] = matr[3]*f + loc3*matr[7] - loc2*matr[8];
  adj_m[4] = matr[4]*f + loc1*matr[8] - loc3*matr[6];
  adj_m[5] = matr[5]*f + loc2*matr[6] - loc1*matr[7];

  adj_m[6] = matr[6]*f + loc2*matr[5] - loc3*matr[4];
  adj_m[7] = matr[7]*f + loc3*matr[3] - loc1*matr[5];
  adj_m[8] = matr[8]*f + loc1*matr[4] - loc2*matr[3];

  loc1 = sqrt3*adj_m[1];
  loc2 = sqrt6*adj_m[2];
  loc3 = loc1 + loc2;
  g_obj[0] = -adj_m[0] - loc3;
  g_obj[1] = adj_m[0] - loc3;
  g_obj[2] = 2.0*loc1 - loc2;
  g_obj[3] = 3.0*loc2;

  loc1 = sqrt3*adj_m[4];
  loc2 = sqrt6*adj_m[5];
  loc3 = loc1 + loc2;
  g_obj[4] = -adj_m[3] - loc3;
  g_obj[5] = adj_m[3] - loc3;
  g_obj[6] = 2.0*loc1 - loc2;
  g_obj[7] = 3.0*loc2;

  loc1 = sqrt3*adj_m[7];
  loc2 = sqrt6*adj_m[8];
  loc3 = loc1 + loc2;
  g_obj[8] = -adj_m[6] - loc3;
  g_obj[9] = adj_m[6] - loc3;
  g_obj[10] = 2.0*loc1 - loc2;
  g_obj[11] = 3.0*loc2;
  return 0;
}

/*****************************************************************************/
/* The Hessian calculation is done by blocks.  Only the upper triangular     */
/* blocks are stored.  The results in the data is in the following order:    */
/*    [d1 b1 b2 d2 b3 d3 ]                                                   */
/* The matrices on the diagonal (d1-d3) each contain 10 elements, while the  */
/* off-diagonal elements (b1-b3) each contain 16 elements.                   */
/*                                                                           */
/* The code requires 598 flops.  The gradient evaluation needs 130 flops     */
/* and the function requires 61 flops.                                       */
/*****************************************************************************/
/* The form of the function, gradient, and Hessian is the following:         */
/*   o(x) = a * f(A(x)) * pow(g(A(x)), b)                                    */
/* where A(x) is the matrix generated from:                                  */
/*           [x1-x0 x2-x0 x3-x0]                                             */
/*    A(x) = [y1-y0 y2-y0 y3-y0] * inv(W)                                    */
/*           [z1-z0 z2-z0 z3-z0]                                             */
/* and f() is the squared Frobenius norm of A(x), and g() is the determinant */
/* of A(x).                                                                  */
/*                                                                           */
/* The gradient is calculated as follows:                                    */
/*   alpha := a*pow(g(A(x)),b)                                               */
/*   beta  := a*b*f(A(x))*pow(g(A(x)),b-1)                                   */
/*                                                                           */
/*   do/dx = (alpha * (df/dA) + beta * (dg/dA)) (dA/dx)                      */
/*                                                                           */
/*   Note: this is the optimal ordering for the gradient vector.             */
/*   Distributing (dA/dx) would result in two matrix vector products as      */
/*   opposed to the single matrix vector product in the above formulation.   */
/*                                                                           */
/*   (df/dA)_i = 2*A_i                                                       */
/*   (dg/dA)_i = A_j*A_k - A_l*A_m for some {j,k,l,m}                        */
/*                                                                           */
/*   d^2o/dx^2 = (dA/dx)' * ((d alpha/dA) * (df/dA) +                        */
/*                           (d  beta/dA) * (dg/dA)                          */
/*			            alpha * (d^2f/dA^2)                      */
/*                                   beta * (d^2g/dA^2)) * (dA/dx)           */
/*                                                                           */
/*   Note: since A(x) is a linear function, there are no terms involving     */
/*   d^2A/dx^2 since this matrix is zero.                                    */
/*                                                                           */
/*   gamma := a*b*pow(g(A(x)),b-1)                                           */
/*   delta := a*b*(b-1)*f(A(x))*pow(g(A(x)),b-2)                             */
/*                                                                           */
/*   d^2o/dx^2 = (dA/dx)' * (gamma*((dg/dA)'*(df/dA) + (df/dA)'*(dg/dA)) +   */
/*                           delta* (dg/dA)'*(dg/dA) +                       */
/*                           alpha*(d^2f/dA^2) +                             */
/*                            beta*(d^2g/dA^2)) * (dA/dx)                    */
/*                                                                           */
/*   Note: (df/dA) and (dg/dA) are row vectors and we only calculate the     */
/*   upper triangular part of the inner matrix.                              */
/*                                                                           */
/*   For regular tetrahedral elements, we have the following:                */
/*                                                                           */
/*           [-1         1	  0         0         ]                      */
/*       M = [-sqrt(3)  -sqrt(3)  2*sqrt(3) 0         ]                      */
/*           [-sqrt(6)  -sqrt(6)  -sqrt(6)  3*sqrt(6) ]                      */
/*                                                                           */
/*           [M 0 0]                                                         */
/*   dA/dx = [0 M 0]                                                         */
/*           [0 0 M]                                                         */
/*                                                                           */
/*   I belive the above is close to optimal for the calculation of the       */
/*   Hessian.  Distributing the (dA/dx) results in larger vector which are   */
/*   detrimental when forming the outer product.  The way the method is      */
/*   written, we only calculate a 9x9 symmetric matrix in the outer product. */
/*                                                                           */
/*   In two dimensions, the inner matrix computed has a nice structure and   */
/*   we can eliminate some of the computation in the inner product.  This    */
/*   does not appear to be the case in more than two dimensions.             */
/*****************************************************************************/

int h_fcn(double *obj, double g_obj[12], double h_obj[78], const double x[12])
{
  static double matr[9], f;
  static double adj_m[9], g;
  static double dg[9], loc0, loc1, loc2, loc3, loc4;
  static double A[12], J_A[6], J_B[9], J_C[9];

  /* Calculate M = A*inv(W). */
  f       = x[1] + x[0];
  matr[0] = x[1] - x[0];
  matr[1] = (2.0*x[2] - f)*sqrt3;
  matr[2] = (3.0*x[3] - x[2] - f)*sqrt6;

  f       = x[5] + x[4];
  matr[3] = x[5] - x[4];
  matr[4] = (2.0*x[6] - f)*sqrt3;
  matr[5] = (3.0*x[7] - x[6] - f)*sqrt6;

  f       = x[9] + x[8];
  matr[6] = x[9] - x[8];
  matr[7] = (2.0*x[10] - f)*sqrt3;
  matr[8] = (3.0*x[11] - x[10] - f)*sqrt6;

  /* Calculate det(M). */
  dg[0] = matr[4]*matr[8] - matr[5]*matr[7];
  dg[1] = matr[5]*matr[6] - matr[3]*matr[8];
  dg[2] = matr[3]*matr[7] - matr[4]*matr[6];
  g = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];
  if (g <= epsilon) { *obj = g; return 1; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] + 
      matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8];

  loc4 = g;

  /* Calculate objective function. */
  /* loc1 = a * pow(g, b); */
  loc1 = a * rcbrt(g*g);
  *obj = f * loc1;

  /* Calculate the derivative of the objective function. */
  f = 2.0*loc1;
  g = b*(*obj)/g; 

  dg[3] = matr[2]*matr[7] - matr[1]*matr[8];
  dg[4] = matr[0]*matr[8] - matr[2]*matr[6];
  dg[5] = matr[1]*matr[6] - matr[0]*matr[7];
  dg[6] = matr[1]*matr[5] - matr[2]*matr[4];
  dg[7] = matr[2]*matr[3] - matr[0]*matr[5];
  dg[8] = matr[0]*matr[4] - matr[1]*matr[3];

  adj_m[0] = matr[0]*f + dg[0]*g;
  adj_m[1] = matr[1]*f + dg[1]*g;
  adj_m[2] = matr[2]*f + dg[2]*g;
  adj_m[3] = matr[3]*f + dg[3]*g;
  adj_m[4] = matr[4]*f + dg[4]*g;
  adj_m[5] = matr[5]*f + dg[5]*g;
  adj_m[6] = matr[6]*f + dg[6]*g;
  adj_m[7] = matr[7]*f + dg[7]*g;
  adj_m[8] = matr[8]*f + dg[8]*g;

  loc1 = sqrt3*adj_m[1];
  loc2 = sqrt6*adj_m[2];
  loc3 = loc1 + loc2;
  g_obj[0] = -adj_m[0] - loc3;
  g_obj[1] = adj_m[0] - loc3;
  g_obj[2] = 2.0*loc1 - loc2;
  g_obj[3] = 3.0*loc2;

  loc1 = sqrt3*adj_m[4];
  loc2 = sqrt6*adj_m[5];
  loc3 = loc1 + loc2;
  g_obj[4] = -adj_m[3] - loc3;
  g_obj[5] = adj_m[3] - loc3;
  g_obj[6] = 2.0*loc1 - loc2;
  g_obj[7] = 3.0*loc2;

  loc1 = sqrt3*adj_m[7];
  loc2 = sqrt6*adj_m[8];
  loc3 = loc1 + loc2;
  g_obj[8] = -adj_m[6] - loc3;
  g_obj[9] = adj_m[6] - loc3;
  g_obj[10] = 2.0*loc1 - loc2;
  g_obj[11] = 3.0*loc2;

  loc0 = g;
  loc1 = f;
  f = f*b/loc4;
  g = g*bm1/loc4;

  /* First block of rows */
  loc2 = matr[0]*f;
  loc3 = dg[0]*f;
  loc4 = dg[0]*g + loc2;

  J_A[0] = loc1 + dg[0]*(loc2 + loc4);
  J_A[1] = loc3*matr[1] + loc4*dg[1];
  J_A[2] = loc3*matr[2] + loc4*dg[2];
  J_B[0] = loc3*matr[3] + loc4*dg[3];
  J_B[1] = loc3*matr[4] + loc4*dg[4];
  J_B[2] = loc3*matr[5] + loc4*dg[5];
  J_C[0] = loc3*matr[6] + loc4*dg[6];
  J_C[1] = loc3*matr[7] + loc4*dg[7];
  J_C[2] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[1]*f;
  loc3 = dg[1]*f;
  loc4 = dg[1]*g + loc2;

  J_A[3] = loc1 + dg[1]*(loc2 + loc4);
  J_A[4] = loc3*matr[2] + loc4*dg[2];
  J_B[3] = loc3*matr[3] + loc4*dg[3];
  J_B[4] = loc3*matr[4] + loc4*dg[4];
  J_B[5] = loc3*matr[5] + loc4*dg[5];
  J_C[3] = loc3*matr[6] + loc4*dg[6];
  J_C[4] = loc3*matr[7] + loc4*dg[7];
  J_C[5] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[2]*f;
  loc3 = dg[2]*f;
  loc4 = dg[2]*g + loc2;

  J_A[5] = loc1 + dg[2]*(loc2 + loc4);
  J_B[6] = loc3*matr[3] + loc4*dg[3];
  J_B[7] = loc3*matr[4] + loc4*dg[4];
  J_B[8] = loc3*matr[5] + loc4*dg[5];
  J_C[6] = loc3*matr[6] + loc4*dg[6];
  J_C[7] = loc3*matr[7] + loc4*dg[7];
  J_C[8] = loc3*matr[8] + loc4*dg[8];

  /* First diagonal block */
  loc2 = sqrt3*J_A[1];
  loc3 = sqrt6*J_A[2];
  loc4 = loc2 + loc3;

  A[0] = -J_A[0] - loc4;
  A[1] =  J_A[0] - loc4;

  loc2 = sqrt3*J_A[3];
  loc3 = sqrt6*J_A[4];
  loc4 = loc2 + loc3;

  A[4] = -J_A[1] - loc4;
  A[5] =  J_A[1] - loc4;
  A[6] = 2.0*loc2 - loc3;

  loc2 = sqrt3*J_A[4];
  loc3 = sqrt6*J_A[5];
  loc4 = loc2 + loc3;

  A[8] = -J_A[2] - loc4;
  A[9] =  J_A[2] - loc4;
  A[10] = 2.0*loc2 - loc3;
  A[11] = 3.0*loc3;

  loc2 = sqrt3*A[4];
  loc3 = sqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[0] = -A[0] - loc4;
  h_obj[1] =  A[0] - loc4;
  h_obj[2] = 2.0*loc2 - loc3;
  h_obj[3] = 3.0*loc3;

  loc2 = sqrt3*A[5];
  loc3 = sqrt6*A[9];

  h_obj[4] = A[1] - loc2 - loc3;
  h_obj[5] = 2.0*loc2 - loc3;
  h_obj[6] = 3.0*loc3;

  loc3 = sqrt6*A[10];
  h_obj[7] = tsqrt3*A[6] - loc3;
  h_obj[8] = 3.0*loc3;

  h_obj[9] = tsqrt6*A[11];

  /* First off-diagonal block */
  loc2 = matr[8]*loc0;
  J_B[1] += loc2;
  J_B[3] -= loc2;

  loc2 = matr[7]*loc0;
  J_B[2] -= loc2;
  J_B[6] += loc2;

  loc2 = matr[6]*loc0;
  J_B[5] += loc2;
  J_B[7] -= loc2;

  loc2 = sqrt3*J_B[3];
  loc3 = sqrt6*J_B[6];
  loc4 = loc2 + loc3;

  A[0] = -J_B[0] - loc4;
  A[1] =  J_B[0] - loc4;
  A[2] = 2.0*loc2 - loc3;
  A[3] = 3.0*loc3;

  loc2 = sqrt3*J_B[4];
  loc3 = sqrt6*J_B[7];
  loc4 = loc2 + loc3;

  A[4] = -J_B[1] - loc4;
  A[5] =  J_B[1] - loc4;
  A[6] = 2.0*loc2 - loc3;
  A[7] = 3.0*loc3;

  loc2 = sqrt3*J_B[5];
  loc3 = sqrt6*J_B[8];
  loc4 = loc2 + loc3;

  A[8] = -J_B[2] - loc4;
  A[9] =  J_B[2] - loc4;
  A[10] = 2.0*loc2 - loc3;
  A[11] = 3.0*loc3;

  loc2 = sqrt3*A[4];
  loc3 = sqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[10] = -A[0] - loc4;
  h_obj[11] =  A[0] - loc4;
  h_obj[12] = 2.0*loc2 - loc3;
  h_obj[13] = 3.0*loc3;

  loc2 = sqrt3*A[5];
  loc3 = sqrt6*A[9];
  loc4 = loc2 + loc3;

  h_obj[14] = -A[1] - loc4;
  h_obj[15] =  A[1] - loc4;
  h_obj[16] = 2.0*loc2 - loc3;
  h_obj[17] = 3.0*loc3;

  loc2 = sqrt3*A[6];
  loc3 = sqrt6*A[10];
  loc4 = loc2 + loc3;

  h_obj[18] = -A[2] - loc4;
  h_obj[19] =  A[2] - loc4;
  h_obj[20] = 2.0*loc2 - loc3;
  h_obj[21] = 3.0*loc3;

  loc2 = sqrt3*A[7];
  loc3 = sqrt6*A[11];
  loc4 = loc2 + loc3;

  h_obj[22] = -A[3] - loc4;
  h_obj[23] =  A[3] - loc4;
  h_obj[24] = 2.0*loc2 - loc3;
  h_obj[25] = 3.0*loc3;

  /* Second off-diagonal block */
  loc2 = matr[5]*loc0;
  J_C[1] -= loc2;
  J_C[3] += loc2;

  loc2 = matr[4]*loc0;
  J_C[2] += loc2;
  J_C[6] -= loc2;

  loc2 = matr[3]*loc0;
  J_C[5] -= loc2;
  J_C[7] += loc2;

  loc2 = sqrt3*J_C[3];
  loc3 = sqrt6*J_C[6];
  loc4 = loc2 + loc3;

  A[0] = -J_C[0] - loc4;
  A[1] =  J_C[0] - loc4;
  A[2] = 2.0*loc2 - loc3;
  A[3] = 3.0*loc3;

  loc2 = sqrt3*J_C[4];
  loc3 = sqrt6*J_C[7];
  loc4 = loc2 + loc3;

  A[4] = -J_C[1] - loc4;
  A[5] =  J_C[1] - loc4;
  A[6] = 2.0*loc2 - loc3;
  A[7] = 3.0*loc3;

  loc2 = sqrt3*J_C[5];
  loc3 = sqrt6*J_C[8];
  loc4 = loc2 + loc3;

  A[8] = -J_C[2] - loc4;
  A[9] =  J_C[2] - loc4;
  A[10] = 2.0*loc2 - loc3;
  A[11] = 3.0*loc3;

  loc2 = sqrt3*A[4];
  loc3 = sqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[26] = -A[0] - loc4;
  h_obj[27] =  A[0] - loc4;
  h_obj[28] = 2.0*loc2 - loc3;
  h_obj[29] = 3.0*loc3;

  loc2 = sqrt3*A[5];
  loc3 = sqrt6*A[9];
  loc4 = loc2 + loc3;

  h_obj[30] = -A[1] - loc4;
  h_obj[31] =  A[1] - loc4;
  h_obj[32] = 2.0*loc2 - loc3;
  h_obj[33] = 3.0*loc3;

  loc2 = sqrt3*A[6];
  loc3 = sqrt6*A[10];
  loc4 = loc2 + loc3;

  h_obj[34] = -A[2] - loc4;
  h_obj[35] =  A[2] - loc4;
  h_obj[36] = 2.0*loc2 - loc3;
  h_obj[37] = 3.0*loc3;

  loc2 = sqrt3*A[7];
  loc3 = sqrt6*A[11];
  loc4 = loc2 + loc3;

  h_obj[38] = -A[3] - loc4;
  h_obj[39] =  A[3] - loc4;
  h_obj[40] = 2.0*loc2 - loc3;
  h_obj[41] = 3.0*loc3;

  /* Second block of rows */
  loc2 = matr[3]*f;
  loc3 = dg[3]*f;
  loc4 = dg[3]*g + loc2;

  J_A[0] = loc1 + dg[3]*(loc2 + loc4);
  J_A[1] = loc3*matr[4] + loc4*dg[4];
  J_A[2] = loc3*matr[5] + loc4*dg[5];
  J_B[0] = loc3*matr[6] + loc4*dg[6];
  J_B[1] = loc3*matr[7] + loc4*dg[7];
  J_B[2] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[4]*f;
  loc3 = dg[4]*f;
  loc4 = dg[4]*g + loc2;

  J_A[3] = loc1 + dg[4]*(loc2 + loc4);
  J_A[4] = loc3*matr[5] + loc4*dg[5];
  J_B[3] = loc3*matr[6] + loc4*dg[6];
  J_B[4] = loc3*matr[7] + loc4*dg[7];
  J_B[5] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[5]*f;
  loc3 = dg[5]*f;
  loc4 = dg[5]*g + loc2;

  J_A[5] = loc1 + dg[5]*(loc2 + loc4);
  J_B[6] = loc3*matr[6] + loc4*dg[6];
  J_B[7] = loc3*matr[7] + loc4*dg[7];
  J_B[8] = loc3*matr[8] + loc4*dg[8];

  /* Second diagonal block */
  loc2 = sqrt3*J_A[1];
  loc3 = sqrt6*J_A[2];
  loc4 = loc2 + loc3;

  A[0] = -J_A[0] - loc4;
  A[1] =  J_A[0] - loc4;

  loc2 = sqrt3*J_A[3];
  loc3 = sqrt6*J_A[4];
  loc4 = loc2 + loc3;

  A[4] = -J_A[1] - loc4;
  A[5] =  J_A[1] - loc4;
  A[6] = 2.0*loc2 - loc3;

  loc2 = sqrt3*J_A[4];
  loc3 = sqrt6*J_A[5];
  loc4 = loc2 + loc3;

  A[8] = -J_A[2] - loc4;
  A[9] =  J_A[2] - loc4;
  A[10] = 2.0*loc2 - loc3;
  A[11] = 3.0*loc3;

  loc2 = sqrt3*A[4];
  loc3 = sqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[42] = -A[0] - loc4;
  h_obj[43] =  A[0] - loc4;
  h_obj[44] = 2.0*loc2 - loc3;
  h_obj[45] = 3.0*loc3;

  loc2 = sqrt3*A[5];
  loc3 = sqrt6*A[9];

  h_obj[46] = A[1] - loc2 - loc3;
  h_obj[47] = 2.0*loc2 - loc3;
  h_obj[48] = 3.0*loc3;

  loc3 = sqrt6*A[10];
  h_obj[49] = tsqrt3*A[6] - loc3;
  h_obj[50] = 3.0*loc3;

  h_obj[51] = tsqrt6*A[11];

  /* Third off-diagonal block */
  loc2 = matr[2]*loc0;
  J_B[1] += loc2;
  J_B[3] -= loc2;

  loc2 = matr[1]*loc0;
  J_B[2] -= loc2;
  J_B[6] += loc2;

  loc2 = matr[0]*loc0;
  J_B[5] += loc2;
  J_B[7] -= loc2;

  loc2 = sqrt3*J_B[3];
  loc3 = sqrt6*J_B[6];
  loc4 = loc2 + loc3;

  A[0] = -J_B[0] - loc4;
  A[1] =  J_B[0] - loc4;
  A[2] = 2.0*loc2 - loc3;
  A[3] = 3.0*loc3;

  loc2 = sqrt3*J_B[4];
  loc3 = sqrt6*J_B[7];
  loc4 = loc2 + loc3;

  A[4] = -J_B[1] - loc4;
  A[5] =  J_B[1] - loc4;
  A[6] = 2.0*loc2 - loc3;
  A[7] = 3.0*loc3;

  loc2 = sqrt3*J_B[5];
  loc3 = sqrt6*J_B[8];
  loc4 = loc2 + loc3;

  A[8] = -J_B[2] - loc4;
  A[9] =  J_B[2] - loc4;
  A[10] = 2.0*loc2 - loc3;
  A[11] = 3.0*loc3;

  loc2 = sqrt3*A[4];
  loc3 = sqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[52] = -A[0] - loc4;
  h_obj[53] =  A[0] - loc4;
  h_obj[54] = 2.0*loc2 - loc3;
  h_obj[55] = 3.0*loc3;

  loc2 = sqrt3*A[5];
  loc3 = sqrt6*A[9];
  loc4 = loc2 + loc3;

  h_obj[56] = -A[1] - loc4;
  h_obj[57] =  A[1] - loc4;
  h_obj[58] = 2.0*loc2 - loc3;
  h_obj[59] = 3.0*loc3;

  loc2 = sqrt3*A[6];
  loc3 = sqrt6*A[10];
  loc4 = loc2 + loc3;

  h_obj[60] = -A[2] - loc4;
  h_obj[61] =  A[2] - loc4;
  h_obj[62] = 2.0*loc2 - loc3;
  h_obj[63] = 3.0*loc3;

  loc2 = sqrt3*A[7];
  loc3 = sqrt6*A[11];
  loc4 = loc2 + loc3;

  h_obj[64] = -A[3] - loc4;
  h_obj[65] =  A[3] - loc4;
  h_obj[66] = 2.0*loc2 - loc3;
  h_obj[67] = 3.0*loc3;

  /* Third block of rows */
  loc2 = matr[6]*f;
  loc3 = dg[6]*f;
  loc4 = dg[6]*g + loc2;

  J_A[0] = loc1 + dg[6]*(loc2 + loc4);
  J_A[1] = loc3*matr[7] + loc4*dg[7];
  J_A[2] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[7]*f;
  loc3 = dg[7]*f;
  loc4 = dg[7]*g + loc2;

  J_A[3] = loc1 + dg[7]*(loc2 + loc4);
  J_A[4] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[8]*f;
  loc4 = dg[8]*g + loc2;

  J_A[5] = loc1 + dg[8]*(loc2 + loc4);

  /* Third diagonal block */
  loc2 = sqrt3*J_A[1];
  loc3 = sqrt6*J_A[2];
  loc4 = loc2 + loc3;

  A[0] = -J_A[0] - loc4;
  A[1] =  J_A[0] - loc4;

  loc2 = sqrt3*J_A[3];
  loc3 = sqrt6*J_A[4];
  loc4 = loc2 + loc3;

  A[4] = -J_A[1] - loc4;
  A[5] =  J_A[1] - loc4;
  A[6] = 2.0*loc2 - loc3;

  loc2 = sqrt3*J_A[4];
  loc3 = sqrt6*J_A[5];
  loc4 = loc2 + loc3;

  A[8] = -J_A[2] - loc4;
  A[9] =  J_A[2] - loc4;
  A[10] = 2.0*loc2 - loc3;
  A[11] = 3.0*loc3;

  loc2 = sqrt3*A[4];
  loc3 = sqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[68] = -A[0] - loc4;
  h_obj[69] =  A[0] - loc4;
  h_obj[70] = 2.0*loc2 - loc3;
  h_obj[71] = 3.0*loc3;

  loc2 = sqrt3*A[5];
  loc3 = sqrt6*A[9];

  h_obj[72] = A[1] - loc2 - loc3;
  h_obj[73] = 2.0*loc2 - loc3;
  h_obj[74] = 3.0*loc3;

  loc3 = sqrt6*A[10];
  h_obj[75] = tsqrt3*A[6] - loc3;
  h_obj[76] = 3.0*loc3;

  h_obj[77] = tsqrt6*A[11];
  return 0;
}

void h_only(double h_obj[78], const double x[12])
{
  static double matr[9], obj;
  static double dg[9], f, g, loc0, loc1, loc2, loc3, loc4;
  static double A[12], J_A[6], J_B[9], J_C[9];

  /* Calculate M = A*inv(W). */
  f       = x[1] + x[0];
  matr[0] = x[1] - x[0];
  matr[1] = (2.0*x[2] - f)*sqrt3;
  matr[2] = (3.0*x[3] - x[2] - f)*sqrt6;

  f       = x[5] + x[4];
  matr[3] = x[5] - x[4];
  matr[4] = (2.0*x[6] - f)*sqrt3;
  matr[5] = (3.0*x[7] - x[6] - f)*sqrt6;

  f       = x[9] + x[8];
  matr[6] = x[9] - x[8];
  matr[7] = (2.0*x[10] - f)*sqrt3;
  matr[8] = (3.0*x[11] - x[10] - f)*sqrt6;

  /* Calculate det(M). */
  dg[0] = matr[4]*matr[8] - matr[5]*matr[7];
  dg[1] = matr[5]*matr[6] - matr[3]*matr[8];
  dg[2] = matr[3]*matr[7] - matr[4]*matr[6];
  g = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] + 
      matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8];

  loc4 = g;

  /* Calculate objective function. */
  /* loc1 = a * pow(g, b); */
  loc1 = a * rcbrt(g*g);
  obj = f * loc1;

  /* Calculate the derivative of the objective function. */
  f = 2.0*loc1;
  g = b*obj/g; 

  dg[3] = matr[2]*matr[7] - matr[1]*matr[8];
  dg[4] = matr[0]*matr[8] - matr[2]*matr[6];
  dg[5] = matr[1]*matr[6] - matr[0]*matr[7];
  dg[6] = matr[1]*matr[5] - matr[2]*matr[4];
  dg[7] = matr[2]*matr[3] - matr[0]*matr[5];
  dg[8] = matr[0]*matr[4] - matr[1]*matr[3];

  loc0 = g;
  loc1 = f;
  f = f*b/loc4;
  g = g*bm1/loc4;

  /* First block of rows */
  loc2 = matr[0]*f;
  loc3 = dg[0]*f;
  loc4 = dg[0]*g + loc2;

  J_A[0] = loc1 + dg[0]*(loc2 + loc4);
  J_A[1] = loc3*matr[1] + loc4*dg[1];
  J_A[2] = loc3*matr[2] + loc4*dg[2];
  J_B[0] = loc3*matr[3] + loc4*dg[3];
  J_B[1] = loc3*matr[4] + loc4*dg[4];
  J_B[2] = loc3*matr[5] + loc4*dg[5];
  J_C[0] = loc3*matr[6] + loc4*dg[6];
  J_C[1] = loc3*matr[7] + loc4*dg[7];
  J_C[2] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[1]*f;
  loc3 = dg[1]*f;
  loc4 = dg[1]*g + loc2;

  J_A[3] = loc1 + dg[1]*(loc2 + loc4);
  J_A[4] = loc3*matr[2] + loc4*dg[2];
  J_B[3] = loc3*matr[3] + loc4*dg[3];
  J_B[4] = loc3*matr[4] + loc4*dg[4];
  J_B[5] = loc3*matr[5] + loc4*dg[5];
  J_C[3] = loc3*matr[6] + loc4*dg[6];
  J_C[4] = loc3*matr[7] + loc4*dg[7];
  J_C[5] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[2]*f;
  loc3 = dg[2]*f;
  loc4 = dg[2]*g + loc2;

  J_A[5] = loc1 + dg[2]*(loc2 + loc4);
  J_B[6] = loc3*matr[3] + loc4*dg[3];
  J_B[7] = loc3*matr[4] + loc4*dg[4];
  J_B[8] = loc3*matr[5] + loc4*dg[5];
  J_C[6] = loc3*matr[6] + loc4*dg[6];
  J_C[7] = loc3*matr[7] + loc4*dg[7];
  J_C[8] = loc3*matr[8] + loc4*dg[8];

  /* First diagonal block */
  loc2 = sqrt3*J_A[1];
  loc3 = sqrt6*J_A[2];
  loc4 = loc2 + loc3;

  A[0] = -J_A[0] - loc4;
  A[1] =  J_A[0] - loc4;

  loc2 = sqrt3*J_A[3];
  loc3 = sqrt6*J_A[4];
  loc4 = loc2 + loc3;

  A[4] = -J_A[1] - loc4;
  A[5] =  J_A[1] - loc4;
  A[6] = 2.0*loc2 - loc3;

  loc2 = sqrt3*J_A[4];
  loc3 = sqrt6*J_A[5];
  loc4 = loc2 + loc3;

  A[8] = -J_A[2] - loc4;
  A[9] =  J_A[2] - loc4;
  A[10] = 2.0*loc2 - loc3;
  A[11] = 3.0*loc3;

  loc2 = sqrt3*A[4];
  loc3 = sqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[0] = -A[0] - loc4;
  h_obj[1] =  A[0] - loc4;
  h_obj[2] = 2.0*loc2 - loc3;
  h_obj[3] = 3.0*loc3;

  loc2 = sqrt3*A[5];
  loc3 = sqrt6*A[9];

  h_obj[4] = A[1] - loc2 - loc3;
  h_obj[5] = 2.0*loc2 - loc3;
  h_obj[6] = 3.0*loc3;

  loc3 = sqrt6*A[10];
  h_obj[7] = tsqrt3*A[6] - loc3;
  h_obj[8] = 3.0*loc3;

  h_obj[9] = tsqrt6*A[11];

  /* First off-diagonal block */
  loc2 = matr[8]*loc0;
  J_B[1] += loc2;
  J_B[3] -= loc2;

  loc2 = matr[7]*loc0;
  J_B[2] -= loc2;
  J_B[6] += loc2;

  loc2 = matr[6]*loc0;
  J_B[5] += loc2;
  J_B[7] -= loc2;

  loc2 = sqrt3*J_B[3];
  loc3 = sqrt6*J_B[6];
  loc4 = loc2 + loc3;

  A[0] = -J_B[0] - loc4;
  A[1] =  J_B[0] - loc4;
  A[2] = 2.0*loc2 - loc3;
  A[3] = 3.0*loc3;

  loc2 = sqrt3*J_B[4];
  loc3 = sqrt6*J_B[7];
  loc4 = loc2 + loc3;

  A[4] = -J_B[1] - loc4;
  A[5] =  J_B[1] - loc4;
  A[6] = 2.0*loc2 - loc3;
  A[7] = 3.0*loc3;

  loc2 = sqrt3*J_B[5];
  loc3 = sqrt6*J_B[8];
  loc4 = loc2 + loc3;

  A[8] = -J_B[2] - loc4;
  A[9] =  J_B[2] - loc4;
  A[10] = 2.0*loc2 - loc3;
  A[11] = 3.0*loc3;

  loc2 = sqrt3*A[4];
  loc3 = sqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[10] = -A[0] - loc4;
  h_obj[11] =  A[0] - loc4;
  h_obj[12] = 2.0*loc2 - loc3;
  h_obj[13] = 3.0*loc3;

  loc2 = sqrt3*A[5];
  loc3 = sqrt6*A[9];
  loc4 = loc2 + loc3;

  h_obj[14] = -A[1] - loc4;
  h_obj[15] =  A[1] - loc4;
  h_obj[16] = 2.0*loc2 - loc3;
  h_obj[17] = 3.0*loc3;

  loc2 = sqrt3*A[6];
  loc3 = sqrt6*A[10];
  loc4 = loc2 + loc3;

  h_obj[18] = -A[2] - loc4;
  h_obj[19] =  A[2] - loc4;
  h_obj[20] = 2.0*loc2 - loc3;
  h_obj[21] = 3.0*loc3;

  loc2 = sqrt3*A[7];
  loc3 = sqrt6*A[11];
  loc4 = loc2 + loc3;

  h_obj[22] = -A[3] - loc4;
  h_obj[23] =  A[3] - loc4;
  h_obj[24] = 2.0*loc2 - loc3;
  h_obj[25] = 3.0*loc3;

  /* Second off-diagonal block */
  loc2 = matr[5]*loc0;
  J_C[1] -= loc2;
  J_C[3] += loc2;

  loc2 = matr[4]*loc0;
  J_C[2] += loc2;
  J_C[6] -= loc2;

  loc2 = matr[3]*loc0;
  J_C[5] -= loc2;
  J_C[7] += loc2;

  loc2 = sqrt3*J_C[3];
  loc3 = sqrt6*J_C[6];
  loc4 = loc2 + loc3;

  A[0] = -J_C[0] - loc4;
  A[1] =  J_C[0] - loc4;
  A[2] = 2.0*loc2 - loc3;
  A[3] = 3.0*loc3;

  loc2 = sqrt3*J_C[4];
  loc3 = sqrt6*J_C[7];
  loc4 = loc2 + loc3;

  A[4] = -J_C[1] - loc4;
  A[5] =  J_C[1] - loc4;
  A[6] = 2.0*loc2 - loc3;
  A[7] = 3.0*loc3;

  loc2 = sqrt3*J_C[5];
  loc3 = sqrt6*J_C[8];
  loc4 = loc2 + loc3;

  A[8] = -J_C[2] - loc4;
  A[9] =  J_C[2] - loc4;
  A[10] = 2.0*loc2 - loc3;
  A[11] = 3.0*loc3;

  loc2 = sqrt3*A[4];
  loc3 = sqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[26] = -A[0] - loc4;
  h_obj[27] =  A[0] - loc4;
  h_obj[28] = 2.0*loc2 - loc3;
  h_obj[29] = 3.0*loc3;

  loc2 = sqrt3*A[5];
  loc3 = sqrt6*A[9];
  loc4 = loc2 + loc3;

  h_obj[30] = -A[1] - loc4;
  h_obj[31] =  A[1] - loc4;
  h_obj[32] = 2.0*loc2 - loc3;
  h_obj[33] = 3.0*loc3;

  loc2 = sqrt3*A[6];
  loc3 = sqrt6*A[10];
  loc4 = loc2 + loc3;

  h_obj[34] = -A[2] - loc4;
  h_obj[35] =  A[2] - loc4;
  h_obj[36] = 2.0*loc2 - loc3;
  h_obj[37] = 3.0*loc3;

  loc2 = sqrt3*A[7];
  loc3 = sqrt6*A[11];
  loc4 = loc2 + loc3;

  h_obj[38] = -A[3] - loc4;
  h_obj[39] =  A[3] - loc4;
  h_obj[40] = 2.0*loc2 - loc3;
  h_obj[41] = 3.0*loc3;

  /* Second block of rows */
  loc2 = matr[3]*f;
  loc3 = dg[3]*f;
  loc4 = dg[3]*g + loc2;

  J_A[0] = loc1 + dg[3]*(loc2 + loc4);
  J_A[1] = loc3*matr[4] + loc4*dg[4];
  J_A[2] = loc3*matr[5] + loc4*dg[5];
  J_B[0] = loc3*matr[6] + loc4*dg[6];
  J_B[1] = loc3*matr[7] + loc4*dg[7];
  J_B[2] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[4]*f;
  loc3 = dg[4]*f;
  loc4 = dg[4]*g + loc2;

  J_A[3] = loc1 + dg[4]*(loc2 + loc4);
  J_A[4] = loc3*matr[5] + loc4*dg[5];
  J_B[3] = loc3*matr[6] + loc4*dg[6];
  J_B[4] = loc3*matr[7] + loc4*dg[7];
  J_B[5] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[5]*f;
  loc3 = dg[5]*f;
  loc4 = dg[5]*g + loc2;

  J_A[5] = loc1 + dg[5]*(loc2 + loc4);
  J_B[6] = loc3*matr[6] + loc4*dg[6];
  J_B[7] = loc3*matr[7] + loc4*dg[7];
  J_B[8] = loc3*matr[8] + loc4*dg[8];

  /* Second diagonal block */
  loc2 = sqrt3*J_A[1];
  loc3 = sqrt6*J_A[2];
  loc4 = loc2 + loc3;

  A[0] = -J_A[0] - loc4;
  A[1] =  J_A[0] - loc4;

  loc2 = sqrt3*J_A[3];
  loc3 = sqrt6*J_A[4];
  loc4 = loc2 + loc3;

  A[4] = -J_A[1] - loc4;
  A[5] =  J_A[1] - loc4;
  A[6] = 2.0*loc2 - loc3;

  loc2 = sqrt3*J_A[4];
  loc3 = sqrt6*J_A[5];
  loc4 = loc2 + loc3;

  A[8] = -J_A[2] - loc4;
  A[9] =  J_A[2] - loc4;
  A[10] = 2.0*loc2 - loc3;
  A[11] = 3.0*loc3;

  loc2 = sqrt3*A[4];
  loc3 = sqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[42] = -A[0] - loc4;
  h_obj[43] =  A[0] - loc4;
  h_obj[44] = 2.0*loc2 - loc3;
  h_obj[45] = 3.0*loc3;

  loc2 = sqrt3*A[5];
  loc3 = sqrt6*A[9];

  h_obj[46] = A[1] - loc2 - loc3;
  h_obj[47] = 2.0*loc2 - loc3;
  h_obj[48] = 3.0*loc3;

  loc3 = sqrt6*A[10];
  h_obj[49] = tsqrt3*A[6] - loc3;
  h_obj[50] = 3.0*loc3;

  h_obj[51] = tsqrt6*A[11];

  /* Third off-diagonal block */
  loc2 = matr[2]*loc0;
  J_B[1] += loc2;
  J_B[3] -= loc2;

  loc2 = matr[1]*loc0;
  J_B[2] -= loc2;
  J_B[6] += loc2;

  loc2 = matr[0]*loc0;
  J_B[5] += loc2;
  J_B[7] -= loc2;

  loc2 = sqrt3*J_B[3];
  loc3 = sqrt6*J_B[6];
  loc4 = loc2 + loc3;

  A[0] = -J_B[0] - loc4;
  A[1] =  J_B[0] - loc4;
  A[2] = 2.0*loc2 - loc3;
  A[3] = 3.0*loc3;

  loc2 = sqrt3*J_B[4];
  loc3 = sqrt6*J_B[7];
  loc4 = loc2 + loc3;

  A[4] = -J_B[1] - loc4;
  A[5] =  J_B[1] - loc4;
  A[6] = 2.0*loc2 - loc3;
  A[7] = 3.0*loc3;

  loc2 = sqrt3*J_B[5];
  loc3 = sqrt6*J_B[8];
  loc4 = loc2 + loc3;

  A[8] = -J_B[2] - loc4;
  A[9] =  J_B[2] - loc4;
  A[10] = 2.0*loc2 - loc3;
  A[11] = 3.0*loc3;

  loc2 = sqrt3*A[4];
  loc3 = sqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[52] = -A[0] - loc4;
  h_obj[53] =  A[0] - loc4;
  h_obj[54] = 2.0*loc2 - loc3;
  h_obj[55] = 3.0*loc3;

  loc2 = sqrt3*A[5];
  loc3 = sqrt6*A[9];
  loc4 = loc2 + loc3;

  h_obj[56] = -A[1] - loc4;
  h_obj[57] =  A[1] - loc4;
  h_obj[58] = 2.0*loc2 - loc3;
  h_obj[59] = 3.0*loc3;

  loc2 = sqrt3*A[6];
  loc3 = sqrt6*A[10];
  loc4 = loc2 + loc3;

  h_obj[60] = -A[2] - loc4;
  h_obj[61] =  A[2] - loc4;
  h_obj[62] = 2.0*loc2 - loc3;
  h_obj[63] = 3.0*loc3;

  loc2 = sqrt3*A[7];
  loc3 = sqrt6*A[11];
  loc4 = loc2 + loc3;

  h_obj[64] = -A[3] - loc4;
  h_obj[65] =  A[3] - loc4;
  h_obj[66] = 2.0*loc2 - loc3;
  h_obj[67] = 3.0*loc3;

  /* Third block of rows */
  loc2 = matr[6]*f;
  loc3 = dg[6]*f;
  loc4 = dg[6]*g + loc2;

  J_A[0] = loc1 + dg[6]*(loc2 + loc4);
  J_A[1] = loc3*matr[7] + loc4*dg[7];
  J_A[2] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[7]*f;
  loc3 = dg[7]*f;
  loc4 = dg[7]*g + loc2;

  J_A[3] = loc1 + dg[7]*(loc2 + loc4);
  J_A[4] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[8]*f;
  loc4 = dg[8]*g + loc2;

  J_A[5] = loc1 + dg[8]*(loc2 + loc4);

  /* Third diagonal block */
  loc2 = sqrt3*J_A[1];
  loc3 = sqrt6*J_A[2];
  loc4 = loc2 + loc3;

  A[0] = -J_A[0] - loc4;
  A[1] =  J_A[0] - loc4;

  loc2 = sqrt3*J_A[3];
  loc3 = sqrt6*J_A[4];
  loc4 = loc2 + loc3;

  A[4] = -J_A[1] - loc4;
  A[5] =  J_A[1] - loc4;
  A[6] = 2.0*loc2 - loc3;

  loc2 = sqrt3*J_A[4];
  loc3 = sqrt6*J_A[5];
  loc4 = loc2 + loc3;

  A[8] = -J_A[2] - loc4;
  A[9] =  J_A[2] - loc4;
  A[10] = 2.0*loc2 - loc3;
  A[11] = 3.0*loc3;

  loc2 = sqrt3*A[4];
  loc3 = sqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[68] = -A[0] - loc4;
  h_obj[69] =  A[0] - loc4;
  h_obj[70] = 2.0*loc2 - loc3;
  h_obj[71] = 3.0*loc3;

  loc2 = sqrt3*A[5];
  loc3 = sqrt6*A[9];

  h_obj[72] = A[1] - loc2 - loc3;
  h_obj[73] = 2.0*loc2 - loc3;
  h_obj[74] = 3.0*loc3;

  loc3 = sqrt6*A[10];
  h_obj[75] = tsqrt3*A[6] - loc3;
  h_obj[76] = 3.0*loc3;

  h_obj[77] = tsqrt6*A[11];
  return;
}
