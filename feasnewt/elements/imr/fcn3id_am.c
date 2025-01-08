#include <math.h>
#include "fcn.h"

/*****************************************************************************/
/* This set of functions reference tetrahedral elements to a right           */
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
/* Constant fixed for the averaging metric for hexahedraral elements.        */
/*****************************************************************************/

#define a       4.16666666666666666666666666666e-02        /*  1.0/24.0     */
#define b      -6.66666666666666666666666666667e-01        /* -2.0/ 3.0      */
#define bm1    -1.66666666666666666666666666667e-00        /* -5.0/ 3.0      */

int o_fcn(double *obj, const double x[12])
{
  static double matr[9], f;
  static double g;

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = x[2] - x[0];
  matr[2] = x[3] - x[0];

  matr[3] = x[5] - x[4];
  matr[4] = x[6] - x[4];
  matr[5] = x[7] - x[4];

  matr[6] = x[9] - x[8];
  matr[7] = x[10] - x[8];
  matr[8] = x[11] - x[8];

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
  (*obj) = a * f * pow(g, b);
  return 0;
}

/*****************************************************************************/
/* Optimal derivative calculation courtesy of Paul Hovland (at least we      */
/* think it is optimal).  The original code provided was modified to         */
/* reduce the number of flops and intermediate variables, and improve the    */
/* locality of reference.                                                    */
/*                                                                           */
/* This requires 94 flops.  The function only requires 43 flops.             */
/*****************************************************************************/

int g_fcn(double *obj, double g_obj[12], const double x[12])
{
  static double matr[9], f;
  static double adj_m[9], g;
  static double loc1, loc2, loc3, loc4;

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = x[2] - x[0];
  matr[2] = x[3] - x[0];

  matr[3] = x[5] - x[4];
  matr[4] = x[6] - x[4];
  matr[5] = x[7] - x[4];

  matr[6] = x[9] - x[8];
  matr[7] = x[10] - x[8];
  matr[8] = x[11] - x[8];

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
  loc4 = a * pow(g, b);
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

  g_obj[0] = -adj_m[0] - adj_m[1] - adj_m[2];
  g_obj[1] =  adj_m[0];
  g_obj[2] =  adj_m[1];
  g_obj[3] =  adj_m[2];

  g_obj[4] = -adj_m[3] - adj_m[4] - adj_m[5];
  g_obj[5] =  adj_m[3];
  g_obj[6] =  adj_m[4];
  g_obj[7] =  adj_m[5];

  g_obj[8]  = -adj_m[6] - adj_m[7] - adj_m[8];
  g_obj[9]  =  adj_m[6];
  g_obj[10] =  adj_m[7];
  g_obj[11] =  adj_m[8];
  return 0;
}

/*****************************************************************************/
/* The Hessian calculation is done by blocks.  Only the upper triangular     */
/* blocks are stored.  The results in the data is in the following order:    */
/*    [d1 b1 b2 d2 b3 d3 ]                                                   */
/* The matrices on the diagonal (d1-d3) each contain 10 elements, while the  */
/* off-diagonal elements (b1-b3) each contain 16 elements.                   */
/*                                                                           */
/* The code requires 364 flops.  The gradient evaluation needs 94 flops      */
/* and the function requires 43 flops.                                       */
/*****************************************************************************/

int h_fcn(double *obj, double g_obj[12], double h_obj[24], const double x[12])
{
  static double matr[9], f;
  static double adj_m[9], g;
  static double dg[9], loc0, loc1, loc2, loc3, loc4;
  static double J_A[6], J_B[9], J_C[9];

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = x[2] - x[0];
  matr[2] = x[3] - x[0];

  matr[3] = x[5] - x[4];
  matr[4] = x[6] - x[4];
  matr[5] = x[7] - x[4];

  matr[6] = x[9] - x[8];
  matr[7] = x[10] - x[8];
  matr[8] = x[11] - x[8];

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
  loc1 = a * pow(g, b);
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

  g_obj[0] = -adj_m[0] - adj_m[1] - adj_m[2];
  g_obj[1] =  adj_m[0];
  g_obj[2] =  adj_m[1];
  g_obj[3] =  adj_m[2];

  g_obj[4] = -adj_m[3] - adj_m[4] - adj_m[5];
  g_obj[5] =  adj_m[3];
  g_obj[6] =  adj_m[4];
  g_obj[7] =  adj_m[5];

  g_obj[8]  = -adj_m[6] - adj_m[7] - adj_m[8];
  g_obj[9]  =  adj_m[6];
  g_obj[10] =  adj_m[7];
  g_obj[11] =  adj_m[8];

  /* Start of the Hessian evaluation */
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
  h_obj[0] = J_A[0] + J_A[3] + J_A[5] + 2.0*(J_A[1] + J_A[2] + J_A[4]);
  h_obj[1] = J_A[0];
  h_obj[2] = J_A[3];
  h_obj[3] = J_A[5];

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

  h_obj[4] = J_B[0] + J_B[1] + J_B[2] +
             J_B[3] + J_B[4] + J_B[5] +
             J_B[6] + J_B[7] + J_B[8];
  h_obj[5] = J_B[0];
  h_obj[6] = J_B[4];
  h_obj[7] = J_B[8];

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

  h_obj[8] = J_C[0] + J_C[1] + J_C[2] +
             J_C[3] + J_C[4] + J_C[5] +
             J_C[6] + J_C[7] + J_C[8];
  h_obj[9] = J_C[0];
  h_obj[10] = J_C[4];
  h_obj[11] = J_C[8];

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
  h_obj[12] = J_A[0] + J_A[3] + J_A[5] + 2.0*(J_A[1] + J_A[2] + J_A[4]);
  h_obj[13] = J_A[0];
  h_obj[14] = J_A[3];
  h_obj[15] = J_A[5];

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

  h_obj[16] = J_B[0] + J_B[1] + J_B[2] +
              J_B[3] + J_B[4] + J_B[5] +
              J_B[6] + J_B[7] + J_B[8];
  h_obj[17] =  J_B[0];
  h_obj[18] =  J_B[4];
  h_obj[19] =  J_B[8];

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
  h_obj[20] = J_A[0] + J_A[3] + J_A[5] + 2.0*(J_A[1] + J_A[2] + J_A[4]);
  h_obj[21] = J_A[0];
  h_obj[22] = J_A[3];
  h_obj[23] = J_A[5];
  return 0;
}

void h_only(double h_obj[24], const double x[12])
{
  static double matr[9], obj;
  static double dg[9], f, g, loc0, loc1, loc2, loc3, loc4;
  static double J_A[6], J_B[9], J_C[9];

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = x[2] - x[0];
  matr[2] = x[3] - x[0];

  matr[3] = x[5] - x[4];
  matr[4] = x[6] - x[4];
  matr[5] = x[7] - x[4];

  matr[6] = x[9] - x[8];
  matr[7] = x[10] - x[8];
  matr[8] = x[11] - x[8];

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
  loc1 = a * pow(g, b);
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
  h_obj[0] = J_A[0] + J_A[3] + J_A[5] + 2.0*(J_A[1] + J_A[2] + J_A[4]);
  h_obj[1] = J_A[0];
  h_obj[2] = J_A[3];
  h_obj[3] = J_A[5];

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

  h_obj[4] = J_B[0] + J_B[1] + J_B[2] +
             J_B[3] + J_B[4] + J_B[5] +
             J_B[6] + J_B[7] + J_B[8];
  h_obj[5] = J_B[0];
  h_obj[6] = J_B[4];
  h_obj[7] = J_B[8];

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

  h_obj[8] = J_C[0] + J_C[1] + J_C[2] +
             J_C[3] + J_C[4] + J_C[5] +
             J_C[6] + J_C[7] + J_C[8];
  h_obj[9] = J_C[0];
  h_obj[10] = J_C[4];
  h_obj[11] = J_C[8];

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
  h_obj[12] = J_A[0] + J_A[3] + J_A[5] + 2.0*(J_A[1] + J_A[2] + J_A[4]);
  h_obj[13] = J_A[0];
  h_obj[14] = J_A[3];
  h_obj[15] = J_A[5];

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

  h_obj[16] = J_B[0] + J_B[1] + J_B[2] +
              J_B[3] + J_B[4] + J_B[5] +
              J_B[6] + J_B[7] + J_B[8];
  h_obj[17] =  J_B[0];
  h_obj[18] =  J_B[4];
  h_obj[19] =  J_B[8];

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
  h_obj[20] = J_A[0] + J_A[3] + J_A[5] + 2.0*(J_A[1] + J_A[2] + J_A[4]);
  h_obj[21] = J_A[0];
  h_obj[22] = J_A[3];
  h_obj[23] = J_A[5];
  return;
}
