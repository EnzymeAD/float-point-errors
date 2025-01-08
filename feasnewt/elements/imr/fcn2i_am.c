#include <math.h>
#include "fcn.h"

/*****************************************************************************/
/* This set of functions reference triangular elements to an isoscles        */
/* triangle.  The input are the coordinates in the following order:          */
/*      [x1 x2 x3 y1 y2 y3]                                                  */
/* A zero return value indicates success, while a nonzero value indicates    */
/* failure.                                                                  */
/*****************************************************************************/
/* Not all compilers substitute out constants (especially the square root).  */
/* Therefore, they are substituted out manually.  The values below were      */
/* calculated on a solaris machine using long doubles. I believe they are    */
/* accurate.                                                                 */
/*****************************************************************************/
/* Constant fixed for the averaging metric for quadrilateral elements.       */
/*****************************************************************************/

#define a       1.25000000000000000000000000000e-01        /*  1.0/8.0       */
#define b      -1.00000000000000000000000000000e-00        /* -1.0/1.0       */
#define bm1    -2.00000000000000000000000000000e-00        /* -2.0/1.0       */

int o_fcn(double *obj, const double x[6])
{
  static double matr[4], f;
  static double g;

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = x[2] - x[0];

  matr[2] = x[4] - x[3];
  matr[3] = x[5] - x[3];

  /* Calculate det(M). */
  g = matr[0]*matr[3] - matr[1]*matr[2];
  if (g <= epsilon) { *obj = g; return 1; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + 
      matr[2]*matr[2] + matr[3]*matr[3];

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
/* This requires 34 flops.  The function only requires 17 flops.             */
/*****************************************************************************/

int g_fcn(double *obj, double g_obj[6], const double x[6])
{
  static double matr[4], f;
  static double adj_m[4], g;
  static double loc1;

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = x[2] - x[0];

  matr[2] = x[4] - x[3];
  matr[3] = x[5] - x[3];

  /* Calculate det(M). */
  g = matr[0]*matr[3] - matr[1]*matr[2];
  if (g <= epsilon) { *obj = g; return 1; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + 
      matr[2]*matr[2] + matr[3]*matr[3];
 
  /* Calculate objective function. */
  loc1 = a * pow(g, b);
  *obj = f * loc1;

  /* Calculate the derivative of the objective function. */
  f = 2.0*loc1;
  g = b*(*obj)/g; 

  adj_m[0] = matr[0]*f + matr[3]*g;
  adj_m[1] = matr[1]*f - matr[2]*g;
  adj_m[2] = matr[2]*f - matr[1]*g;
  adj_m[3] = matr[3]*f + matr[0]*g;

  g_obj[0] = -adj_m[0] - adj_m[1];
  g_obj[1] =  adj_m[0];
  g_obj[2] =  adj_m[1];

  g_obj[3] = -adj_m[2] - adj_m[3];
  g_obj[4] =  adj_m[2];
  g_obj[5] =  adj_m[3];
  return 0;
}

/*****************************************************************************/
/* The Hessian calculation is done by blocks.  Only the upper triangular     */
/* blocks are stored.  The results in the data is in the following order:    */
/*    [d1 b1 d2 ]                                                            */
/* The matrices on the diagonal (d1-d2) each contain 3 elements, while the   */
/* off-diagonal elements (b1) each contain 4 elements.                       */
/*                                                                           */
/* The code requires 96 flops.  The gradient evaluation needs 34 flops       */
/* and the function requires 17 flops.  We can further exploit the structure */
/* to eliminate some of the repeated calculations.                           */
/*****************************************************************************/

int h_fcn(double *obj, double g_obj[6], double h_obj[21], const double x[6])
{
  static double matr[4], f;
  static double adj_m[4], g;
  static double loc0, loc1, loc2, loc3, loc4;
  static double A[2], J_A[3], J_B[4];

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = x[2] - x[0];

  matr[2] = x[4] - x[3];
  matr[3] = x[5] - x[3];

  /* Calculate det(M). */
  g = matr[0]*matr[3] - matr[1]*matr[2];
  if (g <= epsilon) { *obj = g; return 1; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + 
      matr[2]*matr[2] + matr[3]*matr[3];

  loc4 = g;

  /* Calculate objective function. */
  loc1 = a * pow(g, b);
  *obj = f * loc1;

  /* Calculate the derivative of the objective function. */
  f = 2.0*loc1;
  g = b*(*obj)/g; 

  adj_m[0] = matr[0]*f + matr[3]*g;
  adj_m[1] = matr[1]*f - matr[2]*g;
  adj_m[2] = matr[2]*f - matr[1]*g;
  adj_m[3] = matr[3]*f + matr[0]*g;

  g_obj[0] = -adj_m[0] - adj_m[1];
  g_obj[1] =  adj_m[0];
  g_obj[2] =  adj_m[1];

  g_obj[3] = -adj_m[2] - adj_m[3];
  g_obj[4] =  adj_m[2];
  g_obj[5] =  adj_m[3];

  loc0 = g;
  loc1 = f;
  f = f*b/loc4;
  g = g*bm1/loc4;

  /* First block of rows */
  loc2 = matr[0]*f;
  loc3 = matr[3]*f;
  loc4 = matr[3]*g + loc2;

  J_A[0] = loc1 + matr[3]*(loc2 + loc4);
  J_A[1] = loc3*matr[1] - loc4*matr[2];
  J_B[0] = loc3*matr[2] - loc4*matr[1];
  J_B[1] = loc3*matr[3] + loc4*matr[0];

  loc2 =  matr[1]*f;
  loc3 = -matr[2]*f;
  loc4 = -matr[2]*g + loc2;

  J_A[2] = loc1 - matr[2]*(loc2 + loc4);
  J_B[2] = loc3*matr[2] - loc4*matr[1];
  J_B[3] = loc3*matr[3] + loc4*matr[0];

  /* First diagonal block */
  A[0] = -J_A[0] - J_A[1];
  A[1] = -J_A[1] - J_A[2];

  h_obj[0] = -A[0] - A[1];
  h_obj[1] =  A[0];
  h_obj[2] =  A[1];

  h_obj[3] =  J_A[0];
  h_obj[4] =  J_A[1];

  h_obj[5] =  J_A[2];

  /* First off-diagonal block */
  J_B[1] += loc0;
  J_B[2] -= loc0;

  A[0] = -J_B[0] - J_B[2];
  A[1] = -J_B[1] - J_B[3];

  h_obj[6] = -A[0] - A[1];
  h_obj[7] =  A[0];
  h_obj[8] =  A[1];

  h_obj[9]  = -J_B[0] - J_B[1];
  h_obj[10] =  J_B[0];
  h_obj[11] =  J_B[1];

  h_obj[12] = -J_B[2] - J_B[3];
  h_obj[13] =  J_B[2];
  h_obj[14] =  J_B[3];

  /* Second block of rows */
  loc2 = matr[2]*f;
  loc3 = -matr[1]*f;
  loc4 = -matr[1]*g + loc2;

  J_A[0] = loc1 - matr[1]*(loc2 + loc4);
  J_A[1] = loc3*matr[3] + loc4*matr[0];

  loc2 = matr[3]*f;
  loc4 = matr[0]*g + loc2;

  J_A[2] = loc1 + matr[0]*(loc2 + loc4);

  /* Second diagonal block */
  A[0] = -J_A[0] - J_A[1];
  A[1] = -J_A[1] - J_A[2];

  h_obj[15] = -A[0] - A[1];
  h_obj[16] =  A[0];
  h_obj[17] =  A[1];

  h_obj[18] = J_A[0];
  h_obj[19] = J_A[1];

  h_obj[20] = J_A[2];
  return 0;
}

void h_only(double h_obj[21], const double x[6])
{
  static double matr[4], obj, f, g;
  static double loc0, loc1, loc2, loc3, loc4;
  static double A[2], J_A[3], J_B[4];

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = x[2] - x[0];

  matr[2] = x[4] - x[3];
  matr[3] = x[5] - x[3];

  /* Calculate det(M). */
  g = matr[0]*matr[3] - matr[1]*matr[2];

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + 
      matr[2]*matr[2] + matr[3]*matr[3];

  loc4 = g;

  /* Calculate objective function. */
  loc1 = a * pow(g, b);
  obj = f * loc1;

  /* Calculate the derivative of the objective function. */
  f = 2.0*loc1;
  g = b*obj/g; 

  loc0 = g;
  loc1 = f;
  f = f*b/loc4;
  g = g*bm1/loc4;

  /* First block of rows */
  loc2 = matr[0]*f;
  loc3 = matr[3]*f;
  loc4 = matr[3]*g + loc2;

  J_A[0] = loc1 + matr[3]*(loc2 + loc4);
  J_A[1] = loc3*matr[1] - loc4*matr[2];
  J_B[0] = loc3*matr[2] - loc4*matr[1];
  J_B[1] = loc3*matr[3] + loc4*matr[0];

  loc2 =  matr[1]*f;
  loc3 = -matr[2]*f;
  loc4 = -matr[2]*g + loc2;

  J_A[2] = loc1 - matr[2]*(loc2 + loc4);
  J_B[2] = loc3*matr[2] - loc4*matr[1];
  J_B[3] = loc3*matr[3] + loc4*matr[0];

  /* First diagonal block */
  A[0] = -J_A[0] - J_A[1];
  A[1] = -J_A[1] - J_A[2];

  h_obj[0] = -A[0] - A[1];
  h_obj[1] =  A[0];
  h_obj[2] =  A[1];

  h_obj[3] =  J_A[0];
  h_obj[4] =  J_A[1];

  h_obj[5] =  J_A[2];

  /* First off-diagonal block */
  J_B[1] += loc0;
  J_B[2] -= loc0;

  A[0] = -J_B[0] - J_B[2];
  A[1] = -J_B[1] - J_B[3];

  h_obj[6] = -A[0] - A[1];
  h_obj[7] =  A[0];
  h_obj[8] =  A[1];

  h_obj[9]  = -J_B[0] - J_B[1];
  h_obj[10] =  J_B[0];
  h_obj[11] =  J_B[1];

  h_obj[12] = -J_B[2] - J_B[3];
  h_obj[13] =  J_B[2];
  h_obj[14] =  J_B[3];

  /* Second block of rows */
  loc2 = matr[2]*f;
  loc3 = -matr[1]*f;
  loc4 = -matr[1]*g + loc2;

  J_A[0] = loc1 - matr[1]*(loc2 + loc4);
  J_A[1] = loc3*matr[3] + loc4*matr[0];

  loc2 = matr[3]*f;
  loc4 = matr[0]*g + loc2;

  J_A[2] = loc1 + matr[0]*(loc2 + loc4);

  /* Second diagonal block */
  A[0] = -J_A[0] - J_A[1];
  A[1] = -J_A[1] - J_A[2];

  h_obj[15] = -A[0] - A[1];
  h_obj[16] =  A[0];
  h_obj[17] =  A[1];

  h_obj[18] = J_A[0];
  h_obj[19] = J_A[1];

  h_obj[20] = J_A[2];
  return;
}

