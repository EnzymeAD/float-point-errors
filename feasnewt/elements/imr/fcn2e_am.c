#include <math.h>
#include "fcn.h"

/*****************************************************************************/
/* This set of functions reference triangular elements to an equilateral     */
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

#define a       5.00000000000000000000000000000e-01        /*  1.0/2.0       */
#define b      -1.00000000000000000000000000000e-00        /* -1.0/1.0       */
#define bm1    -2.00000000000000000000000000000e-00        /* -2.0/1.0       */

#define sqrt3   5.77350269189625797959429519858e-01        /*  1.0/sqrt(3.0) */

#define tsqrt3  1.15470053837925159591885903972e+00        /*  2.0/sqrt(3.0) */

#define e1a 1.00000000000	/* 1.0 / 1.0             */
#define e1b 0.57735026919	/* 1.0 / sqrt(3.0)       */
#define e2b 0.33333333333	/* 1.0 / 3.0             */

int o_fcn(double *obj, const double x[6])
{
  static double matr[4], f;
  static double g;

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = (2.0*x[2] - x[1] - x[0])*sqrt3;

  matr[2] = x[3] - x[4];	/* Negate to obtain original */
  matr[3] = (2.0*x[5] - x[4] - x[3])*sqrt3;

  /* Calculate det(M). */
  g = matr[0]*matr[3] + matr[1]*matr[2];
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
/* This requires 46 flops.  The function only requires 23 flops.             */
/*****************************************************************************/

int g_fcn(double *obj, double g_obj[6], const double x[6])
{
  static double matr[4], f;
  static double adj_m[4], g;
  static double loc1;

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = (2.0*x[2] - x[1] - x[0])*sqrt3;

  matr[2] = x[3] - x[4];	/* Negate to obtain original */
  matr[3] = (2.0*x[5] - x[4] - x[3])*sqrt3;

  /* Calculate det(M). */
  g = matr[0]*matr[3] + matr[1]*matr[2];
  if (g <= epsilon) { *obj = g; return 1; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + 
      matr[2]*matr[2] + matr[3]*matr[3];
 
  /* Calculate objective function. */
  loc1 = a * pow(g, b);
  *obj = f * loc1;

  /* Calculate the derivative of the objective function. */
  f = 2.0 * loc1;
  g = b * (*obj) / g; 

  adj_m[0] = f*matr[0] + g*matr[3];
  adj_m[1] = -sqrt3*(f*matr[1] + g*matr[2]);
  adj_m[2] = f*matr[2] + g*matr[1];
  adj_m[3] = -sqrt3*(f*matr[3] + g*matr[0]);

  g_obj[0] = adj_m[1] - adj_m[0];
  g_obj[1] = adj_m[1] + adj_m[0];
  g_obj[2] = -2.0*adj_m[1];

  g_obj[3] = adj_m[3] + adj_m[2];
  g_obj[4] = adj_m[3] - adj_m[2];
  g_obj[5] = -2.0*adj_m[3];
  return 0;
}

/*****************************************************************************/
/* The Hessian calculation is done by blocks.  Only the upper triangular     */
/* blocks are stored.  The results in the data is in the following order:    */
/*    [d1 b1 d2 ]                                                            */
/* The matrices on the diagonal (d1-d2) each contain 3 elements, while the   */
/* off-diagonal elements (b1) each contain 4 elements.                       */
/*                                                                           */
/* The code requires 138 flops.  The gradient evaluation needs 46 flops      */
/* and the function requires 23 flops.                                       */
/*****************************************************************************/

int h_fcn(double *obj, double g_obj[6], double h_obj[21], const double x[6])
{
  static double matr[4], f, dobj_df, dobj_dfdg;
  static double adj_m[4], g, dobj_dg, dobj_dgdg;
  static double dg[4], loc1, loc2;
  static double J_A[3], J_B[4], J_C[3];

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = (2.0*x[2] - x[1] - x[0])*sqrt3;

  matr[2] = x[3] - x[4];	/* Negate to get back to original */
  matr[3] = (2.0*x[5] - x[4] - x[3])*sqrt3;

  /* Calculate det(M). */
  dg[0] = matr[3];
  dg[1] = matr[2];
  dg[2] = matr[1];
  dg[3] = matr[0];
  g = matr[0]*dg[0] + matr[1]*dg[1];
  if (g <= epsilon) { *obj = g; return 1; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + 
      matr[2]*matr[2] + matr[3]*matr[3];

  /* Calculate objective function. */
  loc1 = a * pow(g, b);
  *obj = f * loc1;

  /* Calculate required constants */
  g = 1.0 / g;

  dobj_df = 2.0 * loc1;
  dobj_dg = b * f * loc1 * g; 
  dobj_dfdg = b * dobj_df * g;
  dobj_dgdg = bm1 * dobj_dg * g;

  /* Gradient evaluation */
  adj_m[0] = dobj_df*matr[0] + dobj_dg*matr[3];
  adj_m[1] = -sqrt3*(dobj_df*matr[1] + dobj_dg*matr[2]);
  adj_m[2] = dobj_df*matr[2] + dobj_dg*matr[1];
  adj_m[3] = -sqrt3*(dobj_df*matr[3] + dobj_dg*matr[0]);

  g_obj[0] = adj_m[1] - adj_m[0];
  g_obj[1] = adj_m[1] + adj_m[0];
  g_obj[2] = -2.0*adj_m[1];

  g_obj[3] = adj_m[3] + adj_m[2];
  g_obj[4] = adj_m[3] - adj_m[2];
  g_obj[5] = -2.0*adj_m[3];

  /* Start of Hessian evaluation */
  matr[0] *= dobj_dfdg;
  matr[1] *= dobj_dfdg;
  matr[2] *= dobj_dfdg;
  matr[3] *= dobj_dfdg;

  /* First block of rows */
  loc1 = dobj_dgdg*dg[0] + matr[0];
  J_A[0] = dobj_df + dg[0]*(matr[0] + loc1);
  J_A[1] = e1b*(dg[0]*matr[1] + loc1*dg[1]);
  J_B[0] =      dg[0]*matr[2] + loc1*dg[2];
  J_B[1] = e1b*(dg[0]*matr[3] + loc1*dg[3] + dobj_dg);

  loc1 = dobj_dgdg*dg[1] + matr[1];
  J_A[2] = e2b*(dobj_df + dg[1]*(matr[1] + loc1));
  J_B[2] = e1b*(dg[1]*matr[2] + loc1*dg[2] + dobj_dg);
  J_B[3] = e2b*(dg[1]*matr[3] + loc1*dg[3]);

  loc1 = dobj_dgdg*dg[2] + matr[2];
  J_C[0] = dobj_df + dg[2]*(matr[2] + loc1);
  J_C[1] = e1b*(dg[2]*matr[3] + loc1*dg[3]);

  J_C[2] = e2b*(dobj_df + dg[3]*(2.0*matr[3] + dobj_dgdg*dg[3]));

  /* Assembly */

  loc1 = J_A[1] + J_A[0];
  loc2 = J_A[2] + J_A[1];

  h_obj[0] = loc2 + loc1;
  h_obj[1] = loc2 - loc1;
  h_obj[2] = -2.0*loc2;

  loc1 = J_A[1] - J_A[0];
  loc2 = J_A[2] - J_A[1];

  h_obj[3] = loc2 - loc1;
  h_obj[4] = -2.0*loc2;

  h_obj[5] = 4.0*J_A[2];

  loc1 = J_B[2] + J_B[0];
  loc2 = J_B[3] + J_B[1];

  h_obj[6] = loc2 - loc1;
  h_obj[7] = loc2 + loc1;
  h_obj[8] = -2.0*loc2;

  loc1 = J_B[2] - J_B[0];
  loc2 = J_B[3] - J_B[1];

  h_obj[9]  = loc2 - loc1;
  h_obj[10] = loc2 + loc1;
  h_obj[11] = -2.0*loc2;

  // loc1 = -2.0*J_B[2];
  // loc2 = -2.0*J_B[3];

  h_obj[12] =  2.0*(J_B[2] - J_B[3]); // loc2 - loc1;
  h_obj[13] = -2.0*(J_B[2] + J_B[3]); // loc2 + loc1;
  h_obj[14] =  4.0*J_B[3];            // -2.0*loc2;

  loc1 = J_C[1] - J_C[0];
  loc2 = J_C[2] - J_C[1];

  h_obj[15] = loc2 - loc1;
  h_obj[16] = loc2 + loc1;
  h_obj[17] = -2.0*loc2;

  loc1 = J_C[1] + J_C[0];
  loc2 = J_C[2] + J_C[1];

  h_obj[18] = loc2 + loc1;
  h_obj[19] = -2.0*loc2;

  h_obj[20] = 4.0*J_C[2];
  return 0;
}

void h_only(double h_obj[21], const double x[6])
{
  static double matr[4], f, dobj_df, dobj_dfdg;
  static double dg[4], g, dobj_dg, dobj_dgdg;
  static double loc1, loc2;
  static double J_A[3], J_B[4], J_C[3];

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = (2.0*x[2] - x[1] - x[0])*sqrt3;

  matr[2] = x[3] - x[4];	/* Negate to get back to original */
  matr[3] = (2.0*x[5] - x[4] - x[3])*sqrt3;

  /* Calculate det(M). */
  dg[0] = matr[3];
  dg[1] = matr[2];
  dg[2] = matr[1];
  dg[3] = matr[0];
  g = matr[0]*dg[0] + matr[1]*dg[1];

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + 
      matr[2]*matr[2] + matr[3]*matr[3];

  /* Calculate required constants */
  loc1 = a * pow(g, b);
  g = 1.0 / g;

  dobj_df = 2.0 * loc1;
  dobj_dg = b * f * loc1 * g; 
  dobj_dfdg = b * dobj_df * g;
  dobj_dgdg = bm1 * dobj_dg * g;

  /* Start of Hessian evaluation */
  matr[0] *= dobj_dfdg;
  matr[1] *= dobj_dfdg;
  matr[2] *= dobj_dfdg;
  matr[3] *= dobj_dfdg;

  /* First block of rows */
  loc1 = dobj_dgdg*dg[0] + matr[0];
  J_A[0] = dobj_df + dg[0]*(matr[0] + loc1);
  J_A[1] = e1b*(dg[0]*matr[1] + loc1*dg[1]);
  J_B[0] =      dg[0]*matr[2] + loc1*dg[2];
  J_B[1] = e1b*(dg[0]*matr[3] + loc1*dg[3] + dobj_dg);

  loc1 = dobj_dgdg*dg[1] + matr[1];
  J_A[2] = e2b*(dobj_df + dg[1]*(matr[1] + loc1));
  J_B[2] = e1b*(dg[1]*matr[2] + loc1*dg[2] + dobj_dg);
  J_B[3] = e2b*(dg[1]*matr[3] + loc1*dg[3]);

  loc1 = dobj_dgdg*dg[2] + matr[2];
  J_C[0] = dobj_df + dg[2]*(matr[2] + loc1);
  J_C[1] = e1b*(dg[2]*matr[3] + loc1*dg[3]);

  J_C[2] = e2b*(dobj_df + dg[3]*(2.0*matr[3] + dobj_dgdg*dg[3]));

  /* Assembly */

  loc1 = J_A[1] + J_A[0];
  loc2 = J_A[2] + J_A[1];

  h_obj[0] = loc2 + loc1;
  h_obj[1] = loc2 - loc1;
  h_obj[2] = -2.0*loc2;

  loc1 = J_A[1] - J_A[0];
  loc2 = J_A[2] - J_A[1];

  h_obj[3] = loc2 - loc1;
  h_obj[4] = -2.0*loc2;

  h_obj[5] = 4.0*J_A[2];

  loc1 = J_B[2] + J_B[0];
  loc2 = J_B[3] + J_B[1];

  h_obj[6] = loc2 - loc1;
  h_obj[7] = loc2 + loc1;
  h_obj[8] = -2.0*loc2;

  loc1 = J_B[2] - J_B[0];
  loc2 = J_B[3] - J_B[1];

  h_obj[9] = loc2 - loc1;
  h_obj[10] = loc2 + loc1;
  h_obj[11] = -2.0*loc2;

  loc1 = -2.0*J_B[2];
  loc2 = -2.0*J_B[3];

  h_obj[12] = loc2 - loc1;
  h_obj[13] = loc2 + loc1;
  h_obj[14] = -2.0*loc2;

  loc1 = J_C[1] - J_C[0];
  loc2 = J_C[2] - J_C[1];

  h_obj[15] = loc2 - loc1;
  h_obj[16] = loc2 + loc1;
  h_obj[17] = -2.0*loc2;

  loc1 = J_C[1] + J_C[0];
  loc2 = J_C[2] + J_C[1];

  h_obj[18] = loc2 + loc1;
  h_obj[19] = -2.0*loc2;

  h_obj[20] = 4.0*J_C[2];
  return;
}

