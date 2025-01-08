#include <math.h>
#include "fcn.h"

/*****************************************************************************/
/* This set of functions reference triangular elements to an equilateral     */
/* triangle.  The input are the coordinates in the following order:          */
/*      [x1 x2 x3 y1 y2 y3]                                                  */
/* A zero return value indicates success, while a nonzero value indicates    */
/* failure.                                                                  */
/*****************************************************************************/
/* These functions are optimized for local mesh smoothing where we only need */
/* the derivatives with respect to x1, y1.                                   */
/*****************************************************************************/
/* Not all compilers substitute out constants (especially the square root).  */
/* Therefore, they are substituted out manually.  The values below were      */
/* calculated on a solaris machine using long doubles. I believe they are    */
/* accurate.                                                                 */
/*****************************************************************************/

#define sqrt3   5.77350269189625797959429519858e-01        /*  1.0/sqrt(3.0) */
#define tsqrt3  1.15470053837925159591885903972e+00        /*  2.0/sqrt(3.0) */
#define fthirds 1.33333333333333333333333333333e+00        /*  4.0/3.0       */
#define a       5.00000000000000000000000000000e-01        /*  1.0/2.0       */
#define b      -1.00000000000000000000000000000e-00        /* -1.0/1.0       */
#define bm1    -2.00000000000000000000000000000e-00        /* -2.0/1.0       */

/*****************************************************************************/
/* Optimal derivative calculation courtesy of Paul Hovland (at least we      */
/* think it is optimal).  The original code provided was modified to         */
/* reduce the number of flops and intermediate variables, and improve the    */
/* locality of reference.                                                    */
/*                                                                           */
/* This requires 34 flops.  The function only requires 23 flops.             */
/*****************************************************************************/

int g_fcnl_0(double *obj, double g_obj[2], const double x[6])
{
  static double my_x[6];

  my_x[0] = x[1];
  my_x[1] = x[2];
  my_x[2] = x[0];

  my_x[3] = x[4];
  my_x[4] = x[5];
  my_x[5] = x[3];
  return g_fcnl_2(obj, g_obj, my_x);
}

int g_fcnl_1(double *obj, double g_obj[2], const double x[6])
{
  static double my_x[6];

  my_x[0] = x[2];
  my_x[1] = x[0];
  my_x[2] = x[1];

  my_x[3] = x[5];
  my_x[4] = x[3];
  my_x[5] = x[4];
  return g_fcnl_2(obj, g_obj, my_x);
}

int g_fcnl_2(double *obj, double g_obj[2], const double x[6])
{
  static double matr[4], f, g;
  static double loc1;

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = (2.0*x[2] - x[1] - x[0])*sqrt3;

  matr[2] = x[3] - x[4];	/* Negate to get back to original */
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

  g_obj[0] = tsqrt3*(f*matr[1] + g*matr[2]);
  g_obj[1] = tsqrt3*(f*matr[3] + g*matr[0]);
  return 0;
}

/*****************************************************************************/
/* The Hessian calculation is done by blocks.  Only the upper triangular     */
/* blocks are stored.  The results in the data is in the following order:    */
/*    [d1 b1 d2 ]                                                            */
/* The matrices on the diagonal (d1-d2) each contain 3 elements, while the   */
/* off-diagonal elements (b1) each contain 4 elements.                       */
/*                                                                           */
/* The code requires 56 flops.  The gradient evaluation needs 34 flops       */
/* and the function requires 23 flops.                                       */
/*****************************************************************************/

int h_fcnl_0(double *obj, double g_obj[2], double h_obj[3], const double x[6])
{
  static double my_x[6];

  my_x[0] = x[1];
  my_x[1] = x[2];
  my_x[2] = x[0];

  my_x[3] = x[4];
  my_x[4] = x[5];
  my_x[5] = x[3];
  return h_fcnl_2(obj, g_obj, h_obj, my_x);
}

int h_fcnl_1(double *obj, double g_obj[2], double h_obj[3], const double x[6])
{
  static double my_x[6];

  my_x[0] = x[2];
  my_x[1] = x[0];
  my_x[2] = x[1];

  my_x[3] = x[5];
  my_x[4] = x[3];
  my_x[5] = x[4];
  return h_fcnl_2(obj, g_obj, h_obj, my_x);
}

int h_fcnl_2(double *obj, double g_obj[2], double h_obj[3], const double x[6])
{
  static double matr[4], f, dobj_df, dobj_dfdg;
  static double dg[4], g, dobj_dg, dobj_dgdg;
  static double loc1;

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = (2.0*x[2] - x[1] - x[0])*sqrt3;

  matr[2] = x[3] - x[4];	/* Negate to get back to original */
  matr[3] = (2.0*x[5] - x[4] - x[3])*sqrt3;

  /* Calculate det(M). */
  dg[0] = matr[3];
  dg[1] = matr[2];
  g = matr[0]*dg[0] + matr[1]*dg[1];
  if (g <= epsilon) { *obj = g; return 1; }

  dg[3] = matr[0];

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + 
      matr[2]*matr[2] + matr[3]*matr[3];

  /* Calculate objective function. */
  loc1 = a * pow(g, b);
  *obj = f * loc1;

  /* Calculate the derivative of the objective function. */
  g = 1.0 / g;
  dobj_df = 2.0 * loc1;
  dobj_dg = b * (*obj) * g; 
  dobj_dfdg = b * dobj_df * g;
  dobj_dgdg = bm1 * dobj_dg * g;

  g_obj[0] = tsqrt3*(dobj_df*matr[1] + dobj_dg*matr[2]);
  g_obj[1] = tsqrt3*(dobj_df*matr[3] + dobj_dg*matr[0]);

  /* Start of Hessian evaluation */
  matr[1] *= dobj_dfdg;
  matr[3] *= dobj_dfdg;

  /* Assembly */
  loc1 = dobj_dgdg*dg[1] + matr[1];
  h_obj[0] = fthirds*(dobj_df + dg[1]*(matr[1] + loc1));
  h_obj[1] = fthirds*(dg[1]*matr[3] + loc1*dg[3]);
  h_obj[2] = fthirds*(dobj_df + dg[3]*(2.0*matr[3] + dobj_dgdg*dg[3]));
  return 0;
}

