#include <math.h>
#include "fcn.h"

/*****************************************************************************/
/* This set of functions reference tetrahedral elements to an regular        */
/* tetrahedron.  The input are the coordinates in the following order:       */
/*      [x1 x2 x3 x4 y1 y2 y3 y4 z1 z2 z3 z4]                                */
/* A zero return value indicates success, while a nonzero value indicates    */
/* failure.                                                                  */
/*****************************************************************************/
/* These functions are optimized for local mesh smoothing where we only need */
/* the derivatives with respect to x1, y1, and z1.                           */
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

/*****************************************************************************/
/* Optimal derivative calculation courtesy of Paul Hovland (at least we      */
/* think it is optimal).  The original code provided was modified to         */
/* reduce the number of flops and intermediate variables, and improve the    */
/* locality of reference.                                                    */
/*                                                                           */
/* This requires 82 flops compared to 130 for the complete gradient.         */
/*****************************************************************************/

int g_fcnl_0(double *obj, double g_obj[3], const double x[12])
{
#if defined(NAIVE)
  static double g_full[12];

  if (g_fcn(obj, g_full, x)) return 1;
  g_obj[0] = g_full[0];
  g_obj[1] = g_full[4];
  g_obj[2] = g_full[8];

  return 0;
#elif defined(BETTER)
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

  loc1 = sqrt3*adj_m[1];
  loc2 = sqrt6*adj_m[2];
  loc3 = loc1 + loc2;
  g_obj[0] = -adj_m[0] - loc3;

  loc1 = sqrt3*adj_m[4];
  loc2 = sqrt6*adj_m[5];
  loc3 = loc1 + loc2;
  g_obj[1] = -adj_m[3] - loc3;

  loc1 = sqrt3*adj_m[7];
  loc2 = sqrt6*adj_m[8];
  loc3 = loc1 + loc2;
  g_obj[2] = -adj_m[6] - loc3;
  return 0;
#else
  static double my_x[12];

  my_x[0] = x[1];
  my_x[1] = x[3];
  my_x[2] = x[2];
  my_x[3] = x[0];

  my_x[4] = x[5];
  my_x[5] = x[7];
  my_x[6] = x[6];
  my_x[7] = x[4];

  my_x[8] = x[9];
  my_x[9] = x[11];
  my_x[10] = x[10];
  my_x[11] = x[8];
  return g_fcnl_3(obj, g_obj, my_x);
#endif
}

int g_fcnl_1(double *obj, double g_obj[3], const double x[12])
{
#if defined(NAIVE)
  static double g_full[12];

  if (g_fcn(obj, g_full, x)) return 1;
  g_obj[0] = g_full[1];
  g_obj[1] = g_full[5];
  g_obj[2] = g_full[9];

  return 0;
#elif defined(BETTER)
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

  loc1 = sqrt3*adj_m[1];
  loc2 = sqrt6*adj_m[2];
  loc3 = loc1 + loc2;
  g_obj[0] = adj_m[0] - loc3;

  loc1 = sqrt3*adj_m[4];
  loc2 = sqrt6*adj_m[5];
  loc3 = loc1 + loc2;
  g_obj[1] = adj_m[3] - loc3;

  loc1 = sqrt3*adj_m[7];
  loc2 = sqrt6*adj_m[8];
  loc3 = loc1 + loc2;
  g_obj[2] = adj_m[6] - loc3;
  return 0;
#else
  static double my_x[12];

  my_x[0] = x[0];
  my_x[1] = x[2];
  my_x[2] = x[3];
  my_x[3] = x[1];

  my_x[4] = x[4];
  my_x[5] = x[6];
  my_x[6] = x[7];
  my_x[7] = x[5];

  my_x[8] = x[8];
  my_x[9] = x[10];
  my_x[10] = x[11];
  my_x[11] = x[9];
  return g_fcnl_3(obj, g_obj, my_x);
#endif
}

int g_fcnl_2(double *obj, double g_obj[3], const double x[12])
{
#if defined(NAIVE)
  static double g_full[12];

  if (g_fcn(obj, g_full, x)) return 1;
  g_obj[0] = g_full[2];
  g_obj[1] = g_full[6];
  g_obj[2] = g_full[10];

  return 0;
#elif defined(BETTER)
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
  loc4 = a * pow(g, b);
  *obj = f * loc4;

  /* Calculate the derivative of the objective function. */
  f = 2.0*loc4;
  g = b*(*obj)/g; 

  adj_m[1] = matr[1]*f + loc2*g;
  adj_m[2] = matr[2]*f + loc3*g;

  loc1 = matr[0]*g;
  loc2 = matr[1]*g;
  loc3 = matr[2]*g;

  adj_m[4] = matr[4]*f + loc1*matr[8] - loc3*matr[6];
  adj_m[5] = matr[5]*f + loc2*matr[6] - loc1*matr[7];

  adj_m[7] = matr[7]*f + loc3*matr[3] - loc1*matr[5];
  adj_m[8] = matr[8]*f + loc1*matr[4] - loc2*matr[3];

  loc1 = sqrt3*adj_m[1];
  loc2 = sqrt6*adj_m[2];
  g_obj[0] = 2.0*loc1 - loc2;

  loc1 = sqrt3*adj_m[4];
  loc2 = sqrt6*adj_m[5];
  g_obj[1] = 2.0*loc1 - loc2;

  loc1 = sqrt3*adj_m[7];
  loc2 = sqrt6*adj_m[8];
  g_obj[2] = 2.0*loc1 - loc2;
  return 0;
#else
  static double my_x[12];

  my_x[0] = x[1];
  my_x[1] = x[0];
  my_x[2] = x[3];
  my_x[3] = x[2];

  my_x[4] = x[5];
  my_x[5] = x[4];
  my_x[6] = x[7];
  my_x[7] = x[6];

  my_x[8] = x[9];
  my_x[9] = x[8];
  my_x[10] = x[11];
  my_x[11] = x[10];
  return g_fcnl_3(obj, g_obj, my_x);
#endif
}

int g_fcnl_3(double *obj, double g_obj[3], const double x[12])
{
#if defined(NAIVE)
  static double g_full[12];

  if (g_fcn(obj, g_full, x)) return 1;
  g_obj[0] = g_full[3];
  g_obj[1] = g_full[7];
  g_obj[2] = g_full[11];

  return 0;
#else
  static double matr[9], f, g;
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
  loc4 = a * pow(g, b);
  *obj = f * loc4;

  /* Calculate the derivative of the objective function. */
  f = 2.0 * loc4;
  g = b * (*obj) / g; 

  g_obj[0] = tsqrt6*(f*matr[2] + g*loc3);

  loc1 = g*matr[0];
  loc2 = g*matr[1];

  g_obj[1] = tsqrt6*(f*matr[5] + loc2*matr[6] - loc1*matr[7]);
  g_obj[2] = tsqrt6*(f*matr[8] + loc1*matr[4] - loc2*matr[3]);
  return 0;
#endif
}

/*****************************************************************************/
/* The Hessian calculation is done by blocks.  Only the upper triangular     */
/* blocks are stored.  The results in the data is in the following order:    */
/*    [d1 b1 b2 d2 b3 d3 ]                                                   */
/* The matrices on the diagonal (d1-d3) each contain 10 elements, while the  */
/* off-diagonal elements (b1-b3) each contain 16 elements.                   */
/*                                                                           */
/* The code requires 119 flops compared to 575 for the complete hessian.     */
/*****************************************************************************/

int h_fcnl_0(double *obj, double g_obj[3], double h_obj[6], const double x[12])
{
#if defined(NAIVE)
  static double g_full[12];
  static double h_full[78];

  if (h_fcn(obj, g_full, h_full, x)) return 1;
  g_obj[0] = g_full[0];
  g_obj[1] = g_full[4];
  g_obj[2] = g_full[8];

  h_obj[0] = h_full[0];
  h_obj[1] = h_full[10];
  h_obj[2] = h_full[26];
  h_obj[3] = h_full[42];
  h_obj[4] = h_full[52];
  h_obj[5] = h_full[68];

  return 0;
#elif defined(BETTER)
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

  loc1 = sqrt3*adj_m[1];
  loc2 = sqrt6*adj_m[2];
  loc3 = loc1 + loc2;
  g_obj[0] = -adj_m[0] - loc3;

  loc1 = sqrt3*adj_m[4];
  loc2 = sqrt6*adj_m[5];
  loc3 = loc1 + loc2;
  g_obj[1] = -adj_m[3] - loc3;

  loc1 = sqrt3*adj_m[7];
  loc2 = sqrt6*adj_m[8];
  loc3 = loc1 + loc2;
  g_obj[2] = -adj_m[6] - loc3;

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

  loc2 = sqrt3*J_A[3];
  loc3 = sqrt6*J_A[4];
  loc4 = loc2 + loc3;

  A[4] = -J_A[1] - loc4;

  loc2 = sqrt3*J_A[4];
  loc3 = sqrt6*J_A[5];
  loc4 = loc2 + loc3;

  A[8] = -J_A[2] - loc4;

  loc2 = sqrt3*A[4];
  loc3 = sqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[0] = -A[0] - loc4;

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

  loc2 = sqrt3*J_B[4];
  loc3 = sqrt6*J_B[7];
  loc4 = loc2 + loc3;

  A[4] = -J_B[1] - loc4;

  loc2 = sqrt3*J_B[5];
  loc3 = sqrt6*J_B[8];
  loc4 = loc2 + loc3;

  A[8] = -J_B[2] - loc4;

  loc2 = sqrt3*A[4];
  loc3 = sqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[1] = -A[0] - loc4;

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

  loc2 = sqrt3*J_C[4];
  loc3 = sqrt6*J_C[7];
  loc4 = loc2 + loc3;

  A[4] = -J_C[1] - loc4;

  loc2 = sqrt3*J_C[5];
  loc3 = sqrt6*J_C[8];
  loc4 = loc2 + loc3;

  A[8] = -J_C[2] - loc4;

  loc2 = sqrt3*A[4];
  loc3 = sqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[2] = -A[0] - loc4;

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

  loc2 = sqrt3*J_A[3];
  loc3 = sqrt6*J_A[4];
  loc4 = loc2 + loc3;

  A[4] = -J_A[1] - loc4;

  loc2 = sqrt3*J_A[4];
  loc3 = sqrt6*J_A[5];
  loc4 = loc2 + loc3;

  A[8] = -J_A[2] - loc4;

  loc2 = sqrt3*A[4];
  loc3 = sqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[3] = -A[0] - loc4;

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

  loc2 = sqrt3*J_B[4];
  loc3 = sqrt6*J_B[7];
  loc4 = loc2 + loc3;

  A[4] = -J_B[1] - loc4;

  loc2 = sqrt3*J_B[5];
  loc3 = sqrt6*J_B[8];
  loc4 = loc2 + loc3;

  A[8] = -J_B[2] - loc4;

  loc2 = sqrt3*A[4];
  loc3 = sqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[4] = -A[0] - loc4;

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

  loc2 = sqrt3*J_A[3];
  loc3 = sqrt6*J_A[4];
  loc4 = loc2 + loc3;

  A[4] = -J_A[1] - loc4;

  loc2 = sqrt3*J_A[4];
  loc3 = sqrt6*J_A[5];
  loc4 = loc2 + loc3;

  A[8] = -J_A[2] - loc4;

  loc2 = sqrt3*A[4];
  loc3 = sqrt6*A[8];
  loc4 = loc2 + loc3;

  h_obj[5] = -A[0] - loc4;
  return 0;
#else
  static double my_x[12];

  my_x[0] = x[1];
  my_x[1] = x[3];
  my_x[2] = x[2];
  my_x[3] = x[0];

  my_x[4] = x[5];
  my_x[5] = x[7];
  my_x[6] = x[6];
  my_x[7] = x[4];

  my_x[8] = x[9];
  my_x[9] = x[11];
  my_x[10] = x[10];
  my_x[11] = x[8];
  return h_fcnl_3(obj, g_obj, h_obj, my_x);
#endif
}

int h_fcnl_1(double *obj, double g_obj[3], double h_obj[6], const double x[12])
{
#if defined(NAIVE)
  static double g_full[12];
  static double h_full[78];

  if (h_fcn(obj, g_full, h_full, x)) return 1;
  g_obj[0] = g_full[1];
  g_obj[1] = g_full[5];
  g_obj[2] = g_full[9];

  h_obj[0] = h_full[4];
  h_obj[1] = h_full[15];
  h_obj[2] = h_full[31];
  h_obj[3] = h_full[46];
  h_obj[4] = h_full[57];
  h_obj[5] = h_full[72];

  return 0;
#elif defined(BETTER)
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

  loc1 = sqrt3*adj_m[1];
  loc2 = sqrt6*adj_m[2];
  loc3 = loc1 + loc2;
  g_obj[0] = adj_m[0] - loc3;

  loc1 = sqrt3*adj_m[4];
  loc2 = sqrt6*adj_m[5];
  loc3 = loc1 + loc2;
  g_obj[1] = adj_m[3] - loc3;

  loc1 = sqrt3*adj_m[7];
  loc2 = sqrt6*adj_m[8];
  loc3 = loc1 + loc2;
  g_obj[2] = adj_m[6] - loc3;

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

  A[1] =  J_A[0] - loc4;

  loc2 = sqrt3*J_A[3];
  loc3 = sqrt6*J_A[4];
  loc4 = loc2 + loc3;

  A[5] =  J_A[1] - loc4;

  loc2 = sqrt3*J_A[4];
  loc3 = sqrt6*J_A[5];
  loc4 = loc2 + loc3;

  A[9] =  J_A[2] - loc4;

  loc2 = sqrt3*A[5];
  loc3 = sqrt6*A[9];

  h_obj[0] = A[1] - loc2 - loc3;

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

  A[1] =  J_B[0] - loc4;

  loc2 = sqrt3*J_B[4];
  loc3 = sqrt6*J_B[7];
  loc4 = loc2 + loc3;

  A[5] =  J_B[1] - loc4;

  loc2 = sqrt3*J_B[5];
  loc3 = sqrt6*J_B[8];
  loc4 = loc2 + loc3;

  A[9] =  J_B[2] - loc4;

  loc2 = sqrt3*A[5];
  loc3 = sqrt6*A[9];
  loc4 = loc2 + loc3;

  h_obj[1] =  A[1] - loc4;

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

  A[1] =  J_C[0] - loc4;

  loc2 = sqrt3*J_C[4];
  loc3 = sqrt6*J_C[7];
  loc4 = loc2 + loc3;

  A[5] =  J_C[1] - loc4;

  loc2 = sqrt3*J_C[5];
  loc3 = sqrt6*J_C[8];
  loc4 = loc2 + loc3;

  A[9] =  J_C[2] - loc4;

  loc2 = sqrt3*A[5];
  loc3 = sqrt6*A[9];
  loc4 = loc2 + loc3;

  h_obj[2] =  A[1] - loc4;

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

  A[1] =  J_A[0] - loc4;

  loc2 = sqrt3*J_A[3];
  loc3 = sqrt6*J_A[4];
  loc4 = loc2 + loc3;

  A[5] =  J_A[1] - loc4;

  loc2 = sqrt3*J_A[4];
  loc3 = sqrt6*J_A[5];
  loc4 = loc2 + loc3;

  A[9] =  J_A[2] - loc4;

  loc2 = sqrt3*A[5];
  loc3 = sqrt6*A[9];

  h_obj[3] = A[1] - loc2 - loc3;

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

  A[1] =  J_B[0] - loc4;

  loc2 = sqrt3*J_B[4];
  loc3 = sqrt6*J_B[7];
  loc4 = loc2 + loc3;

  A[5] =  J_B[1] - loc4;

  loc2 = sqrt3*J_B[5];
  loc3 = sqrt6*J_B[8];
  loc4 = loc2 + loc3;

  A[9] =  J_B[2] - loc4;

  loc2 = sqrt3*A[5];
  loc3 = sqrt6*A[9];
  loc4 = loc2 + loc3;

  h_obj[4] =  A[1] - loc4;

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

  A[1] =  J_A[0] - loc4;

  loc2 = sqrt3*J_A[3];
  loc3 = sqrt6*J_A[4];
  loc4 = loc2 + loc3;

  A[5] =  J_A[1] - loc4;

  loc2 = sqrt3*J_A[4];
  loc3 = sqrt6*J_A[5];
  loc4 = loc2 + loc3;

  A[9] =  J_A[2] - loc4;

  loc2 = sqrt3*A[5];
  loc3 = sqrt6*A[9];

  h_obj[5] = A[1] - loc2 - loc3;
  return 0;
#else
  static double my_x[12];

  my_x[0] = x[0];
  my_x[1] = x[2];
  my_x[2] = x[3];
  my_x[3] = x[1];

  my_x[4] = x[4];
  my_x[5] = x[6];
  my_x[6] = x[7];
  my_x[7] = x[5];

  my_x[8] = x[8];
  my_x[9] = x[10];
  my_x[10] = x[11];
  my_x[11] = x[9];
  return h_fcnl_3(obj, g_obj, h_obj, my_x);
#endif
}

int h_fcnl_2(double *obj, double g_obj[3], double h_obj[6], const double x[12])
{
#if defined(NAIVE)
  static double g_full[12];
  static double h_full[78];

  if (h_fcn(obj, g_full, h_full, x)) return 1;
  g_obj[0] = g_full[2];
  g_obj[1] = g_full[6];
  g_obj[2] = g_full[10];

  h_obj[0] = h_full[7];
  h_obj[1] = h_full[20];
  h_obj[2] = h_full[36];
  h_obj[3] = h_full[49];
  h_obj[4] = h_full[62];
  h_obj[5] = h_full[75];

  return 0;
#elif defined (BETTER)
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
  loc1 = a * pow(g, b);
  *obj = f * loc1;

  /* Calculate the derivative of the objective function. */
  f = 2.0*loc1;
  g = b*(*obj)/g; 

  dg[4] = matr[0]*matr[8] - matr[2]*matr[6];
  dg[5] = matr[1]*matr[6] - matr[0]*matr[7];
  dg[7] = matr[2]*matr[3] - matr[0]*matr[5];
  dg[8] = matr[0]*matr[4] - matr[1]*matr[3];

  adj_m[1] = matr[1]*f + dg[1]*g;
  adj_m[2] = matr[2]*f + dg[2]*g;
  adj_m[4] = matr[4]*f + dg[4]*g;
  adj_m[5] = matr[5]*f + dg[5]*g;
  adj_m[7] = matr[7]*f + dg[7]*g;
  adj_m[8] = matr[8]*f + dg[8]*g;

  loc1 = sqrt3*adj_m[1];
  loc2 = sqrt6*adj_m[2];
  g_obj[0] = 2.0*loc1 - loc2;

  loc1 = sqrt3*adj_m[4];
  loc2 = sqrt6*adj_m[5];
  g_obj[1] = 2.0*loc1 - loc2;

  loc1 = sqrt3*adj_m[7];
  loc2 = sqrt6*adj_m[8];
  g_obj[2] = 2.0*loc1 - loc2;

  loc0 = g;
  loc1 = f;
  f = f*b/loc4;
  g = g*bm1/loc4;

  /* First block of rows */
  loc2 = matr[1]*f;
  loc3 = dg[1]*f;
  loc4 = dg[1]*g + loc2;

  J_A[3] = loc1 + dg[1]*(loc2 + loc4);
  J_A[4] = loc3*matr[2] + loc4*dg[2];
  J_B[4] = loc3*matr[4] + loc4*dg[4];
  J_B[5] = loc3*matr[5] + loc4*dg[5];
  J_C[4] = loc3*matr[7] + loc4*dg[7];
  J_C[5] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[2]*f;
  loc3 = dg[2]*f;
  loc4 = dg[2]*g + loc2;

  J_A[5] = loc1 + dg[2]*(loc2 + loc4);
  J_B[7] = loc3*matr[4] + loc4*dg[4];
  J_B[8] = loc3*matr[5] + loc4*dg[5];
  J_C[7] = loc3*matr[7] + loc4*dg[7];
  J_C[8] = loc3*matr[8] + loc4*dg[8];

  /* First diagonal block */
  loc2 = sqrt3*J_A[3];
  loc3 = sqrt6*J_A[4];

  A[6] = 2.0*loc2 - loc3;

  loc2 = sqrt3*J_A[4];
  loc3 = sqrt6*J_A[5];

  A[10] = 2.0*loc2 - loc3;

  loc3 = sqrt6*A[10];
  h_obj[0] = tsqrt3*A[6] - loc3;

  /* First off-diagonal block */
  loc2 = matr[6]*loc0;
  J_B[5] += loc2;
  J_B[7] -= loc2;

  loc2 = sqrt3*J_B[4];
  loc3 = sqrt6*J_B[7];

  A[6] = 2.0*loc2 - loc3;

  loc2 = sqrt3*J_B[5];
  loc3 = sqrt6*J_B[8];

  A[10] = 2.0*loc2 - loc3;

  loc2 = sqrt3*A[6];
  loc3 = sqrt6*A[10];

  h_obj[1] = 2.0*loc2 - loc3;

  /* Second off-diagonal block */
  loc2 = matr[3]*loc0;
  J_C[5] -= loc2;
  J_C[7] += loc2;

  loc2 = sqrt3*J_C[4];
  loc3 = sqrt6*J_C[7];

  A[6] = 2.0*loc2 - loc3;

  loc2 = sqrt3*J_C[5];
  loc3 = sqrt6*J_C[8];

  A[10] = 2.0*loc2 - loc3;

  loc2 = sqrt3*A[6];
  loc3 = sqrt6*A[10];

  h_obj[2] = 2.0*loc2 - loc3;

  /* Second block of rows */
  loc2 = matr[4]*f;
  loc3 = dg[4]*f;
  loc4 = dg[4]*g + loc2;

  J_A[3] = loc1 + dg[4]*(loc2 + loc4);
  J_A[4] = loc3*matr[5] + loc4*dg[5];
  J_B[4] = loc3*matr[7] + loc4*dg[7];
  J_B[5] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[5]*f;
  loc3 = dg[5]*f;
  loc4 = dg[5]*g + loc2;

  J_A[5] = loc1 + dg[5]*(loc2 + loc4);
  J_B[7] = loc3*matr[7] + loc4*dg[7];
  J_B[8] = loc3*matr[8] + loc4*dg[8];

  /* Second diagonal block */
  loc2 = sqrt3*J_A[3];
  loc3 = sqrt6*J_A[4];

  A[6] = 2.0*loc2 - loc3;

  loc2 = sqrt3*J_A[4];
  loc3 = sqrt6*J_A[5];

  A[10] = 2.0*loc2 - loc3;

  loc3 = sqrt6*A[10];
  h_obj[3] = tsqrt3*A[6] - loc3;

  /* Third off-diagonal block */
  loc2 = matr[0]*loc0;
  J_B[5] += loc2;
  J_B[7] -= loc2;

  loc2 = sqrt3*J_B[4];
  loc3 = sqrt6*J_B[7];

  A[6] = 2.0*loc2 - loc3;

  loc2 = sqrt3*J_B[5];
  loc3 = sqrt6*J_B[8];

  A[10] = 2.0*loc2 - loc3;

  loc2 = sqrt3*A[6];
  loc3 = sqrt6*A[10];

  h_obj[4] = 2.0*loc2 - loc3;

  /* Third block of rows */
  loc2 = matr[7]*f;
  loc3 = dg[7]*f;
  loc4 = dg[7]*g + loc2;

  J_A[3] = loc1 + dg[7]*(loc2 + loc4);
  J_A[4] = loc3*matr[8] + loc4*dg[8];

  loc2 = matr[8]*f;
  loc4 = dg[8]*g + loc2;

  J_A[5] = loc1 + dg[8]*(loc2 + loc4);

  /* Third diagonal block */

  loc2 = sqrt3*J_A[3];
  loc3 = sqrt6*J_A[4];

  A[6] = 2.0*loc2 - loc3;

  loc2 = sqrt3*J_A[4];
  loc3 = sqrt6*J_A[5];

  A[10] = 2.0*loc2 - loc3;

  loc3 = sqrt6*A[10];
  h_obj[5] = tsqrt3*A[6] - loc3;
  return 0;
#else
  static double my_x[12];

  my_x[0] = x[1];
  my_x[1] = x[0];
  my_x[2] = x[3];
  my_x[3] = x[2];

  my_x[4] = x[5];
  my_x[5] = x[4];
  my_x[6] = x[7];
  my_x[7] = x[6];

  my_x[8] = x[9];
  my_x[9] = x[8];
  my_x[10] = x[11];
  my_x[11] = x[10];
  return h_fcnl_3(obj, g_obj, h_obj, my_x);
#endif
}

int h_fcnl_3(double *obj, double g_obj[3], double h_obj[6], const double x[12])
{
#if defined(NAIVE)
  static double g_full[12];
  static double h_full[78];

  if (h_fcn(obj, g_full, h_full, x)) return 1;
  g_obj[0] = g_full[3];
  g_obj[1] = g_full[7];
  g_obj[2] = g_full[11];

  h_obj[0] = h_full[9];
  h_obj[1] = h_full[25];
  h_obj[2] = h_full[41];
  h_obj[3] = h_full[51];
  h_obj[4] = h_full[67];
  h_obj[5] = h_full[77];

  return 0;
#else
  static double matr[9], f, dobj_df, dobj_dfdg;
  static double dg[9], g, dobj_dg, dobj_dgdg;
  static double loc1;

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

  dg[5] = matr[1]*matr[6] - matr[0]*matr[7];
  dg[8] = matr[0]*matr[4] - matr[1]*matr[3];

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] + 
      matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8];

  /* Calculate objective function. */
  loc1 = a*pow(g, b);
  *obj = f*loc1;

  /* Calculate constants required */
  g = 1.0 / g;
  dobj_df = 2.0 * loc1;
  dobj_dg = b * (*obj) * g; 
  dobj_dfdg = b * dobj_df * g;
  dobj_dgdg = bm1 * dobj_dg * g;

  /* Calculate the derivative of the objective function. */
  g_obj[0] = tsqrt6*(dobj_df*matr[2] + dobj_dg*dg[2]);
  g_obj[1] = tsqrt6*(dobj_df*matr[5] + dobj_dg*dg[5]);
  g_obj[2] = tsqrt6*(dobj_df*matr[8] + dobj_dg*dg[8]);

  /* Start of Hessian calculation */
  matr[2] *= dobj_dfdg;
  matr[5] *= dobj_dfdg;
  matr[8] *= dobj_dfdg;

  /* First block of rows */
  loc1 = dobj_dgdg*dg[2] + matr[2];
  h_obj[0] = 1.5*(dobj_df + dg[2]*(matr[2] + loc1));
  h_obj[1] = 1.5*(dg[2]*matr[5] + loc1*dg[5]);
  h_obj[2] = 1.5*(dg[2]*matr[8] + loc1*dg[8]);

  /* Second block of rows */
  loc1 = dobj_dgdg*dg[5] + matr[5];
  h_obj[3] = 1.5*(dobj_df + dg[5]*(matr[5] + loc1));
  h_obj[4] = 1.5*(dg[5]*matr[8] + loc1*dg[8]);

  /* Third block of rows */
  h_obj[5] = 1.5*(dobj_df + dg[8]*(2.0*matr[8] + dobj_dgdg*dg[8]));
  return 0;
#endif
}
