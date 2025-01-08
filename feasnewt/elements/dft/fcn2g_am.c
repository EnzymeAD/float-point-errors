#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "fcn.h"

#define a       5.00000000000000000000000000000e-01        /*  1.0/2.0       */
#define b      -1.00000000000000000000000000000e-00        /* -1.0/1.0       */
#define bm1    -2.00000000000000000000000000000e-00        /* -2.0/1.0       */

#define d       1.0e-4					   /*  delta         */
#define fd2     4.0e-8					   /*  4.0*delta^2   */

#define sqrt3   5.77350269189625797959429519858e-01        /*  1.0/sqrt(3.0) */

/* Assume that w is upper triangular. The t variable is the multiplication   */
/* of the rotation matrices required to obtain this matrix.                  */

const double w[4] = {1, -sqrt3, 0, 2*sqrt3};
const double t[4] = {1, 0, 0, 1};

int o_fcn(double *obj, const double x[6])
{
  /* 43 operations */

  static double matr[4], f, t1, t2;
  static double matd[4], g;

  /* Calculate M = A*inv(W). */
  f       = x[1] - x[0];
  matr[0] = f*w[0];
  matr[1] = f*w[1] + (x[2] - x[0])*w[3];

  f       = x[4] - x[3];
  matr[2] = f*w[0];
  matr[3] = f*w[1] + (x[5] - x[3])*w[3];

  /* Calculate det(M). */
  t1 = matr[0]*matr[3] - matr[1]*matr[2];
  t2 = sqrt(t1*t1 + fd2);
  g = t1 + t2;

  matd[0] = matr[0] - t[0];
  matd[1] = matr[1] - t[1];
  matd[2] = matr[2] - t[2];
  matd[3] = matr[3] - t[3];

  /* Calculate norm(M). */
  f = matd[0]*matd[0] + matd[1]*matd[1] +
      matd[2]*matd[2] + matd[3]*matd[3];

  /* Calculate objective function. */
  (*obj) = a * f * pow(g, b);
  return 0;
}

/*****************************************************************************/
/* Optimal derivative calculation.                                           */
/*****************************************************************************/

int g_fcn(double *obj, double g_obj[6], const double x[6])
{
  static double matr[4], f, t1, t2;
  static double matd[4], g;
  static double adj_m[4], loc1;

  /* Calculate M = A*inv(W). */
  f       = x[1] - x[0];
  matr[0] = f*w[0];
  matr[1] = f*w[1] + (x[2] - x[0])*w[3];
  
  f       = x[4] - x[3];
  matr[2] = f*w[0];
  matr[3] = f*w[1] + (x[5] - x[3])*w[3];

  /* Calculate det(M). */
  t1 = matr[0]*matr[3] - matr[1]*matr[2];
  t2 = sqrt(t1*t1 + fd2);
  g = t1 + t2;

  matd[0] = matr[0] - t[0];
  matd[1] = matr[1] - t[1];
  matd[2] = matr[2] - t[2];
  matd[3] = matr[3] - t[3];

  /* Calculate norm(M). */
  f = matd[0]*matd[0] + matd[1]*matd[1] +
      matd[2]*matd[2] + matd[3]*matd[3];

  /* Calculate objective function. */
  loc1 = a * pow(g, b);
  *obj = f * loc1;

  /* Calculate the derivative of the objective function. */
  f = 2.0 * loc1;
  g = b * (*obj) / t2;

  adj_m[0] = f*matd[0] + g*matr[3];
  adj_m[1] = f*matd[1] - g*matr[2];
  adj_m[2] = f*matd[2] - g*matr[1];
  adj_m[3] = f*matd[3] + g*matr[0];

  g_obj[1] = w[0]*adj_m[0] + w[1]*adj_m[1];
  g_obj[2] =                 w[3]*adj_m[1];
  g_obj[0] = -g_obj[1] - g_obj[2];

  g_obj[4] = w[0]*adj_m[2] + w[1]*adj_m[3];
  g_obj[5] =                 w[3]*adj_m[3];
  g_obj[3] = -g_obj[4] - g_obj[5];
  return 0;
}

/*****************************************************************************/
/* The Hessian calculation is done by blocks.  Only the upper triangular     */
/* blocks are stored.  The results in the data is in the following order:    */
/*    [d1 b1 b2 d2 b3 d3 ]                                                   */
/* The matrices on the diagonal (d1-d3) each contain 10 elements, while the  */
/* off-diagonal elements (b1-b3) each contain 16 elements.                   */
/*****************************************************************************/

static void assemble_diag(double h_obj[6], const double J_A[3])
{
  static double A[6];

  A[1]  =  J_A[0]*w[0] + J_A[1]*w[1];
  A[2]  =                J_A[1]*w[3];
  A[0]  = -A[1] - A[2];

  A[4]  =  J_A[1]*w[0] + J_A[2]*w[1];
  A[5]  =                J_A[2]*w[3];
  A[3]  = -A[4] - A[5];

  h_obj[1] =  A[0]*w[0] + A[3]*w[1];
  h_obj[2] =              A[3]*w[3];
  h_obj[0] = -h_obj[1] - h_obj[2];

  h_obj[3] =  A[1]*w[0] + A[4]*w[1];
  h_obj[4] =              A[4]*w[3];

  h_obj[5] =              A[5]*w[3];
  return;
}

static void assemble_offdiag(double h_obj[9], const double J_B[4])
{
  static double A[6];

  A[1]  =  J_B[0]*w[0] + J_B[1]*w[1];
  A[2]  =                J_B[1]*w[3];
  A[0]  = -A[1] - A[2];

  A[4]  =  J_B[2]*w[0] + J_B[3]*w[1];
  A[5]  =                J_B[3]*w[3];
  A[3]  = -A[4] - A[5];

  h_obj[3] = A[0]*w[0] + A[3]*w[1];
  h_obj[4] = A[1]*w[0] + A[4]*w[1];
  h_obj[5] = A[2]*w[0] + A[5]*w[1];

  h_obj[6] = A[3]*w[3];
  h_obj[7] = A[4]*w[3];
  h_obj[8] = A[5]*w[3];

  h_obj[0] = -h_obj[3] - h_obj[6];
  h_obj[1] = -h_obj[4] - h_obj[7];
  h_obj[2] = -h_obj[5] - h_obj[8];
  return;
}

int h_fcn(double *obj, double g_obj[6], double h_obj[21], const double x[6])
{
  static double matr[4], f, t1, t2;
  static double matd[4], g, t3, loc1;
  static double adj_m[4], dg[4], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;
  static double J_A[3], J_B[4], J_D[3];

  /* Calculate M = A*inv(W). */
  f       = x[1] - x[0];
  matr[0] = f*w[0];
  matr[1] = f*w[1] + (x[2] - x[0])*w[3];
  
  f       = x[4] - x[3];
  matr[2] = f*w[0];
  matr[3] = f*w[1] + (x[5] - x[3])*w[3];

  /* Calculate det(M). */
  dg[0] = matr[3];
  dg[1] = -matr[2];
  dg[2] = -matr[1];
  dg[3] = matr[0];

  t1 = matr[0]*matr[3] - matr[1]*matr[2];
  t2 = t1*t1 + fd2;
  t3 = sqrt(t2);
  g = t1 + t3;

  matd[0] = matr[0] - t[0];
  matd[1] = matr[1] - t[1];
  matd[2] = matr[2] - t[2];
  matd[3] = matr[3] - t[3];

  /* Calculate norm(M). */
  f = matd[0]*matd[0] + matd[1]*matd[1] + 
      matd[2]*matd[2] + matd[3]*matd[3];

  /* Calculate objective function. */
  loc1 = a * pow(g, b);
  *obj = f * loc1;

  /* Calculate constants required */

  /* Calculate the derivative of the objective function. */
  t3 = 1.0 / t3;
  dobj_df = 2.0 * loc1;
  dobj_dg = b * (*obj) * t3; 
  dobj_dfdg = b * dobj_df * t3;
  dobj_dgdg = dobj_dg * (bm1*t3 + fd2/(t2*g));

  adj_m[0] = dobj_df*matd[0] + dobj_dg*matr[3];
  adj_m[1] = dobj_df*matd[1] - dobj_dg*matr[2];
  adj_m[2] = dobj_df*matd[2] - dobj_dg*matr[1];
  adj_m[3] = dobj_df*matd[3] + dobj_dg*matr[0];

  g_obj[1] = w[0]*adj_m[0] + w[1]*adj_m[1];
  g_obj[2] =                 w[3]*adj_m[1];
  g_obj[0] = -g_obj[1] - g_obj[2];

  g_obj[4] = w[0]*adj_m[2] + w[1]*adj_m[3];
  g_obj[5] =                 w[3]*adj_m[3];
  g_obj[3] = -g_obj[4] - g_obj[5];

  /* Start of the Hessian evaluation */
  matd[0] *= dobj_dfdg;
  matd[1] *= dobj_dfdg;
  matd[2] *= dobj_dfdg;
  matd[3] *= dobj_dfdg;

  /* Blocks for the Hessian construction */
  loc1 = dobj_dgdg*dg[0] + matd[0];
  J_A[0] = dobj_df + dg[0]*(matd[0] + loc1);
  J_A[1] = dg[0]*matd[1] + loc1*dg[1];
  J_B[0] = dg[0]*matd[2] + loc1*dg[2];
  J_B[1] = dg[0]*matd[3] + loc1*dg[3] + dobj_dg;

  loc1 = dobj_dgdg*dg[1] + matd[1];
  J_A[2] = dobj_df + dg[1]*(matd[1] + loc1);
  J_B[2] = dg[1]*matd[2] + loc1*dg[2] - dobj_dg;
  J_B[3] = dg[1]*matd[3] + loc1*dg[3];

  loc1 = dobj_dgdg*dg[2] + matd[2];
  J_D[0] = dobj_df + dg[2]*(matd[2] + loc1);
  J_D[1] = dg[2]*matd[3] + loc1*dg[3];

  J_D[2] = dobj_df + dg[3]*(2.0*matd[3] + dobj_dgdg*dg[3]);

  assemble_diag(h_obj, J_A);
  assemble_diag(h_obj + 15, J_D);
  assemble_offdiag(h_obj + 6, J_B);
  return 0;
}

void h_only(double h_obj[20], const double x[6])
{
  static double matr[4], f, t1, t2;
  static double matd[4], g, t3, loc1;
  static double dg[4], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;
  static double J_A[3], J_B[4], J_D[3];

  /* Calculate M = A*inv(W). */
  f       = x[1] - x[0];
  matr[0] = f*w[0];
  matr[1] = f*w[1] + (x[2] - x[0])*w[3];
  
  f       = x[4] - x[3];
  matr[2] = f*w[0];
  matr[3] = f*w[1] + (x[5] - x[3])*w[3];

  /* Calculate det(M). */
  dg[0] = matr[3];
  dg[1] = -matr[2];
  dg[2] = -matr[1];
  dg[3] = matr[0];

  t1 = matr[0]*matr[3] - matr[1]*matr[2];
  t2 = t1*t1 + fd2;
  t3 = sqrt(t2);
  g = t1 + t3;

  matd[0] = matr[0] - t[0];
  matd[1] = matr[1] - t[1];
  matd[2] = matr[2] - t[2];
  matd[3] = matr[3] - t[3];

  /* Calculate norm(M). */
  f = matd[0]*matd[0] + matd[1]*matd[1] + 
      matd[2]*matd[2] + matd[3]*matd[3];

  /* Calculate objective function. */
  loc1 = a * pow(g, b);

  /* Calculate constants required */

  /* Calculate the derivative of the objective function. */
  t3 = 1.0 / t3;
  dobj_df = 2.0 * loc1;
  dobj_dg = b * f * loc1 * t3; 
  dobj_dfdg = b * dobj_df * t3;
  dobj_dgdg = dobj_dg * (bm1*t3 + fd2/(t2*g));

  /* Start of the Hessian evaluation */
  matd[0] *= dobj_dfdg;
  matd[1] *= dobj_dfdg;
  matd[2] *= dobj_dfdg;
  matd[3] *= dobj_dfdg;

  /* Blocks for the Hessian construction */
  loc1 = dobj_dgdg*dg[0] + matd[0];
  J_A[0] = dobj_df + dg[0]*(matd[0] + loc1);
  J_A[1] = dg[0]*matd[1] + loc1*dg[1];
  J_B[0] = dg[0]*matd[2] + loc1*dg[2];
  J_B[1] = dg[0]*matd[3] + loc1*dg[3] + dobj_dg;

  loc1 = dobj_dgdg*dg[1] + matd[1];
  J_A[2] = dobj_df + dg[1]*(matd[1] + loc1);
  J_B[2] = dg[1]*matd[2] + loc1*dg[2] - dobj_dg;
  J_B[3] = dg[1]*matd[3] + loc1*dg[3];

  loc1 = dobj_dgdg*dg[2] + matd[2];
  J_D[0] = dobj_df + dg[2]*(matd[2] + loc1);
  J_D[1] = dg[2]*matd[3] + loc1*dg[3];

  J_D[2] = dobj_df + dg[3]*(2.0*matd[3] + dobj_dgdg*dg[3]);

  assemble_diag(h_obj, J_A);
  assemble_diag(h_obj + 15, J_D);
  assemble_offdiag(h_obj + 6, J_B);
  return;
}
