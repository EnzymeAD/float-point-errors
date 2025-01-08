#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "fcn.h"

#define a       5.00000000000000000000000000000e-01        /*  1.0/2.0       */
#define b      -6.66666666666666666666666666667e-01        /* -2.0/3.0       */
#define bm1    -1.66666666666666666666666666667e-00        /* -5.0/3.0       */

#define d       1.0e-4					   /*  delta         */
#define fd2     4.0e-8					   /*  4.0*delta^2   */

#define sqrt3   5.77350269189625797959429519858e-01        /*  1.0/sqrt(3.0) */
#define sqrt6   4.08248290463863052509822647505e-01        /*  1.0/sqrt(6.0) */

/* Assume that w is upper triangular. The t variable is the multiplication   */
/* of the rotation matrices required to obtain this matrix.                  */

const double my_w[6] = {1, -sqrt3, -sqrt6, 2*sqrt3, -sqrt6, 3*sqrt6};
const double my_t[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};

int o_fcn(double *obj, const double x[12], 
          const double w[6], const double t[9])
{
  /* 83 operations */

  static double matr[9], f, t1, t2;
  static double matd[9], g;

  /* Calculate M = A*inv(W). */
  f       = x[1] - x[0];
  g       = x[2] - x[0];
  t1      = x[3] - x[0];
  matr[0] = f*w[0];
  matr[1] = f*w[1] + g*w[3];
  matr[2] = f*w[2] + g*w[4] + t1*w[5];

  f       = x[5] - x[4];
  g       = x[6] - x[4];
  t1      = x[7] - x[4];
  matr[3] = f*w[0];
  matr[4] = f*w[1] + g*w[3];
  matr[5] = f*w[2] + g*w[4] + t1*w[5];

  f       = x[9] - x[8];
  g       = x[10] - x[8];
  t1      = x[11] - x[8];
  matr[6] = f*w[0];
  matr[7] = f*w[1] + g*w[3];
  matr[8] = f*w[2] + g*w[4] + t1*w[5];

  /* Calculate det(M). */
  t1 = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) +
       matr[1]*(matr[5]*matr[6] - matr[3]*matr[8]) +
       matr[2]*(matr[3]*matr[7] - matr[4]*matr[6]);
  t2 = sqrt(t1*t1 + fd2);
  g = t1 + t2;

  matd[0] = matr[0] - t[0];
  matd[1] = matr[1] - t[1];
  matd[2] = matr[2] - t[2];
  matd[3] = matr[3] - t[3];
  matd[4] = matr[4] - t[4];
  matd[5] = matr[5] - t[5];
  matd[6] = matr[6] - t[6];
  matd[7] = matr[7] - t[7];
  matd[8] = matr[8] - t[8];

  /* Calculate norm(M). */
  f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
      matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
      matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

  /* Calculate objective function. */
  (*obj) = a * f * pow(g, b);
  return 0;
}

/*****************************************************************************/
/* Optimal derivative calculation.                                           */
/*****************************************************************************/

int g_fcn(double *obj, double g_obj[12], const double x[12],
	  const double w[6], const double t[9])
{
  /* 161 operations (83 for function) */

  static double matr[9], f, t1, t2;
  static double matd[9], g;
  static double adj_m[9], loc1, loc2, loc3, loc4;

  /* Calculate M = A*inv(W). */
  f       = x[1] - x[0];
  g       = x[2] - x[0];
  t1      = x[3] - x[0];
  matr[0] = f*w[0];
  matr[1] = f*w[1] + g*w[3];
  matr[2] = f*w[2] + g*w[4] + t1*w[5];

  f       = x[5] - x[4];
  g       = x[6] - x[4];
  t1      = x[7] - x[4];
  matr[3] = f*w[0];
  matr[4] = f*w[1] + g*w[3];
  matr[5] = f*w[2] + g*w[4] + t1*w[5];

  f       = x[9] - x[8];
  g       = x[10] - x[8];
  t1      = x[11] - x[8];
  matr[6] = f*w[0];
  matr[7] = f*w[1] + g*w[3];
  matr[8] = f*w[2] + g*w[4] + t1*w[5];

  /* Calculate det(M). */
  loc1 = matr[4]*matr[8] - matr[5]*matr[7];
  loc2 = matr[5]*matr[6] - matr[3]*matr[8];
  loc3 = matr[3]*matr[7] - matr[4]*matr[6];
  t1 = matr[0]*loc1 + matr[1]*loc2 + matr[2]*loc3;
  t2 = sqrt(t1*t1 + fd2);
  g = t1 + t2;

  matd[0] = matr[0] - t[0];
  matd[1] = matr[1] - t[1];
  matd[2] = matr[2] - t[2];
  matd[3] = matr[3] - t[3];
  matd[4] = matr[4] - t[4];
  matd[5] = matr[5] - t[5];
  matd[6] = matr[6] - t[6];
  matd[7] = matr[7] - t[7];
  matd[8] = matr[8] - t[8];

  /* Calculate norm(M). */
  f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
      matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
      matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

  /* Calculate objective function. */
  loc4 = a * pow(g, b);
  *obj = f * loc4;

  /* Calculate the derivative of the objective function. */
  f = 2.0 * loc4;
  g = b * (*obj) / t2;

  adj_m[0] = f*matd[0] + g*loc1;
  adj_m[1] = f*matd[1] + g*loc2;
  adj_m[2] = f*matd[2] + g*loc3;

  loc1 = g*matr[0];
  loc2 = g*matr[1];
  loc3 = g*matr[2];

  adj_m[3] = f*matd[3] + loc3*matr[7] - loc2*matr[8];
  adj_m[4] = f*matd[4] + loc1*matr[8] - loc3*matr[6];
  adj_m[5] = f*matd[5] + loc2*matr[6] - loc1*matr[7];

  adj_m[6] = f*matd[6] + loc2*matr[5] - loc3*matr[4];
  adj_m[7] = f*matd[7] + loc3*matr[3] - loc1*matr[5];
  adj_m[8] = f*matd[8] + loc1*matr[4] - loc2*matr[3];

  g_obj[1] = w[0]*adj_m[0] + w[1]*adj_m[1] + w[2]*adj_m[2];
  g_obj[2] =                 w[3]*adj_m[1] + w[4]*adj_m[2];
  g_obj[3] =                                 w[5]*adj_m[2];
  g_obj[0] = -g_obj[1] - g_obj[2] - g_obj[3];

  g_obj[5] = w[0]*adj_m[3] + w[1]*adj_m[4] + w[2]*adj_m[5];
  g_obj[6] =                 w[3]*adj_m[4] + w[4]*adj_m[5];
  g_obj[7] =                                 w[5]*adj_m[5];
  g_obj[4] = -g_obj[5] - g_obj[6] - g_obj[7];

  g_obj[9]  = w[0]*adj_m[6] + w[1]*adj_m[7] + w[2]*adj_m[8];
  g_obj[10] =                 w[3]*adj_m[7] + w[4]*adj_m[8];
  g_obj[11] =                                 w[5]*adj_m[8];
  g_obj[8]  = -g_obj[9] - g_obj[10] - g_obj[11];
  return 0;
}

/*****************************************************************************/
/* The Hessian calculation is done by blocks.  Only the upper triangular     */
/* blocks are stored.  The results in the data is in the following order:    */
/*    [d1 b1 b2 d2 b3 d3 ]                                                   */
/* The matrices on the diagonal (d1-d3) each contain 10 elements, while the  */
/* off-diagonal elements (b1-b3) each contain 16 elements.                   */
/*****************************************************************************/

#ifndef USE_INLINE
static void assemble_diag(double h_obj[10], const double J_A[6],
			  const double w[6])
{
  static double A[12];

  A[1]  =  J_A[0]*w[0] + J_A[1]*w[1] + J_A[2]*w[2];
  A[2]  =                J_A[1]*w[3] + J_A[2]*w[4];
  A[3]  =                              J_A[2]*w[5];
  A[0]  = -A[1] - A[2] - A[3];

  A[5]  =  J_A[1]*w[0] + J_A[3]*w[1] + J_A[4]*w[2];
  A[6]  =                J_A[3]*w[3] + J_A[4]*w[4];
  A[7]  =                              J_A[4]*w[5];
  A[4]  = -A[5] - A[6] - A[7];

  A[9]  =  J_A[2]*w[0] + J_A[4]*w[1] + J_A[5]*w[2];
  A[10] =                J_A[4]*w[3] + J_A[5]*w[4];
  A[11] =                              J_A[5]*w[5];
  A[8]  = -A[9] - A[10] - A[11];

  h_obj[1] =  A[0]*w[0] + A[4]*w[1] + A[8]*w[2];
  h_obj[2] =              A[4]*w[3] + A[8]*w[4];
  h_obj[3] =                          A[8]*w[5];
  h_obj[0] = -h_obj[1] - h_obj[2] - h_obj[3];

  h_obj[4] =  A[1]*w[0] + A[5]*w[1] + A[9]*w[2];
  h_obj[5] =              A[5]*w[3] + A[9]*w[4];
  h_obj[6] =                          A[9]*w[5];

  h_obj[7] =              A[6]*w[3] + A[10]*w[4];
  h_obj[8] =                          A[10]*w[5];

  h_obj[9] =                          A[11]*w[5];
  return;
}

static void assemble_offdiag(double h_obj[16], const double J_B[9],
			     const double w[6])
{
  static double A[12];

  A[1]  =  J_B[0]*w[0] + J_B[1]*w[1] + J_B[2]*w[2];
  A[2]  =                J_B[1]*w[3] + J_B[2]*w[4];
  A[3]  =                              J_B[2]*w[5];
  A[0]  = -A[1] - A[2] - A[3];

  A[5]  =  J_B[3]*w[0] + J_B[4]*w[1] + J_B[5]*w[2];
  A[6]  =                J_B[4]*w[3] + J_B[5]*w[4];
  A[7]  =                              J_B[5]*w[5];
  A[4]  = -A[5] - A[6] - A[7];

  A[9]  =  J_B[6]*w[0] + J_B[7]*w[1] + J_B[8]*w[2];
  A[10] =                J_B[7]*w[3] + J_B[8]*w[4];
  A[11] =                              J_B[8]*w[5];
  A[8]  = -A[9] - A[10] - A[11];

  h_obj[4] =  A[0]*w[0] + A[4]*w[1] + A[8]*w[2];
  h_obj[5] =  A[1]*w[0] + A[5]*w[1] + A[9]*w[2];
  h_obj[6] =  A[2]*w[0] + A[6]*w[1] + A[10]*w[2];
  h_obj[7] =  A[3]*w[0] + A[7]*w[1] + A[11]*w[2];

  h_obj[8]  = A[4]*w[3] + A[8]*w[4];
  h_obj[9]  = A[5]*w[3] + A[9]*w[4];
  h_obj[10] = A[6]*w[3] + A[10]*w[4];
  h_obj[11] = A[7]*w[3] + A[11]*w[4];

  h_obj[12] = A[8]*w[5];
  h_obj[13] = A[9]*w[5];
  h_obj[14] = A[10]*w[5];
  h_obj[15] = A[11]*w[5];

  h_obj[0] = -h_obj[4] - h_obj[8] - h_obj[12];
  h_obj[1] = -h_obj[5] - h_obj[9] - h_obj[13];
  h_obj[2] = -h_obj[6] - h_obj[10] - h_obj[14];
  h_obj[3] = -h_obj[7] - h_obj[11] - h_obj[15];
  return;
}

int h_fcn(double *obj, double g_obj[12], double h_obj[78], const double x[12],
	  const double w[6], const double t[9])
{
  /* 767 operations (161 for gradient, 83 for function) */

  static double matr[9], f, t1, t2;
  static double matd[9], g, t3, loc1;
  static double adj_m[9], dg[9], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;
  static double J_A[6], J_B[10], J_C[10], J_D[6], J_E[10], J_F[6];

  /* Calculate M = A*inv(W). */
  f       = x[1] - x[0];
  g       = x[2] - x[0];
  t1      = x[3] - x[0];
  matr[0] = f*w[0];
  matr[1] = f*w[1] + g*w[3];
  matr[2] = f*w[2] + g*w[4] + t1*w[5];

  f       = x[5] - x[4];
  g       = x[6] - x[4];
  t1      = x[7] - x[4];
  matr[3] = f*w[0];
  matr[4] = f*w[1] + g*w[3];
  matr[5] = f*w[2] + g*w[4] + t1*w[5];

  f       = x[9] - x[8];
  g       = x[10] - x[8];
  t1      = x[11] - x[8];
  matr[6] = f*w[0];
  matr[7] = f*w[1] + g*w[3];
  matr[8] = f*w[2] + g*w[4] + t1*w[5];

  /* Calculate det(M). */
  dg[0] = matr[4]*matr[8] - matr[5]*matr[7];
  dg[1] = matr[5]*matr[6] - matr[3]*matr[8];
  dg[2] = matr[3]*matr[7] - matr[4]*matr[6];
  dg[3] = matr[2]*matr[7] - matr[1]*matr[8];
  dg[4] = matr[0]*matr[8] - matr[2]*matr[6];
  dg[5] = matr[1]*matr[6] - matr[0]*matr[7];
  dg[6] = matr[1]*matr[5] - matr[2]*matr[4];
  dg[7] = matr[2]*matr[3] - matr[0]*matr[5];
  dg[8] = matr[0]*matr[4] - matr[1]*matr[3];

  t1 = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];
  t2 = t1*t1 + fd2;
  t3 = sqrt(t2);
  g = t1 + t3;

  matd[0] = matr[0] - t[0];
  matd[1] = matr[1] - t[1];
  matd[2] = matr[2] - t[2];
  matd[3] = matr[3] - t[3];
  matd[4] = matr[4] - t[4];
  matd[5] = matr[5] - t[5];
  matd[6] = matr[6] - t[6];
  matd[7] = matr[7] - t[7];
  matd[8] = matr[8] - t[8];

  /* Calculate norm(M). */
  f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] + 
      matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
      matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

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

  adj_m[0] = dobj_df*matd[0] + dobj_dg*dg[0];
  adj_m[1] = dobj_df*matd[1] + dobj_dg*dg[1];
  adj_m[2] = dobj_df*matd[2] + dobj_dg*dg[2];
  adj_m[3] = dobj_df*matd[3] + dobj_dg*dg[3];
  adj_m[4] = dobj_df*matd[4] + dobj_dg*dg[4];
  adj_m[5] = dobj_df*matd[5] + dobj_dg*dg[5];
  adj_m[6] = dobj_df*matd[6] + dobj_dg*dg[6];
  adj_m[7] = dobj_df*matd[7] + dobj_dg*dg[7];
  adj_m[8] = dobj_df*matd[8] + dobj_dg*dg[8];

  g_obj[1] = w[0]*adj_m[0] + w[1]*adj_m[1] + w[2]*adj_m[2];
  g_obj[2] =                 w[3]*adj_m[1] + w[4]*adj_m[2];
  g_obj[3] =                                 w[5]*adj_m[2];
  g_obj[0] = -g_obj[1] - g_obj[2] - g_obj[3];
 
  g_obj[5] = w[0]*adj_m[3] + w[1]*adj_m[4] + w[2]*adj_m[5];
  g_obj[6] =                 w[3]*adj_m[4] + w[4]*adj_m[5];
  g_obj[7] =                                 w[5]*adj_m[5];
  g_obj[4] = -g_obj[5] - g_obj[6] - g_obj[7];

  g_obj[9]  = w[0]*adj_m[6] + w[1]*adj_m[7] + w[2]*adj_m[8];
  g_obj[10] =                 w[3]*adj_m[7] + w[4]*adj_m[8];
  g_obj[11] =                                 w[5]*adj_m[8];
  g_obj[8]  = -g_obj[9] - g_obj[10] - g_obj[11];

  /* Start of the Hessian evaluation */
  adj_m[0] = dobj_dg*matr[0]; matd[0] *= dobj_dfdg;
  adj_m[1] = dobj_dg*matr[1]; matd[1] *= dobj_dfdg;
  adj_m[2] = dobj_dg*matr[2]; matd[2] *= dobj_dfdg;
  adj_m[3] = dobj_dg*matr[3]; matd[3] *= dobj_dfdg;
  adj_m[4] = dobj_dg*matr[4]; matd[4] *= dobj_dfdg;
  adj_m[5] = dobj_dg*matr[5]; matd[5] *= dobj_dfdg;
  adj_m[6] = dobj_dg*matr[6]; matd[6] *= dobj_dfdg;
  adj_m[7] = dobj_dg*matr[7]; matd[7] *= dobj_dfdg;
  adj_m[8] = dobj_dg*matr[8]; matd[8] *= dobj_dfdg;

  /* Blocks for the Hessian construction */
  loc1 = dobj_dgdg*dg[0] + matd[0];
  J_A[0] = dobj_df + dg[0]*(matd[0] + loc1);
  J_A[1] = dg[0]*matd[1] + loc1*dg[1];
  J_A[2] = dg[0]*matd[2] + loc1*dg[2];
  J_B[0] = dg[0]*matd[3] + loc1*dg[3];
  J_B[1] = dg[0]*matd[4] + loc1*dg[4] + adj_m[8];
  J_B[2] = dg[0]*matd[5] + loc1*dg[5] - adj_m[7];
  J_C[0] = dg[0]*matd[6] + loc1*dg[6];
  J_C[1] = dg[0]*matd[7] + loc1*dg[7] - adj_m[5];
  J_C[2] = dg[0]*matd[8] + loc1*dg[8] + adj_m[4];

  loc1 = dobj_dgdg*dg[1] + matd[1];
  J_A[3] = dobj_df + dg[1]*(matd[1] + loc1);
  J_A[4] = dg[1]*matd[2] + loc1*dg[2];
  J_B[3] = dg[1]*matd[3] + loc1*dg[3] - adj_m[8];
  J_B[4] = dg[1]*matd[4] + loc1*dg[4];
  J_B[5] = dg[1]*matd[5] + loc1*dg[5] + adj_m[6];
  J_C[3] = dg[1]*matd[6] + loc1*dg[6] + adj_m[5];
  J_C[4] = dg[1]*matd[7] + loc1*dg[7];
  J_C[5] = dg[1]*matd[8] + loc1*dg[8] - adj_m[3];

  loc1 = dobj_dgdg*dg[2] + matd[2];
  J_A[5] = dobj_df + dg[2]*(matd[2] + loc1);
  J_B[6] = dg[2]*matd[3] + loc1*dg[3] + adj_m[7];
  J_B[7] = dg[2]*matd[4] + loc1*dg[4] - adj_m[6];
  J_B[8] = dg[2]*matd[5] + loc1*dg[5];
  J_C[6] = dg[2]*matd[6] + loc1*dg[6] - adj_m[4];
  J_C[7] = dg[2]*matd[7] + loc1*dg[7] + adj_m[3];
  J_C[8] = dg[2]*matd[8] + loc1*dg[8];

  loc1 = dobj_dgdg*dg[3] + matd[3];
  J_D[0] = dobj_df + dg[3]*(matd[3] + loc1);
  J_D[1] = dg[3]*matd[4] + loc1*dg[4];
  J_D[2] = dg[3]*matd[5] + loc1*dg[5];
  J_E[0] = dg[3]*matd[6] + loc1*dg[6];
  J_E[1] = dg[3]*matd[7] + loc1*dg[7] + adj_m[2];
  J_E[2] = dg[3]*matd[8] + loc1*dg[8] - adj_m[1];

  loc1 = dobj_dgdg*dg[4] + matd[4];
  J_D[3] = dobj_df + dg[4]*(matd[4] + loc1);
  J_D[4] = dg[4]*matd[5] + loc1*dg[5];
  J_E[3] = dg[4]*matd[6] + loc1*dg[6] - adj_m[2];
  J_E[4] = dg[4]*matd[7] + loc1*dg[7];
  J_E[5] = dg[4]*matd[8] + loc1*dg[8] + adj_m[0];

  loc1 = dobj_dgdg*dg[5] + matd[5];
  J_D[5] = dobj_df + dg[5]*(matd[5] + loc1);
  J_E[6] = dg[5]*matd[6] + loc1*dg[6] + adj_m[1];
  J_E[7] = dg[5]*matd[7] + loc1*dg[7] - adj_m[0];
  J_E[8] = dg[5]*matd[8] + loc1*dg[8];

  loc1 = dobj_dgdg*dg[6] + matd[6];
  J_F[0] = dobj_df + dg[6]*(matd[6] + loc1);
  J_F[1] = dg[6]*matd[7] + loc1*dg[7];
  J_F[2] = dg[6]*matd[8] + loc1*dg[8];

  loc1 = dobj_dgdg*dg[7] + matd[7];
  J_F[3] = dobj_df + dg[7]*(matd[7] + loc1);
  J_F[4] = dg[7]*matd[8] + loc1*dg[8];

  J_F[5] = dobj_df + dg[8]*(2.0*matd[8] + dobj_dgdg*dg[8]);

  assemble_diag(h_obj, J_A, w);
  assemble_offdiag(h_obj + 10, J_B, w);
  assemble_offdiag(h_obj + 26, J_C, w);

  assemble_diag(h_obj + 42, J_D, w);
  assemble_offdiag(h_obj + 52, J_E, w);

  assemble_diag(h_obj + 68, J_F, w);
  return 0;
}

void h_only(double h_obj[78], const double x[12],
	    const double w[6], const double t[9])
{
  /* 707 operations (161 for gradient, 83 for function) */

  static double matr[9], f, t1, t2;
  static double matd[9], g, t3, loc1;
  static double adj_m[9], dg[9], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;
  static double J_A[6], J_B[10], J_C[10], J_D[6], J_E[10], J_F[6];

  /* Calculate M = A*inv(W). */
  f       = x[1] - x[0];
  g       = x[2] - x[0];
  t1      = x[3] - x[0];
  matr[0] = f*w[0];
  matr[1] = f*w[1] + g*w[3];
  matr[2] = f*w[2] + g*w[4] + t1*w[5];

  f       = x[5] - x[4];
  g       = x[6] - x[4];
  t1      = x[7] - x[4];
  matr[3] = f*w[0];
  matr[4] = f*w[1] + g*w[3];
  matr[5] = f*w[2] + g*w[4] + t1*w[5];

  f       = x[9] - x[8];
  g       = x[10] - x[8];
  t1      = x[11] - x[8];
  matr[6] = f*w[0];
  matr[7] = f*w[1] + g*w[3];
  matr[8] = f*w[2] + g*w[4] + t1*w[5];

  /* Calculate det(M). */
  dg[0] = matr[4]*matr[8] - matr[5]*matr[7];
  dg[1] = matr[5]*matr[6] - matr[3]*matr[8];
  dg[2] = matr[3]*matr[7] - matr[4]*matr[6];
  dg[3] = matr[2]*matr[7] - matr[1]*matr[8];
  dg[4] = matr[0]*matr[8] - matr[2]*matr[6];
  dg[5] = matr[1]*matr[6] - matr[0]*matr[7];
  dg[6] = matr[1]*matr[5] - matr[2]*matr[4];
  dg[7] = matr[2]*matr[3] - matr[0]*matr[5];
  dg[8] = matr[0]*matr[4] - matr[1]*matr[3];

  t1 = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];
  t2 = t1*t1 + fd2;
  t3 = sqrt(t2);
  g = t1 + t3;

  matd[0] = matr[0] - t[0];
  matd[1] = matr[1] - t[1];
  matd[2] = matr[2] - t[2];
  matd[3] = matr[3] - t[3];
  matd[4] = matr[4] - t[4];
  matd[5] = matr[5] - t[5];
  matd[6] = matr[6] - t[6];
  matd[7] = matr[7] - t[7];
  matd[8] = matr[8] - t[8];

  /* Calculate norm(M). */
  f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] + 
      matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
      matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

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
  adj_m[0] = dobj_dg*matr[0]; matd[0] *= dobj_dfdg;
  adj_m[1] = dobj_dg*matr[1]; matd[1] *= dobj_dfdg;
  adj_m[2] = dobj_dg*matr[2]; matd[2] *= dobj_dfdg;
  adj_m[3] = dobj_dg*matr[3]; matd[3] *= dobj_dfdg;
  adj_m[4] = dobj_dg*matr[4]; matd[4] *= dobj_dfdg;
  adj_m[5] = dobj_dg*matr[5]; matd[5] *= dobj_dfdg;
  adj_m[6] = dobj_dg*matr[6]; matd[6] *= dobj_dfdg;
  adj_m[7] = dobj_dg*matr[7]; matd[7] *= dobj_dfdg;
  adj_m[8] = dobj_dg*matr[8]; matd[8] *= dobj_dfdg;

  /* Blocks for the Hessian construction */
  loc1 = dobj_dgdg*dg[0] + matd[0];
  J_A[0] = dobj_df + dg[0]*(matd[0] + loc1);
  J_A[1] = dg[0]*matd[1] + loc1*dg[1];
  J_A[2] = dg[0]*matd[2] + loc1*dg[2];
  J_B[0] = dg[0]*matd[3] + loc1*dg[3];
  J_B[1] = dg[0]*matd[4] + loc1*dg[4] + adj_m[8];
  J_B[2] = dg[0]*matd[5] + loc1*dg[5] - adj_m[7];
  J_C[0] = dg[0]*matd[6] + loc1*dg[6];
  J_C[1] = dg[0]*matd[7] + loc1*dg[7] - adj_m[5];
  J_C[2] = dg[0]*matd[8] + loc1*dg[8] + adj_m[4];

  loc1 = dobj_dgdg*dg[1] + matd[1];
  J_A[3] = dobj_df + dg[1]*(matd[1] + loc1);
  J_A[4] = dg[1]*matd[2] + loc1*dg[2];
  J_B[3] = dg[1]*matd[3] + loc1*dg[3] - adj_m[8];
  J_B[4] = dg[1]*matd[4] + loc1*dg[4];
  J_B[5] = dg[1]*matd[5] + loc1*dg[5] + adj_m[6];
  J_C[3] = dg[1]*matd[6] + loc1*dg[6] + adj_m[5];
  J_C[4] = dg[1]*matd[7] + loc1*dg[7];
  J_C[5] = dg[1]*matd[8] + loc1*dg[8] - adj_m[3];

  loc1 = dobj_dgdg*dg[2] + matd[2];
  J_A[5] = dobj_df + dg[2]*(matd[2] + loc1);
  J_B[6] = dg[2]*matd[3] + loc1*dg[3] + adj_m[7];
  J_B[7] = dg[2]*matd[4] + loc1*dg[4] - adj_m[6];
  J_B[8] = dg[2]*matd[5] + loc1*dg[5];
  J_C[6] = dg[2]*matd[6] + loc1*dg[6] - adj_m[4];
  J_C[7] = dg[2]*matd[7] + loc1*dg[7] + adj_m[3];
  J_C[8] = dg[2]*matd[8] + loc1*dg[8];

  loc1 = dobj_dgdg*dg[3] + matd[3];
  J_D[0] = dobj_df + dg[3]*(matd[3] + loc1);
  J_D[1] = dg[3]*matd[4] + loc1*dg[4];
  J_D[2] = dg[3]*matd[5] + loc1*dg[5];
  J_E[0] = dg[3]*matd[6] + loc1*dg[6];
  J_E[1] = dg[3]*matd[7] + loc1*dg[7] + adj_m[2];
  J_E[2] = dg[3]*matd[8] + loc1*dg[8] - adj_m[1];

  loc1 = dobj_dgdg*dg[4] + matd[4];
  J_D[3] = dobj_df + dg[4]*(matd[4] + loc1);
  J_D[4] = dg[4]*matd[5] + loc1*dg[5];
  J_E[3] = dg[4]*matd[6] + loc1*dg[6] - adj_m[2];
  J_E[4] = dg[4]*matd[7] + loc1*dg[7];
  J_E[5] = dg[4]*matd[8] + loc1*dg[8] + adj_m[0];

  loc1 = dobj_dgdg*dg[5] + matd[5];
  J_D[5] = dobj_df + dg[5]*(matd[5] + loc1);
  J_E[6] = dg[5]*matd[6] + loc1*dg[6] + adj_m[1];
  J_E[7] = dg[5]*matd[7] + loc1*dg[7] - adj_m[0];
  J_E[8] = dg[5]*matd[8] + loc1*dg[8];

  loc1 = dobj_dgdg*dg[6] + matd[6];
  J_F[0] = dobj_df + dg[6]*(matd[6] + loc1);
  J_F[1] = dg[6]*matd[7] + loc1*dg[7];
  J_F[2] = dg[6]*matd[8] + loc1*dg[8];

  loc1 = dobj_dgdg*dg[7] + matd[7];
  J_F[3] = dobj_df + dg[7]*(matd[7] + loc1);
  J_F[4] = dg[7]*matd[8] + loc1*dg[8];

  J_F[5] = dobj_df + dg[8]*(2.0*matd[8] + dobj_dgdg*dg[8]);

  assemble_diag(h_obj, J_A, w);
  assemble_offdiag(h_obj + 10, J_B, w);
  assemble_offdiag(h_obj + 26, J_C, w);

  assemble_diag(h_obj + 42, J_D, w);
  assemble_offdiag(h_obj + 52, J_E, w);

  assemble_diag(h_obj + 68, J_F, w);
  return;
}
#else
int h_fcn(double *obj, double g_obj[12], double h_obj[78], const double x[12],
	  const double w[6], const double t[9])
{
  /* 767 operations (161 for gradient, 83 for function) */

  static double matr[9], f, t1, t2;
  static double matd[9], g, t3, loc1;
  static double adj_m[9], dg[9], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;
  static double J_A[6], J_B[10], J_C[10], J_D[6], J_E[10], J_F[6];
  static double A[12];

  /* Calculate M = A*inv(W). */
  f       = x[1] - x[0];
  g       = x[2] - x[0];
  t1      = x[3] - x[0];
  matr[0] = f*w[0];
  matr[1] = f*w[1] + g*w[3];
  matr[2] = f*w[2] + g*w[4] + t1*w[5];

  f       = x[5] - x[4];
  g       = x[6] - x[4];
  t1      = x[7] - x[4];
  matr[3] = f*w[0];
  matr[4] = f*w[1] + g*w[3];
  matr[5] = f*w[2] + g*w[4] + t1*w[5];

  f       = x[9] - x[8];
  g       = x[10] - x[8];
  t1      = x[11] - x[8];
  matr[6] = f*w[0];
  matr[7] = f*w[1] + g*w[3];
  matr[8] = f*w[2] + g*w[4] + t1*w[5];

  /* Calculate det(M). */
  dg[0] = matr[4]*matr[8] - matr[5]*matr[7];
  dg[1] = matr[5]*matr[6] - matr[3]*matr[8];
  dg[2] = matr[3]*matr[7] - matr[4]*matr[6];
  dg[3] = matr[2]*matr[7] - matr[1]*matr[8];
  dg[4] = matr[0]*matr[8] - matr[2]*matr[6];
  dg[5] = matr[1]*matr[6] - matr[0]*matr[7];
  dg[6] = matr[1]*matr[5] - matr[2]*matr[4];
  dg[7] = matr[2]*matr[3] - matr[0]*matr[5];
  dg[8] = matr[0]*matr[4] - matr[1]*matr[3];

  t1 = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];
  t2 = t1*t1 + fd2;
  t3 = sqrt(t2);
  g = t1 + t3;

  matd[0] = matr[0] - t[0];
  matd[1] = matr[1] - t[1];
  matd[2] = matr[2] - t[2];
  matd[3] = matr[3] - t[3];
  matd[4] = matr[4] - t[4];
  matd[5] = matr[5] - t[5];
  matd[6] = matr[6] - t[6];
  matd[7] = matr[7] - t[7];
  matd[8] = matr[8] - t[8];

  /* Calculate norm(M). */
  f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] + 
      matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
      matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

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

  adj_m[0] = dobj_df*matd[0] + dobj_dg*dg[0];
  adj_m[1] = dobj_df*matd[1] + dobj_dg*dg[1];
  adj_m[2] = dobj_df*matd[2] + dobj_dg*dg[2];
  adj_m[3] = dobj_df*matd[3] + dobj_dg*dg[3];
  adj_m[4] = dobj_df*matd[4] + dobj_dg*dg[4];
  adj_m[5] = dobj_df*matd[5] + dobj_dg*dg[5];
  adj_m[6] = dobj_df*matd[6] + dobj_dg*dg[6];
  adj_m[7] = dobj_df*matd[7] + dobj_dg*dg[7];
  adj_m[8] = dobj_df*matd[8] + dobj_dg*dg[8];

  g_obj[1] = w[0]*adj_m[0] + w[1]*adj_m[1] + w[2]*adj_m[2];
  g_obj[2] =                 w[3]*adj_m[1] + w[4]*adj_m[2];
  g_obj[3] =                                 w[5]*adj_m[2];
  g_obj[0] = -g_obj[1] - g_obj[2] - g_obj[3];
 
  g_obj[5] = w[0]*adj_m[3] + w[1]*adj_m[4] + w[2]*adj_m[5];
  g_obj[6] =                 w[3]*adj_m[4] + w[4]*adj_m[5];
  g_obj[7] =                                 w[5]*adj_m[5];
  g_obj[4] = -g_obj[5] - g_obj[6] - g_obj[7];

  g_obj[9]  = w[0]*adj_m[6] + w[1]*adj_m[7] + w[2]*adj_m[8];
  g_obj[10] =                 w[3]*adj_m[7] + w[4]*adj_m[8];
  g_obj[11] =                                 w[5]*adj_m[8];
  g_obj[8]  = -g_obj[9] - g_obj[10] - g_obj[11];

  /* Start of the Hessian evaluation */
  adj_m[0] = dobj_dg*matr[0]; matd[0] *= dobj_dfdg;
  adj_m[1] = dobj_dg*matr[1]; matd[1] *= dobj_dfdg;
  adj_m[2] = dobj_dg*matr[2]; matd[2] *= dobj_dfdg;
  adj_m[3] = dobj_dg*matr[3]; matd[3] *= dobj_dfdg;
  adj_m[4] = dobj_dg*matr[4]; matd[4] *= dobj_dfdg;
  adj_m[5] = dobj_dg*matr[5]; matd[5] *= dobj_dfdg;
  adj_m[6] = dobj_dg*matr[6]; matd[6] *= dobj_dfdg;
  adj_m[7] = dobj_dg*matr[7]; matd[7] *= dobj_dfdg;
  adj_m[8] = dobj_dg*matr[8]; matd[8] *= dobj_dfdg;

  /* Blocks for the Hessian construction */
  loc1 = dobj_dgdg*dg[0] + matd[0];
  J_A[0] = dobj_df + dg[0]*(matd[0] + loc1);
  J_A[1] = dg[0]*matd[1] + loc1*dg[1];
  J_A[2] = dg[0]*matd[2] + loc1*dg[2];
  J_B[0] = dg[0]*matd[3] + loc1*dg[3];
  J_B[1] = dg[0]*matd[4] + loc1*dg[4] + adj_m[8];
  J_B[2] = dg[0]*matd[5] + loc1*dg[5] - adj_m[7];
  J_C[0] = dg[0]*matd[6] + loc1*dg[6];
  J_C[1] = dg[0]*matd[7] + loc1*dg[7] - adj_m[5];
  J_C[2] = dg[0]*matd[8] + loc1*dg[8] + adj_m[4];

  loc1 = dobj_dgdg*dg[1] + matd[1];
  J_A[3] = dobj_df + dg[1]*(matd[1] + loc1);
  J_A[4] = dg[1]*matd[2] + loc1*dg[2];
  J_B[3] = dg[1]*matd[3] + loc1*dg[3] - adj_m[8];
  J_B[4] = dg[1]*matd[4] + loc1*dg[4];
  J_B[5] = dg[1]*matd[5] + loc1*dg[5] + adj_m[6];
  J_C[3] = dg[1]*matd[6] + loc1*dg[6] + adj_m[5];
  J_C[4] = dg[1]*matd[7] + loc1*dg[7];
  J_C[5] = dg[1]*matd[8] + loc1*dg[8] - adj_m[3];

  loc1 = dobj_dgdg*dg[2] + matd[2];
  J_A[5] = dobj_df + dg[2]*(matd[2] + loc1);
  J_B[6] = dg[2]*matd[3] + loc1*dg[3] + adj_m[7];
  J_B[7] = dg[2]*matd[4] + loc1*dg[4] - adj_m[6];
  J_B[8] = dg[2]*matd[5] + loc1*dg[5];
  J_C[6] = dg[2]*matd[6] + loc1*dg[6] - adj_m[4];
  J_C[7] = dg[2]*matd[7] + loc1*dg[7] + adj_m[3];
  J_C[8] = dg[2]*matd[8] + loc1*dg[8];

  loc1 = dobj_dgdg*dg[3] + matd[3];
  J_D[0] = dobj_df + dg[3]*(matd[3] + loc1);
  J_D[1] = dg[3]*matd[4] + loc1*dg[4];
  J_D[2] = dg[3]*matd[5] + loc1*dg[5];
  J_E[0] = dg[3]*matd[6] + loc1*dg[6];
  J_E[1] = dg[3]*matd[7] + loc1*dg[7] + adj_m[2];
  J_E[2] = dg[3]*matd[8] + loc1*dg[8] - adj_m[1];

  loc1 = dobj_dgdg*dg[4] + matd[4];
  J_D[3] = dobj_df + dg[4]*(matd[4] + loc1);
  J_D[4] = dg[4]*matd[5] + loc1*dg[5];
  J_E[3] = dg[4]*matd[6] + loc1*dg[6] - adj_m[2];
  J_E[4] = dg[4]*matd[7] + loc1*dg[7];
  J_E[5] = dg[4]*matd[8] + loc1*dg[8] + adj_m[0];

  loc1 = dobj_dgdg*dg[5] + matd[5];
  J_D[5] = dobj_df + dg[5]*(matd[5] + loc1);
  J_E[6] = dg[5]*matd[6] + loc1*dg[6] + adj_m[1];
  J_E[7] = dg[5]*matd[7] + loc1*dg[7] - adj_m[0];
  J_E[8] = dg[5]*matd[8] + loc1*dg[8];

  loc1 = dobj_dgdg*dg[6] + matd[6];
  J_F[0] = dobj_df + dg[6]*(matd[6] + loc1);
  J_F[1] = dg[6]*matd[7] + loc1*dg[7];
  J_F[2] = dg[6]*matd[8] + loc1*dg[8];

  loc1 = dobj_dgdg*dg[7] + matd[7];
  J_F[3] = dobj_df + dg[7]*(matd[7] + loc1);
  J_F[4] = dg[7]*matd[8] + loc1*dg[8];

  J_F[5] = dobj_df + dg[8]*(2.0*matd[8] + dobj_dgdg*dg[8]);

  /* assemble_diag(h_obj, J_A, w); */
  A[1]  =  J_A[0]*w[0] + J_A[1]*w[1] + J_A[2]*w[2];
  A[2]  =                J_A[1]*w[3] + J_A[2]*w[4];
  A[3]  =                              J_A[2]*w[5];
  A[0]  = -A[1] - A[2] - A[3];

  A[5]  =  J_A[1]*w[0] + J_A[3]*w[1] + J_A[4]*w[2];
  A[6]  =                J_A[3]*w[3] + J_A[4]*w[4];
  A[7]  =                              J_A[4]*w[5];
  A[4]  = -A[5] - A[6] - A[7];

  A[9]  =  J_A[2]*w[0] + J_A[4]*w[1] + J_A[5]*w[2];
  A[10] =                J_A[4]*w[3] + J_A[5]*w[4];
  A[11] =                              J_A[5]*w[5];
  A[8]  = -A[9] - A[10] - A[11];

  h_obj[1] =  A[0]*w[0] + A[4]*w[1] + A[8]*w[2];
  h_obj[2] =              A[4]*w[3] + A[8]*w[4];
  h_obj[3] =                          A[8]*w[5];
  h_obj[0] = -h_obj[1] - h_obj[2] - h_obj[3];

  h_obj[4] =  A[1]*w[0] + A[5]*w[1] + A[9]*w[2];
  h_obj[5] =              A[5]*w[3] + A[9]*w[4];
  h_obj[6] =                          A[9]*w[5];

  h_obj[7] =              A[6]*w[3] + A[10]*w[4];
  h_obj[8] =                          A[10]*w[5];

  h_obj[9] =                          A[11]*w[5];

  /* assemble_offdiag(h_obj + 10, J_B, w); */
  A[1]  =  J_B[0]*w[0] + J_B[1]*w[1] + J_B[2]*w[2];
  A[2]  =                J_B[1]*w[3] + J_B[2]*w[4];
  A[3]  =                              J_B[2]*w[5];
  A[0]  = -A[1] - A[2] - A[3];

  A[5]  =  J_B[3]*w[0] + J_B[4]*w[1] + J_B[5]*w[2];
  A[6]  =                J_B[4]*w[3] + J_B[5]*w[4];
  A[7]  =                              J_B[5]*w[5];
  A[4]  = -A[5] - A[6] - A[7];

  A[9]  =  J_B[6]*w[0] + J_B[7]*w[1] + J_B[8]*w[2];
  A[10] =                J_B[7]*w[3] + J_B[8]*w[4];
  A[11] =                              J_B[8]*w[5];
  A[8]  = -A[9] - A[10] - A[11];

  h_obj[14] =  A[0]*w[0] + A[4]*w[1] + A[8]*w[2];
  h_obj[15] =  A[1]*w[0] + A[5]*w[1] + A[9]*w[2];
  h_obj[16] =  A[2]*w[0] + A[6]*w[1] + A[10]*w[2];
  h_obj[17] =  A[3]*w[0] + A[7]*w[1] + A[11]*w[2];

  h_obj[18] = A[4]*w[3] + A[8]*w[4];
  h_obj[19] = A[5]*w[3] + A[9]*w[4];
  h_obj[20] = A[6]*w[3] + A[10]*w[4];
  h_obj[21] = A[7]*w[3] + A[11]*w[4];

  h_obj[22] = A[8]*w[5];
  h_obj[23] = A[9]*w[5];
  h_obj[24] = A[10]*w[5];
  h_obj[25] = A[11]*w[5];

  h_obj[10] = -h_obj[14] - h_obj[18] - h_obj[22];
  h_obj[11] = -h_obj[15] - h_obj[19] - h_obj[23];
  h_obj[12] = -h_obj[16] - h_obj[20] - h_obj[24];
  h_obj[13] = -h_obj[17] - h_obj[21] - h_obj[25];

  /* assemble_offdiag(h_obj + 26, J_C, w); */
  A[1]  =  J_C[0]*w[0] + J_C[1]*w[1] + J_C[2]*w[2];
  A[2]  =                J_C[1]*w[3] + J_C[2]*w[4];
  A[3]  =                              J_C[2]*w[5];
  A[0]  = -A[1] - A[2] - A[3];

  A[5]  =  J_C[3]*w[0] + J_C[4]*w[1] + J_C[5]*w[2];
  A[6]  =                J_C[4]*w[3] + J_C[5]*w[4];
  A[7]  =                              J_C[5]*w[5];
  A[4]  = -A[5] - A[6] - A[7];

  A[9]  =  J_C[6]*w[0] + J_C[7]*w[1] + J_C[8]*w[2];
  A[10] =                J_C[7]*w[3] + J_C[8]*w[4];
  A[11] =                              J_C[8]*w[5];
  A[8]  = -A[9] - A[10] - A[11];

  h_obj[30] =  A[0]*w[0] + A[4]*w[1] + A[8]*w[2];
  h_obj[31] =  A[1]*w[0] + A[5]*w[1] + A[9]*w[2];
  h_obj[32] =  A[2]*w[0] + A[6]*w[1] + A[10]*w[2];
  h_obj[33] =  A[3]*w[0] + A[7]*w[1] + A[11]*w[2];

  h_obj[34] = A[4]*w[3] + A[8]*w[4];
  h_obj[35] = A[5]*w[3] + A[9]*w[4];
  h_obj[36] = A[6]*w[3] + A[10]*w[4];
  h_obj[37] = A[7]*w[3] + A[11]*w[4];

  h_obj[38] = A[8]*w[5];
  h_obj[39] = A[9]*w[5];
  h_obj[40] = A[10]*w[5];
  h_obj[41] = A[11]*w[5];

  h_obj[26] = -h_obj[30] - h_obj[34] - h_obj[38];
  h_obj[27] = -h_obj[31] - h_obj[35] - h_obj[39];
  h_obj[28] = -h_obj[32] - h_obj[36] - h_obj[40];
  h_obj[29] = -h_obj[33] - h_obj[37] - h_obj[41];

  /* assemble_diag(h_obj + 42, J_D, w); */
  A[1]  =  J_D[0]*w[0] + J_D[1]*w[1] + J_D[2]*w[2];
  A[2]  =                J_D[1]*w[3] + J_D[2]*w[4];
  A[3]  =                              J_D[2]*w[5];
  A[0]  = -A[1] - A[2] - A[3];

  A[5]  =  J_D[1]*w[0] + J_D[3]*w[1] + J_D[4]*w[2];
  A[6]  =                J_D[3]*w[3] + J_D[4]*w[4];
  A[7]  =                              J_D[4]*w[5];
  A[4]  = -A[5] - A[6] - A[7];

  A[9]  =  J_D[2]*w[0] + J_D[4]*w[1] + J_D[5]*w[2];
  A[10] =                J_D[4]*w[3] + J_D[5]*w[4];
  A[11] =                              J_D[5]*w[5];
  A[8]  = -A[9] - A[10] - A[11];

  h_obj[43] =  A[0]*w[0] + A[4]*w[1] + A[8]*w[2];
  h_obj[44] =              A[4]*w[3] + A[8]*w[4];
  h_obj[45] =                          A[8]*w[5];
  h_obj[42] = -h_obj[43] - h_obj[44] - h_obj[45];

  h_obj[46] =  A[1]*w[0] + A[5]*w[1] + A[9]*w[2];
  h_obj[47] =              A[5]*w[3] + A[9]*w[4];
  h_obj[48] =                          A[9]*w[5];

  h_obj[49] =              A[6]*w[3] + A[10]*w[4];
  h_obj[50] =                          A[10]*w[5];

  h_obj[51] =                          A[11]*w[5];

  /* assemble_offdiag(h_obj + 52, J_E, w); */
  A[1]  =  J_E[0]*w[0] + J_E[1]*w[1] + J_E[2]*w[2];
  A[2]  =                J_E[1]*w[3] + J_E[2]*w[4];
  A[3]  =                              J_E[2]*w[5];
  A[0]  = -A[1] - A[2] - A[3];

  A[5]  =  J_E[3]*w[0] + J_E[4]*w[1] + J_E[5]*w[2];
  A[6]  =                J_E[4]*w[3] + J_E[5]*w[4];
  A[7]  =                              J_E[5]*w[5];
  A[4]  = -A[5] - A[6] - A[7];

  A[9]  =  J_E[6]*w[0] + J_E[7]*w[1] + J_E[8]*w[2];
  A[10] =                J_E[7]*w[3] + J_E[8]*w[4];
  A[11] =                              J_E[8]*w[5];
  A[8]  = -A[9] - A[10] - A[11];

  h_obj[56] =  A[0]*w[0] + A[4]*w[1] + A[8]*w[2];
  h_obj[57] =  A[1]*w[0] + A[5]*w[1] + A[9]*w[2];
  h_obj[58] =  A[2]*w[0] + A[6]*w[1] + A[10]*w[2];
  h_obj[59] =  A[3]*w[0] + A[7]*w[1] + A[11]*w[2];

  h_obj[60] = A[4]*w[3] + A[8]*w[4];
  h_obj[61] = A[5]*w[3] + A[9]*w[4];
  h_obj[62] = A[6]*w[3] + A[10]*w[4];
  h_obj[63] = A[7]*w[3] + A[11]*w[4];

  h_obj[64] = A[8]*w[5];
  h_obj[65] = A[9]*w[5];
  h_obj[66] = A[10]*w[5];
  h_obj[67] = A[11]*w[5];

  h_obj[52] = -h_obj[56] - h_obj[60] - h_obj[64];
  h_obj[53] = -h_obj[57] - h_obj[61] - h_obj[65];
  h_obj[54] = -h_obj[58] - h_obj[62] - h_obj[66];
  h_obj[55] = -h_obj[59] - h_obj[63] - h_obj[67];
  
  /* assemble_diag(h_obj + 68, J_F, w); */
  A[1]  =  J_F[0]*w[0] + J_F[1]*w[1] + J_F[2]*w[2];
  A[2]  =                J_F[1]*w[3] + J_F[2]*w[4];
  A[3]  =                              J_F[2]*w[5];
  A[0]  = -A[1] - A[2] - A[3];

  A[5]  =  J_F[1]*w[0] + J_F[3]*w[1] + J_F[4]*w[2];
  A[6]  =                J_F[3]*w[3] + J_F[4]*w[4];
  A[7]  =                              J_F[4]*w[5];
  A[4]  = -A[5] - A[6] - A[7];

  A[9]  =  J_F[2]*w[0] + J_F[4]*w[1] + J_F[5]*w[2];
  A[10] =                J_F[4]*w[3] + J_F[5]*w[4];
  A[11] =                              J_F[5]*w[5];
  A[8]  = -A[9] - A[10] - A[11];

  h_obj[69] =  A[0]*w[0] + A[4]*w[1] + A[8]*w[2];
  h_obj[70] =              A[4]*w[3] + A[8]*w[4];
  h_obj[71] =                          A[8]*w[5];
  h_obj[68] = -h_obj[69] - h_obj[70] - h_obj[71];

  h_obj[72] =  A[1]*w[0] + A[5]*w[1] + A[9]*w[2];
  h_obj[73] =              A[5]*w[3] + A[9]*w[4];
  h_obj[74] =                          A[9]*w[5];

  h_obj[75] =              A[6]*w[3] + A[10]*w[4];
  h_obj[76] =                          A[10]*w[5];

  h_obj[77] =                          A[11]*w[5];
  return 0;
}

void h_only(double h_obj[78], const double x[12],
	    const double w[6], const double t[9])
{
  /* 707 operations (161 for gradient, 83 for function) */

  static double matr[9], f, t1, t2;
  static double matd[9], g, t3, loc1;
  static double adj_m[9], dg[9], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;
  static double J_A[6], J_B[10], J_C[10], J_D[6], J_E[10], J_F[6];
  static double A[12];

  /* Calculate M = A*inv(W). */
  f       = x[1] - x[0];
  g       = x[2] - x[0];
  t1      = x[3] - x[0];
  matr[0] = f*w[0];
  matr[1] = f*w[1] + g*w[3];
  matr[2] = f*w[2] + g*w[4] + t1*w[5];

  f       = x[5] - x[4];
  g       = x[6] - x[4];
  t1      = x[7] - x[4];
  matr[3] = f*w[0];
  matr[4] = f*w[1] + g*w[3];
  matr[5] = f*w[2] + g*w[4] + t1*w[5];

  f       = x[9] - x[8];
  g       = x[10] - x[8];
  t1      = x[11] - x[8];
  matr[6] = f*w[0];
  matr[7] = f*w[1] + g*w[3];
  matr[8] = f*w[2] + g*w[4] + t1*w[5];

  /* Calculate det(M). */
  dg[0] = matr[4]*matr[8] - matr[5]*matr[7];
  dg[1] = matr[5]*matr[6] - matr[3]*matr[8];
  dg[2] = matr[3]*matr[7] - matr[4]*matr[6];
  dg[3] = matr[2]*matr[7] - matr[1]*matr[8];
  dg[4] = matr[0]*matr[8] - matr[2]*matr[6];
  dg[5] = matr[1]*matr[6] - matr[0]*matr[7];
  dg[6] = matr[1]*matr[5] - matr[2]*matr[4];
  dg[7] = matr[2]*matr[3] - matr[0]*matr[5];
  dg[8] = matr[0]*matr[4] - matr[1]*matr[3];

  t1 = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];
  t2 = t1*t1 + fd2;
  t3 = sqrt(t2);
  g = t1 + t3;

  matd[0] = matr[0] - t[0];
  matd[1] = matr[1] - t[1];
  matd[2] = matr[2] - t[2];
  matd[3] = matr[3] - t[3];
  matd[4] = matr[4] - t[4];
  matd[5] = matr[5] - t[5];
  matd[6] = matr[6] - t[6];
  matd[7] = matr[7] - t[7];
  matd[8] = matr[8] - t[8];

  /* Calculate norm(M). */
  f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] + 
      matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
      matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

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
  adj_m[0] = dobj_dg*matr[0]; matd[0] *= dobj_dfdg;
  adj_m[1] = dobj_dg*matr[1]; matd[1] *= dobj_dfdg;
  adj_m[2] = dobj_dg*matr[2]; matd[2] *= dobj_dfdg;
  adj_m[3] = dobj_dg*matr[3]; matd[3] *= dobj_dfdg;
  adj_m[4] = dobj_dg*matr[4]; matd[4] *= dobj_dfdg;
  adj_m[5] = dobj_dg*matr[5]; matd[5] *= dobj_dfdg;
  adj_m[6] = dobj_dg*matr[6]; matd[6] *= dobj_dfdg;
  adj_m[7] = dobj_dg*matr[7]; matd[7] *= dobj_dfdg;
  adj_m[8] = dobj_dg*matr[8]; matd[8] *= dobj_dfdg;

  /* Blocks for the Hessian construction */
  loc1 = dobj_dgdg*dg[0] + matd[0];
  J_A[0] = dobj_df + dg[0]*(matd[0] + loc1);
  J_A[1] = dg[0]*matd[1] + loc1*dg[1];
  J_A[2] = dg[0]*matd[2] + loc1*dg[2];
  J_B[0] = dg[0]*matd[3] + loc1*dg[3];
  J_B[1] = dg[0]*matd[4] + loc1*dg[4] + adj_m[8];
  J_B[2] = dg[0]*matd[5] + loc1*dg[5] - adj_m[7];
  J_C[0] = dg[0]*matd[6] + loc1*dg[6];
  J_C[1] = dg[0]*matd[7] + loc1*dg[7] - adj_m[5];
  J_C[2] = dg[0]*matd[8] + loc1*dg[8] + adj_m[4];

  loc1 = dobj_dgdg*dg[1] + matd[1];
  J_A[3] = dobj_df + dg[1]*(matd[1] + loc1);
  J_A[4] = dg[1]*matd[2] + loc1*dg[2];
  J_B[3] = dg[1]*matd[3] + loc1*dg[3] - adj_m[8];
  J_B[4] = dg[1]*matd[4] + loc1*dg[4];
  J_B[5] = dg[1]*matd[5] + loc1*dg[5] + adj_m[6];
  J_C[3] = dg[1]*matd[6] + loc1*dg[6] + adj_m[5];
  J_C[4] = dg[1]*matd[7] + loc1*dg[7];
  J_C[5] = dg[1]*matd[8] + loc1*dg[8] - adj_m[3];

  loc1 = dobj_dgdg*dg[2] + matd[2];
  J_A[5] = dobj_df + dg[2]*(matd[2] + loc1);
  J_B[6] = dg[2]*matd[3] + loc1*dg[3] + adj_m[7];
  J_B[7] = dg[2]*matd[4] + loc1*dg[4] - adj_m[6];
  J_B[8] = dg[2]*matd[5] + loc1*dg[5];
  J_C[6] = dg[2]*matd[6] + loc1*dg[6] - adj_m[4];
  J_C[7] = dg[2]*matd[7] + loc1*dg[7] + adj_m[3];
  J_C[8] = dg[2]*matd[8] + loc1*dg[8];

  loc1 = dobj_dgdg*dg[3] + matd[3];
  J_D[0] = dobj_df + dg[3]*(matd[3] + loc1);
  J_D[1] = dg[3]*matd[4] + loc1*dg[4];
  J_D[2] = dg[3]*matd[5] + loc1*dg[5];
  J_E[0] = dg[3]*matd[6] + loc1*dg[6];
  J_E[1] = dg[3]*matd[7] + loc1*dg[7] + adj_m[2];
  J_E[2] = dg[3]*matd[8] + loc1*dg[8] - adj_m[1];

  loc1 = dobj_dgdg*dg[4] + matd[4];
  J_D[3] = dobj_df + dg[4]*(matd[4] + loc1);
  J_D[4] = dg[4]*matd[5] + loc1*dg[5];
  J_E[3] = dg[4]*matd[6] + loc1*dg[6] - adj_m[2];
  J_E[4] = dg[4]*matd[7] + loc1*dg[7];
  J_E[5] = dg[4]*matd[8] + loc1*dg[8] + adj_m[0];

  loc1 = dobj_dgdg*dg[5] + matd[5];
  J_D[5] = dobj_df + dg[5]*(matd[5] + loc1);
  J_E[6] = dg[5]*matd[6] + loc1*dg[6] + adj_m[1];
  J_E[7] = dg[5]*matd[7] + loc1*dg[7] - adj_m[0];
  J_E[8] = dg[5]*matd[8] + loc1*dg[8];

  loc1 = dobj_dgdg*dg[6] + matd[6];
  J_F[0] = dobj_df + dg[6]*(matd[6] + loc1);
  J_F[1] = dg[6]*matd[7] + loc1*dg[7];
  J_F[2] = dg[6]*matd[8] + loc1*dg[8];

  loc1 = dobj_dgdg*dg[7] + matd[7];
  J_F[3] = dobj_df + dg[7]*(matd[7] + loc1);
  J_F[4] = dg[7]*matd[8] + loc1*dg[8];

  J_F[5] = dobj_df + dg[8]*(2.0*matd[8] + dobj_dgdg*dg[8]);

  /* assemble_diag(h_obj, J_A, w); */
  A[1]  =  J_A[0]*w[0] + J_A[1]*w[1] + J_A[2]*w[2];
  A[2]  =                J_A[1]*w[3] + J_A[2]*w[4];
  A[3]  =                              J_A[2]*w[5];
  A[0]  = -A[1] - A[2] - A[3];

  A[5]  =  J_A[1]*w[0] + J_A[3]*w[1] + J_A[4]*w[2];
  A[6]  =                J_A[3]*w[3] + J_A[4]*w[4];
  A[7]  =                              J_A[4]*w[5];
  A[4]  = -A[5] - A[6] - A[7];

  A[9]  =  J_A[2]*w[0] + J_A[4]*w[1] + J_A[5]*w[2];
  A[10] =                J_A[4]*w[3] + J_A[5]*w[4];
  A[11] =                              J_A[5]*w[5];
  A[8]  = -A[9] - A[10] - A[11];

  h_obj[1] =  A[0]*w[0] + A[4]*w[1] + A[8]*w[2];
  h_obj[2] =              A[4]*w[3] + A[8]*w[4];
  h_obj[3] =                          A[8]*w[5];
  h_obj[0] = -h_obj[1] - h_obj[2] - h_obj[3];

  h_obj[4] =  A[1]*w[0] + A[5]*w[1] + A[9]*w[2];
  h_obj[5] =              A[5]*w[3] + A[9]*w[4];
  h_obj[6] =                          A[9]*w[5];

  h_obj[7] =              A[6]*w[3] + A[10]*w[4];
  h_obj[8] =                          A[10]*w[5];

  h_obj[9] =                          A[11]*w[5];

  /* assemble_offdiag(h_obj + 10, J_B, w); */
  A[1]  =  J_B[0]*w[0] + J_B[1]*w[1] + J_B[2]*w[2];
  A[2]  =                J_B[1]*w[3] + J_B[2]*w[4];
  A[3]  =                              J_B[2]*w[5];
  A[0]  = -A[1] - A[2] - A[3];

  A[5]  =  J_B[3]*w[0] + J_B[4]*w[1] + J_B[5]*w[2];
  A[6]  =                J_B[4]*w[3] + J_B[5]*w[4];
  A[7]  =                              J_B[5]*w[5];
  A[4]  = -A[5] - A[6] - A[7];

  A[9]  =  J_B[6]*w[0] + J_B[7]*w[1] + J_B[8]*w[2];
  A[10] =                J_B[7]*w[3] + J_B[8]*w[4];
  A[11] =                              J_B[8]*w[5];
  A[8]  = -A[9] - A[10] - A[11];

  h_obj[14] =  A[0]*w[0] + A[4]*w[1] + A[8]*w[2];
  h_obj[15] =  A[1]*w[0] + A[5]*w[1] + A[9]*w[2];
  h_obj[16] =  A[2]*w[0] + A[6]*w[1] + A[10]*w[2];
  h_obj[17] =  A[3]*w[0] + A[7]*w[1] + A[11]*w[2];

  h_obj[18] = A[4]*w[3] + A[8]*w[4];
  h_obj[19] = A[5]*w[3] + A[9]*w[4];
  h_obj[20] = A[6]*w[3] + A[10]*w[4];
  h_obj[21] = A[7]*w[3] + A[11]*w[4];

  h_obj[22] = A[8]*w[5];
  h_obj[23] = A[9]*w[5];
  h_obj[24] = A[10]*w[5];
  h_obj[25] = A[11]*w[5];

  h_obj[10] = -h_obj[14] - h_obj[18] - h_obj[22];
  h_obj[11] = -h_obj[15] - h_obj[19] - h_obj[23];
  h_obj[12] = -h_obj[16] - h_obj[20] - h_obj[24];
  h_obj[13] = -h_obj[17] - h_obj[21] - h_obj[25];

  /* assemble_offdiag(h_obj + 26, J_C, w); */
  A[1]  =  J_C[0]*w[0] + J_C[1]*w[1] + J_C[2]*w[2];
  A[2]  =                J_C[1]*w[3] + J_C[2]*w[4];
  A[3]  =                              J_C[2]*w[5];
  A[0]  = -A[1] - A[2] - A[3];

  A[5]  =  J_C[3]*w[0] + J_C[4]*w[1] + J_C[5]*w[2];
  A[6]  =                J_C[4]*w[3] + J_C[5]*w[4];
  A[7]  =                              J_C[5]*w[5];
  A[4]  = -A[5] - A[6] - A[7];

  A[9]  =  J_C[6]*w[0] + J_C[7]*w[1] + J_C[8]*w[2];
  A[10] =                J_C[7]*w[3] + J_C[8]*w[4];
  A[11] =                              J_C[8]*w[5];
  A[8]  = -A[9] - A[10] - A[11];

  h_obj[30] =  A[0]*w[0] + A[4]*w[1] + A[8]*w[2];
  h_obj[31] =  A[1]*w[0] + A[5]*w[1] + A[9]*w[2];
  h_obj[32] =  A[2]*w[0] + A[6]*w[1] + A[10]*w[2];
  h_obj[33] =  A[3]*w[0] + A[7]*w[1] + A[11]*w[2];

  h_obj[34] = A[4]*w[3] + A[8]*w[4];
  h_obj[35] = A[5]*w[3] + A[9]*w[4];
  h_obj[36] = A[6]*w[3] + A[10]*w[4];
  h_obj[37] = A[7]*w[3] + A[11]*w[4];

  h_obj[38] = A[8]*w[5];
  h_obj[39] = A[9]*w[5];
  h_obj[40] = A[10]*w[5];
  h_obj[41] = A[11]*w[5];

  h_obj[26] = -h_obj[30] - h_obj[34] - h_obj[38];
  h_obj[27] = -h_obj[31] - h_obj[35] - h_obj[39];
  h_obj[28] = -h_obj[32] - h_obj[36] - h_obj[40];
  h_obj[29] = -h_obj[33] - h_obj[37] - h_obj[41];

  /* assemble_diag(h_obj + 42, J_D, w); */
  A[1]  =  J_D[0]*w[0] + J_D[1]*w[1] + J_D[2]*w[2];
  A[2]  =                J_D[1]*w[3] + J_D[2]*w[4];
  A[3]  =                              J_D[2]*w[5];
  A[0]  = -A[1] - A[2] - A[3];

  A[5]  =  J_D[1]*w[0] + J_D[3]*w[1] + J_D[4]*w[2];
  A[6]  =                J_D[3]*w[3] + J_D[4]*w[4];
  A[7]  =                              J_D[4]*w[5];
  A[4]  = -A[5] - A[6] - A[7];

  A[9]  =  J_D[2]*w[0] + J_D[4]*w[1] + J_D[5]*w[2];
  A[10] =                J_D[4]*w[3] + J_D[5]*w[4];
  A[11] =                              J_D[5]*w[5];
  A[8]  = -A[9] - A[10] - A[11];

  h_obj[43] =  A[0]*w[0] + A[4]*w[1] + A[8]*w[2];
  h_obj[44] =              A[4]*w[3] + A[8]*w[4];
  h_obj[45] =                          A[8]*w[5];
  h_obj[42] = -h_obj[43] - h_obj[44] - h_obj[45];

  h_obj[46] =  A[1]*w[0] + A[5]*w[1] + A[9]*w[2];
  h_obj[47] =              A[5]*w[3] + A[9]*w[4];
  h_obj[48] =                          A[9]*w[5];

  h_obj[49] =              A[6]*w[3] + A[10]*w[4];
  h_obj[50] =                          A[10]*w[5];

  h_obj[51] =                          A[11]*w[5];

  /* assemble_offdiag(h_obj + 52, J_E, w); */
  A[1]  =  J_E[0]*w[0] + J_E[1]*w[1] + J_E[2]*w[2];
  A[2]  =                J_E[1]*w[3] + J_E[2]*w[4];
  A[3]  =                              J_E[2]*w[5];
  A[0]  = -A[1] - A[2] - A[3];

  A[5]  =  J_E[3]*w[0] + J_E[4]*w[1] + J_E[5]*w[2];
  A[6]  =                J_E[4]*w[3] + J_E[5]*w[4];
  A[7]  =                              J_E[5]*w[5];
  A[4]  = -A[5] - A[6] - A[7];

  A[9]  =  J_E[6]*w[0] + J_E[7]*w[1] + J_E[8]*w[2];
  A[10] =                J_E[7]*w[3] + J_E[8]*w[4];
  A[11] =                              J_E[8]*w[5];
  A[8]  = -A[9] - A[10] - A[11];

  h_obj[56] =  A[0]*w[0] + A[4]*w[1] + A[8]*w[2];
  h_obj[57] =  A[1]*w[0] + A[5]*w[1] + A[9]*w[2];
  h_obj[58] =  A[2]*w[0] + A[6]*w[1] + A[10]*w[2];
  h_obj[59] =  A[3]*w[0] + A[7]*w[1] + A[11]*w[2];

  h_obj[60] = A[4]*w[3] + A[8]*w[4];
  h_obj[61] = A[5]*w[3] + A[9]*w[4];
  h_obj[62] = A[6]*w[3] + A[10]*w[4];
  h_obj[63] = A[7]*w[3] + A[11]*w[4];

  h_obj[64] = A[8]*w[5];
  h_obj[65] = A[9]*w[5];
  h_obj[66] = A[10]*w[5];
  h_obj[67] = A[11]*w[5];

  h_obj[52] = -h_obj[56] - h_obj[60] - h_obj[64];
  h_obj[53] = -h_obj[57] - h_obj[61] - h_obj[65];
  h_obj[54] = -h_obj[58] - h_obj[62] - h_obj[66];
  h_obj[55] = -h_obj[59] - h_obj[63] - h_obj[67];
  
  /* assemble_diag(h_obj + 68, J_F, w); */
  A[1]  =  J_F[0]*w[0] + J_F[1]*w[1] + J_F[2]*w[2];
  A[2]  =                J_F[1]*w[3] + J_F[2]*w[4];
  A[3]  =                              J_F[2]*w[5];
  A[0]  = -A[1] - A[2] - A[3];

  A[5]  =  J_F[1]*w[0] + J_F[3]*w[1] + J_F[4]*w[2];
  A[6]  =                J_F[3]*w[3] + J_F[4]*w[4];
  A[7]  =                              J_F[4]*w[5];
  A[4]  = -A[5] - A[6] - A[7];

  A[9]  =  J_F[2]*w[0] + J_F[4]*w[1] + J_F[5]*w[2];
  A[10] =                J_F[4]*w[3] + J_F[5]*w[4];
  A[11] =                              J_F[5]*w[5];
  A[8]  = -A[9] - A[10] - A[11];

  h_obj[69] =  A[0]*w[0] + A[4]*w[1] + A[8]*w[2];
  h_obj[70] =              A[4]*w[3] + A[8]*w[4];
  h_obj[71] =                          A[8]*w[5];
  h_obj[68] = -h_obj[69] - h_obj[70] - h_obj[71];

  h_obj[72] =  A[1]*w[0] + A[5]*w[1] + A[9]*w[2];
  h_obj[73] =              A[5]*w[3] + A[9]*w[4];
  h_obj[74] =                          A[9]*w[5];

  h_obj[75] =              A[6]*w[3] + A[10]*w[4];
  h_obj[76] =                          A[10]*w[5];

  h_obj[77] =                          A[11]*w[5];
  return;
}
#endif

#ifdef MAIN
int main() {
  double x[12] = {0, 1, 0, 0,  
                  0, 0, 1, 0, 
                  0, 0, 0, 1};
  double obj;
  double g_obj[12]; 
  double h_obj[78]; 

  o_fcn(&obj, x);
  g_fcn(&obj, g_obj, x);
  h_fcn(&obj, g_obj, h_obj, x);

  printf("%5.4e\n", obj);
  printf("\n%5.4e\n", g_obj[0]);
  printf("%5.4e\n", g_obj[4]);
  printf("%5.4e\n", g_obj[8]);

  printf("\n%5.4e\n", g_obj[1]);
  printf("%5.4e\n", g_obj[5]);
  printf("%5.4e\n", g_obj[9]);

  printf("\n%5.4e\n", g_obj[2]);
  printf("%5.4e\n", g_obj[6]);
  printf("%5.4e\n", g_obj[10]);

  printf("\n%5.4e\n", g_obj[3]);
  printf("%5.4e\n", g_obj[7]);
  printf("%5.4e\n", g_obj[11]);

  printf("\n");
  printf("% 5.4e\t% 5.4e\t% 5.4e\n", h_obj[0], h_obj[10], h_obj[26]);
  printf("% 5.4e\t% 5.4e\t% 5.4e\n",      0.0, h_obj[42], h_obj[52]);
  printf("% 5.4e\t% 5.4e\t% 5.4e\n",      0.0,       0.0, h_obj[68]);

  printf("\n");
  printf("% 5.4e\t% 5.4e\t% 5.4e\n", h_obj[4], h_obj[15], h_obj[31]);
  printf("% 5.4e\t% 5.4e\t% 5.4e\n",      0.0, h_obj[46], h_obj[57]);
  printf("% 5.4e\t% 5.4e\t% 5.4e\n",      0.0,       0.0, h_obj[72]);

  printf("\n");
  printf("% 5.4e\t% 5.4e\t% 5.4e\n", h_obj[7], h_obj[20], h_obj[36]);
  printf("% 5.4e\t% 5.4e\t% 5.4e\n",      0.0, h_obj[49], h_obj[62]);
  printf("% 5.4e\t% 5.4e\t% 5.4e\n",      0.0,       0.0, h_obj[75]);

  printf("\n");
  printf("% 5.4e\t% 5.4e\t% 5.4e\n", h_obj[9], h_obj[25], h_obj[41]);
  printf("% 5.4e\t% 5.4e\t% 5.4e\n",      0.0, h_obj[51], h_obj[67]);
  printf("% 5.4e\t% 5.4e\t% 5.4e\n",      0.0,       0.0, h_obj[77]);
  return 0;
}
#endif
