#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "fcn.h"

#define a       5.00000000000000000000000000000e-01        /*  1.0/2.0       */
#define b      -1.00000000000000000000000000000e-00        /* -1.0/1.0       */
#define bm1    -2.00000000000000000000000000000e-00        /* -2.0/1.0       */

#define sqrt3   5.77350269189625797959429519858e-01        /*  1.0/sqrt(3.0) */

const double my_w[3] = {1, -sqrt3, 2*sqrt3};
const double my_t[4] = {1, 0, 0, 1};

/* Assume that w is upper triangular. The t variable is the multiplication   */
/* of the rotation matrices required to obtain this matrix.                  */

int o_fcn(double *obj, const double x[6], 
          const double w[3], const double t[4])
{
  static double matr[4], f, g;

  /* Calculate M = A*inv(W). */
  f       = x[1] - x[0];
  g       = x[2] - x[0];
  matr[0] = f*w[0];
  matr[1] = f*w[1] + g*w[2];

  f       = x[4] - x[3];
  g       = x[5] - x[3];
  matr[2] = f*w[0];
  matr[3] = f*w[1] + g*w[2];

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
/* Optimal derivative calculation.                                           */
/*****************************************************************************/

int g_fcn(double *obj, double g_obj[6], const double x[6],
	  const double w[3], const double t[4])
{
  static double matr[4], f;
  static double adj_m[4], g;
  static double loc1;

  /* Calculate M = A*inv(W). */
  f       = x[1] - x[0];
  g       = x[2] - x[0];
  matr[0] = f*w[0];
  matr[1] = f*w[1] + g*w[2];

  f       = x[4] - x[3];
  g       = x[5] - x[3];
  matr[2] = f*w[0];
  matr[3] = f*w[1] + g*w[2];

  /* Calculate det(M). */
  g = matr[0]*matr[3] - matr[1]*matr[2];
  if (g <= epsilon) { *obj = g; return 1; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] +
      matr[2]*matr[2] + matr[3]*matr[3];

  /* Calculate objective function. */
  loc4 = a * pow(g, b);
  *obj = f * loc4;

  /* Calculate the derivative of the objective function. */
  f = 2.0 * loc4;
  g = b * (*obj) / g;

  adj_m[0] = f*matr[0] + g*matr[3];
  adj_m[1] = f*matr[1] - g*matr[2];
  adj_m[2] = f*matr[3] - g*matr[1];
  adj_m[3] = f*matr[4] + g*matr[0];

  g_obj[1] = w[0]*adj_m[0] + w[1]*adj_m[1];
  g_obj[2] =                 w[2]*adj_m[1];
  g_obj[0] = -g_obj[1] - g_obj[2];

  g_obj[4] = w[0]*adj_m[2] + w[1]*adj_m[3];
  g_obj[5] =                 w[2]*adj_m[3];
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

int h_fcn(double *obj, double g_obj[12], double h_obj[78], const double x[12],
	  const double w[6], const double t[9])
{
  static double matr[9], f, g, t1, loc1;
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
  g = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];
  if (g <= epsilon) { *obj = g; return 1; }

  dg[3] = matr[2]*matr[7] - matr[1]*matr[8];
  dg[4] = matr[0]*matr[8] - matr[2]*matr[6];
  dg[5] = matr[1]*matr[6] - matr[0]*matr[7];
  dg[6] = matr[1]*matr[5] - matr[2]*matr[4];
  dg[7] = matr[2]*matr[3] - matr[0]*matr[5];
  dg[8] = matr[0]*matr[4] - matr[1]*matr[3];

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] + 
      matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8];

  /* Calculate objective function. */
  loc1 = a * pow(g, b);
  *obj = f * loc1;

  /* Calculate constants required */
  g = 1.0 / g;
  dobj_df = 2.0 * loc1;
  dobj_dg = b * (*obj) * g; 
  dobj_dfdg = b * dobj_df * g;
  dobj_dgdg = bm1 * dobj_dg * g;

  adj_m[0] = dobj_df*matr[0] + dobj_dg*dg[0];
  adj_m[1] = dobj_df*matr[1] + dobj_dg*dg[1];
  adj_m[2] = dobj_df*matr[2] + dobj_dg*dg[2];
  adj_m[3] = dobj_df*matr[3] + dobj_dg*dg[3];
  adj_m[4] = dobj_df*matr[4] + dobj_dg*dg[4];
  adj_m[5] = dobj_df*matr[5] + dobj_dg*dg[5];
  adj_m[6] = dobj_df*matr[6] + dobj_dg*dg[6];
  adj_m[7] = dobj_df*matr[7] + dobj_dg*dg[7];
  adj_m[8] = dobj_df*matr[8] + dobj_dg*dg[8];

  /* Calculate the derivative of the objective function. */
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
  adj_m[0] = dobj_dg*matr[0]; matr[0] *= dobj_dfdg;
  adj_m[1] = dobj_dg*matr[1]; matr[1] *= dobj_dfdg;
  adj_m[2] = dobj_dg*matr[2]; matr[2] *= dobj_dfdg;
  adj_m[3] = dobj_dg*matr[3]; matr[3] *= dobj_dfdg;
  adj_m[4] = dobj_dg*matr[4]; matr[4] *= dobj_dfdg;
  adj_m[5] = dobj_dg*matr[5]; matr[5] *= dobj_dfdg;
  adj_m[6] = dobj_dg*matr[6]; matr[6] *= dobj_dfdg;
  adj_m[7] = dobj_dg*matr[7]; matr[7] *= dobj_dfdg;
  adj_m[8] = dobj_dg*matr[8]; matr[8] *= dobj_dfdg;

  /* Blocks for the Hessian construction */
  loc1 = dobj_dgdg*dg[0] + matr[0];
  J_A[0] = dobj_df + dg[0]*(matr[0] + loc1);
  J_A[1] = dg[0]*matr[1] + loc1*dg[1];
  J_A[2] = dg[0]*matr[2] + loc1*dg[2];
  J_B[0] = dg[0]*matr[3] + loc1*dg[3];
  J_B[1] = dg[0]*matr[4] + loc1*dg[4] + adj_m[8];
  J_B[2] = dg[0]*matr[5] + loc1*dg[5] - adj_m[7];
  J_C[0] = dg[0]*matr[6] + loc1*dg[6];
  J_C[1] = dg[0]*matr[7] + loc1*dg[7] - adj_m[5];
  J_C[2] = dg[0]*matr[8] + loc1*dg[8] + adj_m[4];

  loc1 = dobj_dgdg*dg[1] + matr[1];
  J_A[3] = dobj_df + dg[1]*(matr[1] + loc1);
  J_A[4] = dg[1]*matr[2] + loc1*dg[2];
  J_B[3] = dg[1]*matr[3] + loc1*dg[3] - adj_m[8];
  J_B[4] = dg[1]*matr[4] + loc1*dg[4];
  J_B[5] = dg[1]*matr[5] + loc1*dg[5] + adj_m[6];
  J_C[3] = dg[1]*matr[6] + loc1*dg[6] + adj_m[5];
  J_C[4] = dg[1]*matr[7] + loc1*dg[7];
  J_C[5] = dg[1]*matr[8] + loc1*dg[8] - adj_m[3];

  loc1 = dobj_dgdg*dg[2] + matr[2];
  J_A[5] = dobj_df + dg[2]*(matr[2] + loc1);
  J_B[6] = dg[2]*matr[3] + loc1*dg[3] + adj_m[7];
  J_B[7] = dg[2]*matr[4] + loc1*dg[4] - adj_m[6];
  J_B[8] = dg[2]*matr[5] + loc1*dg[5];
  J_C[6] = dg[2]*matr[6] + loc1*dg[6] - adj_m[4];
  J_C[7] = dg[2]*matr[7] + loc1*dg[7] + adj_m[3];
  J_C[8] = dg[2]*matr[8] + loc1*dg[8];

  loc1 = dobj_dgdg*dg[3] + matr[3];
  J_D[0] = dobj_df + dg[3]*(matr[3] + loc1);
  J_D[1] = dg[3]*matr[4] + loc1*dg[4];
  J_D[2] = dg[3]*matr[5] + loc1*dg[5];
  J_E[0] = dg[3]*matr[6] + loc1*dg[6];
  J_E[1] = dg[3]*matr[7] + loc1*dg[7] + adj_m[2];
  J_E[2] = dg[3]*matr[8] + loc1*dg[8] - adj_m[1];

  loc1 = dobj_dgdg*dg[4] + matr[4];
  J_D[3] = dobj_df + dg[4]*(matr[4] + loc1);
  J_D[4] = dg[4]*matr[5] + loc1*dg[5];
  J_E[3] = dg[4]*matr[6] + loc1*dg[6] - adj_m[2];
  J_E[4] = dg[4]*matr[7] + loc1*dg[7];
  J_E[5] = dg[4]*matr[8] + loc1*dg[8] + adj_m[0];

  loc1 = dobj_dgdg*dg[5] + matr[5];
  J_D[5] = dobj_df + dg[5]*(matr[5] + loc1);
  J_E[6] = dg[5]*matr[6] + loc1*dg[6] + adj_m[1];
  J_E[7] = dg[5]*matr[7] + loc1*dg[7] - adj_m[0];
  J_E[8] = dg[5]*matr[8] + loc1*dg[8];

  loc1 = dobj_dgdg*dg[6] + matr[6];
  J_F[0] = dobj_df + dg[6]*(matr[6] + loc1);
  J_F[1] = dg[6]*matr[7] + loc1*dg[7];
  J_F[2] = dg[6]*matr[8] + loc1*dg[8];

  loc1 = dobj_dgdg*dg[7] + matr[7];
  J_F[3] = dobj_df + dg[7]*(matr[7] + loc1);
  J_F[4] = dg[7]*matr[8] + loc1*dg[8];

  J_F[5] = dobj_df + dg[8]*(2.0*matr[8] + dobj_dgdg*dg[8]);

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
  static double matr[9], f, g, t1, loc1;
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
  g = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] + 
      matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8];

  /* Calculate constants required */
  loc1 = a * pow(g, b);
  g = 1.0 / g;

  dobj_df = 2.0 * loc1;
  dobj_dg = b * f * loc1 * g; 
  dobj_dfdg = b * dobj_df * g;
  dobj_dgdg = bm1 * dobj_dg * g;

  /* Start of the Hessian evaluation */
  adj_m[0] = dobj_dg*matr[0]; matr[0] *= dobj_dfdg;
  adj_m[1] = dobj_dg*matr[1]; matr[1] *= dobj_dfdg;
  adj_m[2] = dobj_dg*matr[2]; matr[2] *= dobj_dfdg;
  adj_m[3] = dobj_dg*matr[3]; matr[3] *= dobj_dfdg;
  adj_m[4] = dobj_dg*matr[4]; matr[4] *= dobj_dfdg;
  adj_m[5] = dobj_dg*matr[5]; matr[5] *= dobj_dfdg;
  adj_m[6] = dobj_dg*matr[6]; matr[6] *= dobj_dfdg;
  adj_m[7] = dobj_dg*matr[7]; matr[7] *= dobj_dfdg;
  adj_m[8] = dobj_dg*matr[8]; matr[8] *= dobj_dfdg;

  /* Blocks for the Hessian construction */
  loc1 = dobj_dgdg*dg[0] + matr[0];
  J_A[0] = dobj_df + dg[0]*(matr[0] + loc1);
  J_A[1] = dg[0]*matr[1] + loc1*dg[1];
  J_A[2] = dg[0]*matr[2] + loc1*dg[2];
  J_B[0] = dg[0]*matr[3] + loc1*dg[3];
  J_B[1] = dg[0]*matr[4] + loc1*dg[4] + adj_m[8];
  J_B[2] = dg[0]*matr[5] + loc1*dg[5] - adj_m[7];
  J_C[0] = dg[0]*matr[6] + loc1*dg[6];
  J_C[1] = dg[0]*matr[7] + loc1*dg[7] - adj_m[5];
  J_C[2] = dg[0]*matr[8] + loc1*dg[8] + adj_m[4];

  loc1 = dobj_dgdg*dg[1] + matr[1];
  J_A[3] = dobj_df + dg[1]*(matr[1] + loc1);
  J_A[4] = dg[1]*matr[2] + loc1*dg[2];
  J_B[3] = dg[1]*matr[3] + loc1*dg[3] - adj_m[8];
  J_B[4] = dg[1]*matr[4] + loc1*dg[4];
  J_B[5] = dg[1]*matr[5] + loc1*dg[5] + adj_m[6];
  J_C[3] = dg[1]*matr[6] + loc1*dg[6] + adj_m[5];
  J_C[4] = dg[1]*matr[7] + loc1*dg[7];
  J_C[5] = dg[1]*matr[8] + loc1*dg[8] - adj_m[3];

  loc1 = dobj_dgdg*dg[2] + matr[2];
  J_A[5] = dobj_df + dg[2]*(matr[2] + loc1);
  J_B[6] = dg[2]*matr[3] + loc1*dg[3] + adj_m[7];
  J_B[7] = dg[2]*matr[4] + loc1*dg[4] - adj_m[6];
  J_B[8] = dg[2]*matr[5] + loc1*dg[5];
  J_C[6] = dg[2]*matr[6] + loc1*dg[6] - adj_m[4];
  J_C[7] = dg[2]*matr[7] + loc1*dg[7] + adj_m[3];
  J_C[8] = dg[2]*matr[8] + loc1*dg[8];

  loc1 = dobj_dgdg*dg[3] + matr[3];
  J_D[0] = dobj_df + dg[3]*(matr[3] + loc1);
  J_D[1] = dg[3]*matr[4] + loc1*dg[4];
  J_D[2] = dg[3]*matr[5] + loc1*dg[5];
  J_E[0] = dg[3]*matr[6] + loc1*dg[6];
  J_E[1] = dg[3]*matr[7] + loc1*dg[7] + adj_m[2];
  J_E[2] = dg[3]*matr[8] + loc1*dg[8] - adj_m[1];

  loc1 = dobj_dgdg*dg[4] + matr[4];
  J_D[3] = dobj_df + dg[4]*(matr[4] + loc1);
  J_D[4] = dg[4]*matr[5] + loc1*dg[5];
  J_E[3] = dg[4]*matr[6] + loc1*dg[6] - adj_m[2];
  J_E[4] = dg[4]*matr[7] + loc1*dg[7];
  J_E[5] = dg[4]*matr[8] + loc1*dg[8] + adj_m[0];

  loc1 = dobj_dgdg*dg[5] + matr[5];
  J_D[5] = dobj_df + dg[5]*(matr[5] + loc1);
  J_E[6] = dg[5]*matr[6] + loc1*dg[6] + adj_m[1];
  J_E[7] = dg[5]*matr[7] + loc1*dg[7] - adj_m[0];
  J_E[8] = dg[5]*matr[8] + loc1*dg[8];

  loc1 = dobj_dgdg*dg[6] + matr[6];
  J_F[0] = dobj_df + dg[6]*(matr[6] + loc1);
  J_F[1] = dg[6]*matr[7] + loc1*dg[7];
  J_F[2] = dg[6]*matr[8] + loc1*dg[8];

  loc1 = dobj_dgdg*dg[7] + matr[7];
  J_F[3] = dobj_df + dg[7]*(matr[7] + loc1);
  J_F[4] = dg[7]*matr[8] + loc1*dg[8];

  J_F[5] = dobj_df + dg[8]*(2.0*matr[8] + dobj_dgdg*dg[8]);

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

