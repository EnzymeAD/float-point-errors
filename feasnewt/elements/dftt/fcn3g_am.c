#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "fcn.h"

#define a       5.00000000000000000000000000000e-01        /*  1.0/2.0       */
#define b      -1.33333333333333333333333333333e-00        /* -4.0/3.0       */
#define bm1    -2.33333333333333333333333333333e-00        /* -7.0/3.0       */

#define d	1.0e-4					   /* delta          */
#define fd2     4.0e-8                                     /*  4.0*delta^2   */

#define sqrt3   5.77350269189625797959429519858e-01        /*  1.0/sqrt(3.0) */
#define sqrt6   4.08248290463863052509822647505e-01        /*  1.0/sqrt(6.0) */

/* Assume that w is upper triangular. The t variable is the multiplication   */
/* of the rotation matrices required to obtain this matrix.                  */

const double w[9] = {1, -sqrt3, -sqrt6, 0, 2*sqrt3, -sqrt6, 0, 0, 3*sqrt6};
const double t[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};

int o_fcn(double *obj, const double x[12])
{
  /* 104 operations */

  static double matr[9], f, t1, t2;
  static double fmat[6], g;

  /* Calculate M = A*inv(W). */
  f       = x[1] - x[0];
  g       = x[2] - x[0];
  t1      = x[3] - x[0];
  matr[0] = f*w[0];
  matr[1] = f*w[1] + g*w[4];
  matr[2] = f*w[2] + g*w[5] + t1*w[8];

  f       = x[5] - x[4];
  g       = x[6] - x[4];
  t1      = x[7] - x[4];
  matr[3] = f*w[0];
  matr[4] = f*w[1] + g*w[4];
  matr[5] = f*w[2] + g*w[5] + t1*w[8];

  f       = x[9] - x[8];
  g       = x[10] - x[8];
  t1      = x[11] - x[8];
  matr[6] = f*w[0];
  matr[7] = f*w[1] + g*w[4];
  matr[8] = f*w[2] + g*w[5] + t1*w[8];

  /* Calculate det(M). */
  t1 = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) +
       matr[1]*(matr[5]*matr[6] - matr[3]*matr[8]) +
       matr[2]*(matr[3]*matr[7] - matr[4]*matr[6]);
  t2 = sqrt(t1*t1 + fd2);
  g = t1 + t2;

  /* Calculate norm(M). */
  fmat[0] = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] - 1.0;
  fmat[1] = matr[0]*matr[3] + matr[1]*matr[4] + matr[2]*matr[5];
  fmat[2] = matr[0]*matr[6] + matr[1]*matr[7] + matr[2]*matr[8];

  fmat[3] = matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] - 1.0;
  fmat[4] = matr[3]*matr[6] + matr[4]*matr[7] + matr[5]*matr[8];

  fmat[5] = matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8] - 1.0;

  f = fmat[0]*fmat[0] + 2.0*fmat[1]*fmat[1] + 2.0*fmat[2]*fmat[2] +
                            fmat[3]*fmat[3] + 2.0*fmat[4]*fmat[4] +
                                                  fmat[5]*fmat[5];

  /* Calculate objective function. */
  (*obj) = a * f * pow(g, b);
  return 0;
}

/*****************************************************************************/
/* Derivative calculation.                                                   */
/*****************************************************************************/

int g_fcn(double *obj, double g_obj[12], const double x[12])
{
  /* 227 operations (104 for function) */

  static double matr[9], f, t1, t2;
  static double fmat[6], g;
  static double adj_m[9], df[9], loc1, loc2, loc3, loc4;

  /* Calculate M = A*inv(W). */
  f       = x[1] - x[0];
  g       = x[2] - x[0];
  t1      = x[3] - x[0];
  matr[0] = f*w[0];
  matr[1] = f*w[1] + g*w[4];
  matr[2] = f*w[2] + g*w[5] + t1*w[8];

  f       = x[5] - x[4];
  g       = x[6] - x[4];
  t1      = x[7] - x[4];
  matr[3] = f*w[0];
  matr[4] = f*w[1] + g*w[4];
  matr[5] = f*w[2] + g*w[5] + t1*w[8];

  f       = x[9] - x[8];
  g       = x[10] - x[8];
  t1      = x[11] - x[8];
  matr[6] = f*w[0];
  matr[7] = f*w[1] + g*w[4];
  matr[8] = f*w[2] + g*w[5] + t1*w[8];

  /* Calculate det(M). */
  loc1 = matr[4]*matr[8] - matr[5]*matr[7];
  loc2 = matr[5]*matr[6] - matr[3]*matr[8];
  loc3 = matr[3]*matr[7] - matr[4]*matr[6];
  t1 = matr[0]*loc1 + matr[1]*loc2 + matr[2]*loc3;
  t2 = sqrt(t1*t1 + fd2);
  g = t1 + t2;

  /* Calculate norm(M). */
  fmat[0] = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] - 1.0;
  fmat[1] = matr[0]*matr[3] + matr[1]*matr[4] + matr[2]*matr[5];
  fmat[2] = matr[0]*matr[6] + matr[1]*matr[7] + matr[2]*matr[8];

  fmat[3] = matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] - 1.0;
  fmat[4] = matr[3]*matr[6] + matr[4]*matr[7] + matr[5]*matr[8];

  fmat[5] = matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8] - 1.0;

  f = fmat[0]*fmat[0] + 2.0*fmat[1]*fmat[1] + 2.0*fmat[2]*fmat[2] +
                            fmat[3]*fmat[3] + 2.0*fmat[4]*fmat[4] +
                                                  fmat[5]*fmat[5];

  /* Calculate objective function. */
  loc4 = a * pow(g, b);
  *obj = f * loc4;

  /* Calculate the derivative of the objective function. */
  f = 4.0 * loc4;                        /* Constant on nabla f */
  g = b * (*obj) / t2;                   /* Constant on nabla g */

  /* Compute d fmat by d mat */
  df[0] = fmat[0]*matr[0] + fmat[1]*matr[3] + fmat[2]*matr[6];
  df[1] = fmat[0]*matr[1] + fmat[1]*matr[4] + fmat[2]*matr[7];
  df[2] = fmat[0]*matr[2] + fmat[1]*matr[5] + fmat[2]*matr[8];

  df[3] = fmat[1]*matr[0] + fmat[3]*matr[3] + fmat[4]*matr[6];
  df[4] = fmat[1]*matr[1] + fmat[3]*matr[4] + fmat[4]*matr[7];
  df[5] = fmat[1]*matr[2] + fmat[3]*matr[5] + fmat[4]*matr[8];

  df[6] = fmat[2]*matr[0] + fmat[4]*matr[3] + fmat[5]*matr[6];
  df[7] = fmat[2]*matr[1] + fmat[4]*matr[4] + fmat[5]*matr[7];
  df[8] = fmat[2]*matr[2] + fmat[4]*matr[5] + fmat[5]*matr[8];

  adj_m[0] = f*df[0] + g*loc1;
  adj_m[1] = f*df[1] + g*loc2;
  adj_m[2] = f*df[2] + g*loc3;

  loc1 = g*matr[0];
  loc2 = g*matr[1];
  loc3 = g*matr[2];

  adj_m[3] = f*df[3] + loc3*matr[7] - loc2*matr[8];
  adj_m[4] = f*df[4] + loc1*matr[8] - loc3*matr[6];
  adj_m[5] = f*df[5] + loc2*matr[6] - loc1*matr[7];

  adj_m[6] = f*df[6] + loc2*matr[5] - loc3*matr[4];
  adj_m[7] = f*df[7] + loc3*matr[3] - loc1*matr[5];
  adj_m[8] = f*df[8] + loc1*matr[4] - loc2*matr[3];

  g_obj[1] = w[0]*adj_m[0] + w[1]*adj_m[1] + w[2]*adj_m[2];
  g_obj[2] =                 w[4]*adj_m[1] + w[5]*adj_m[2];
  g_obj[3] =                                 w[8]*adj_m[2];
  g_obj[0] = -g_obj[1] - g_obj[2] - g_obj[3];

  g_obj[5] = w[0]*adj_m[3] + w[1]*adj_m[4] + w[2]*adj_m[5];
  g_obj[6] =                 w[4]*adj_m[4] + w[5]*adj_m[5];
  g_obj[7] =                                 w[8]*adj_m[5];
  g_obj[4] = -g_obj[5] - g_obj[6] - g_obj[7];

  g_obj[9]  = w[0]*adj_m[6] + w[1]*adj_m[7] + w[2]*adj_m[8];
  g_obj[10] =                 w[4]*adj_m[7] + w[5]*adj_m[8];
  g_obj[11] =                                 w[8]*adj_m[8];
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

static void assemble_diag(double h_obj[10], const double J_A[6])
{ 
  static double A[12];
  
  A[1]  =  J_A[0]*w[0] + J_A[1]*w[1] + J_A[2]*w[2];
  A[2]  =                J_A[1]*w[4] + J_A[2]*w[5];
  A[3]  =                              J_A[2]*w[8];
  A[0]  = -A[1] - A[2] - A[3];

  A[5]  =  J_A[1]*w[0] + J_A[3]*w[1] + J_A[4]*w[2];
  A[6]  =                J_A[3]*w[4] + J_A[4]*w[5];
  A[7]  =                              J_A[4]*w[8];
  A[4]  = -A[5] - A[6] - A[7];

  A[9]  =  J_A[2]*w[0] + J_A[4]*w[1] + J_A[5]*w[2];
  A[10] =                J_A[4]*w[4] + J_A[5]*w[5];
  A[11] =                              J_A[5]*w[8];
  A[8]  = -A[9] - A[10] - A[11];

  h_obj[1] =  A[0]*w[0] + A[4]*w[1] + A[8]*w[2];
  h_obj[2] =              A[4]*w[4] + A[8]*w[5];
  h_obj[3] =                          A[8]*w[8];
  h_obj[0] = -h_obj[1] - h_obj[2] - h_obj[3];

  h_obj[4] =  A[1]*w[0] + A[5]*w[1] + A[9]*w[2];
  h_obj[5] =              A[5]*w[4] + A[9]*w[5];
  h_obj[6] =                          A[9]*w[8];

  h_obj[7] =              A[6]*w[4] + A[10]*w[5];
  h_obj[8] =                          A[10]*w[8];

  h_obj[9] =                          A[11]*w[8];
  return;
}

static void assemble_offdiag(double h_obj[16], const double J_B[9])
{
  static double A[12];

  A[1]  =  J_B[0]*w[0] + J_B[1]*w[1] + J_B[2]*w[2];
  A[2]  =                J_B[1]*w[4] + J_B[2]*w[5];
  A[3]  =                              J_B[2]*w[8];
  A[0]  = -A[1] - A[2] - A[3];

  A[5]  =  J_B[3]*w[0] + J_B[4]*w[1] + J_B[5]*w[2];
  A[6]  =                J_B[4]*w[4] + J_B[5]*w[5];
  A[7]  =                              J_B[5]*w[8];
  A[4]  = -A[5] - A[6] - A[7];

  A[9]  =  J_B[6]*w[0] + J_B[7]*w[1] + J_B[8]*w[2];
  A[10] =                J_B[7]*w[4] + J_B[8]*w[5];
  A[11] =                              J_B[8]*w[8];
  A[8]  = -A[9] - A[10] - A[11];

  h_obj[4] =  A[0]*w[0] + A[4]*w[1] + A[8]*w[2];
  h_obj[5] =  A[1]*w[0] + A[5]*w[1] + A[9]*w[2];
  h_obj[6] =  A[2]*w[0] + A[6]*w[1] + A[10]*w[2];
  h_obj[7] =  A[3]*w[0] + A[7]*w[1] + A[11]*w[2];

  h_obj[8]  = A[4]*w[4] + A[8]*w[5];
  h_obj[9]  = A[5]*w[4] + A[9]*w[5];
  h_obj[10] = A[6]*w[4] + A[10]*w[5];
  h_obj[11] = A[7]*w[4] + A[11]*w[5];

  h_obj[12] = A[8]*w[8];
  h_obj[13] = A[9]*w[8];
  h_obj[14] = A[10]*w[8];
  h_obj[15] = A[11]*w[8];

  h_obj[0] = -h_obj[4] - h_obj[8] - h_obj[12];
  h_obj[1] = -h_obj[5] - h_obj[9] - h_obj[13];
  h_obj[2] = -h_obj[6] - h_obj[10] - h_obj[14];
  h_obj[3] = -h_obj[7] - h_obj[11] - h_obj[15];
  return;
}

int h_fcn(double *obj, double g_obj[12], double h_obj[78], const double x[12])
{
  /* 988 operations (227 for gradient, 104 for function) */

  static double matr[9], f, t1, t2;
  static double fmat[6], ftmat[6], g, t3;
  static double adj_m[9], df[9], dg[9], loc1, loc2;
  static double dobj_df, dobj_dg, dobj_dfdg, dobj_dgdg;
  static double J_A[6], J_B[9], J_C[9], J_D[6], J_E[9], J_F[6];
  static double aux[45];

  /* Calculate M = A*inv(W). */
  f       = x[1] - x[0];
  g       = x[2] - x[0];
  t1      = x[3] - x[0];
  matr[0] = f*w[0];
  matr[1] = f*w[1] + g*w[4];
  matr[2] = f*w[2] + g*w[5] + t1*w[8];

  f       = x[5] - x[4];
  g       = x[6] - x[4];
  t1      = x[7] - x[4];
  matr[3] = f*w[0];
  matr[4] = f*w[1] + g*w[4];
  matr[5] = f*w[2] + g*w[5] + t1*w[8];

  f       = x[9] - x[8];
  g       = x[10] - x[8];
  t1      = x[11] - x[8];
  matr[6] = f*w[0];
  matr[7] = f*w[1] + g*w[4];
  matr[8] = f*w[2] + g*w[5] + t1*w[8];

  /* Calculate products for M*M' */
  aux[ 0] = matr[0]*matr[0];
  aux[ 1] = matr[0]*matr[1];
  aux[ 2] = matr[0]*matr[2];
  aux[ 3] = matr[0]*matr[3];
  aux[ 4] = matr[0]*matr[4];
  aux[ 5] = matr[0]*matr[5];
  aux[ 6] = matr[0]*matr[6];
  aux[ 7] = matr[0]*matr[7];
  aux[ 8] = matr[0]*matr[8];
  aux[ 9] = matr[1]*matr[1];
  aux[10] = matr[1]*matr[2];
  aux[11] = matr[1]*matr[3];
  aux[12] = matr[1]*matr[4];
  aux[13] = matr[1]*matr[5];
  aux[14] = matr[1]*matr[6];
  aux[15] = matr[1]*matr[7];
  aux[16] = matr[1]*matr[8];
  aux[17] = matr[2]*matr[2];
  aux[18] = matr[2]*matr[3];
  aux[19] = matr[2]*matr[4];
  aux[20] = matr[2]*matr[5];
  aux[21] = matr[2]*matr[6];
  aux[22] = matr[2]*matr[7];
  aux[23] = matr[2]*matr[8];
  aux[24] = matr[3]*matr[3];
  aux[25] = matr[3]*matr[4];
  aux[26] = matr[3]*matr[5];
  aux[27] = matr[3]*matr[6];
  aux[28] = matr[3]*matr[7];
  aux[29] = matr[3]*matr[8];
  aux[30] = matr[4]*matr[4];
  aux[31] = matr[4]*matr[5];
  aux[32] = matr[4]*matr[6];
  aux[33] = matr[4]*matr[7];
  aux[34] = matr[4]*matr[8];
  aux[35] = matr[5]*matr[5];
  aux[36] = matr[5]*matr[6];
  aux[37] = matr[5]*matr[7];
  aux[38] = matr[5]*matr[8];
  aux[39] = matr[6]*matr[6];
  aux[40] = matr[6]*matr[7];
  aux[41] = matr[6]*matr[8];
  aux[42] = matr[7]*matr[7];
  aux[43] = matr[7]*matr[8];
  aux[44] = matr[8]*matr[8];

  /* Calculate det(M). */
  dg[0] = aux[34] - aux[37];
  dg[1] = aux[36] - aux[29];
  dg[2] = aux[28] - aux[32];
  dg[3] = aux[22] - aux[16];
  dg[4] = aux[ 8] - aux[21];
  dg[5] = aux[14] - aux[ 7];
  dg[6] = aux[13] - aux[19];
  dg[7] = aux[18] - aux[ 5];
  dg[8] = aux[ 4] - aux[11];

  t1 = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];
  t2 = t1*t1 + fd2;
  t3 = sqrt(t2);
  g = t1 + t3;

  fmat[0] = aux[ 0] + aux[ 9] + aux[17] - 1.0;
  fmat[1] = aux[ 3] + aux[12] + aux[20];
  fmat[2] = aux[ 6] + aux[15] + aux[23];

  fmat[3] = aux[24] + aux[30] + aux[35] - 1.0;
  fmat[4] = aux[27] + aux[33] + aux[38];

  fmat[5] = aux[39] + aux[42] + aux[44] - 1.0;

  f = fmat[0]*fmat[0] + 2.0*fmat[1]*fmat[1] + 2.0*fmat[2]*fmat[2] +
                            fmat[3]*fmat[3] + 2.0*fmat[4]*fmat[4] +
                                                  fmat[5]*fmat[5];

  /* Calculate objective function. */
  loc1 = a * pow(g, b);
  *obj = f * loc1;

  /* Calculate constants required */

  /* Calculate the derivative of the objective function. */
  t3 = 1.0 / t3;
  dobj_df = 4.0 * loc1;
  dobj_dg = b * (*obj) * t3;
  dobj_dfdg = b * dobj_df * t3;
  dobj_dgdg = dobj_dg * (bm1*t3 + fd2/(t2*g));

  df[0] = fmat[0]*matr[0] + fmat[1]*matr[3] + fmat[2]*matr[6];
  df[1] = fmat[0]*matr[1] + fmat[1]*matr[4] + fmat[2]*matr[7];
  df[2] = fmat[0]*matr[2] + fmat[1]*matr[5] + fmat[2]*matr[8];

  df[3] = fmat[1]*matr[0] + fmat[3]*matr[3] + fmat[4]*matr[6];
  df[4] = fmat[1]*matr[1] + fmat[3]*matr[4] + fmat[4]*matr[7];
  df[5] = fmat[1]*matr[2] + fmat[3]*matr[5] + fmat[4]*matr[8];

  df[6] = fmat[2]*matr[0] + fmat[4]*matr[3] + fmat[5]*matr[6];
  df[7] = fmat[2]*matr[1] + fmat[4]*matr[4] + fmat[5]*matr[7];
  df[8] = fmat[2]*matr[2] + fmat[4]*matr[5] + fmat[5]*matr[8];

  adj_m[0] = dobj_df*df[0] + dobj_dg*dg[0];
  adj_m[1] = dobj_df*df[1] + dobj_dg*dg[1];
  adj_m[2] = dobj_df*df[2] + dobj_dg*dg[2];
  adj_m[3] = dobj_df*df[3] + dobj_dg*dg[3];
  adj_m[4] = dobj_df*df[4] + dobj_dg*dg[4];
  adj_m[5] = dobj_df*df[5] + dobj_dg*dg[5];
  adj_m[6] = dobj_df*df[6] + dobj_dg*dg[6];
  adj_m[7] = dobj_df*df[7] + dobj_dg*dg[7];
  adj_m[8] = dobj_df*df[8] + dobj_dg*dg[8];

  g_obj[1] = w[0]*adj_m[0] + w[1]*adj_m[1] + w[2]*adj_m[2];
  g_obj[2] =                 w[4]*adj_m[1] + w[5]*adj_m[2];
  g_obj[3] =                                 w[8]*adj_m[2];
  g_obj[0] = -g_obj[1] - g_obj[2] - g_obj[3];

  g_obj[5] = w[0]*adj_m[3] + w[1]*adj_m[4] + w[2]*adj_m[5];
  g_obj[6] =                 w[4]*adj_m[4] + w[5]*adj_m[5];
  g_obj[7] =                                 w[8]*adj_m[5];
  g_obj[4] = -g_obj[5] - g_obj[6] - g_obj[7];

  g_obj[9]  = w[0]*adj_m[6] + w[1]*adj_m[7] + w[2]*adj_m[8];
  g_obj[10] =                 w[4]*adj_m[7] + w[5]*adj_m[8];
  g_obj[11] =                                 w[8]*adj_m[8];
  g_obj[8]  = -g_obj[9] - g_obj[10] - g_obj[11];

  /* Start of the Hessian evaluation */
  ftmat[0] = aux[ 0] + aux[24] + aux[39];
  ftmat[1] = aux[ 1] + aux[25] + aux[40];
  ftmat[2] = aux[ 2] + aux[26] + aux[41];

  ftmat[3] = aux[ 9] + aux[30] + aux[42];
  ftmat[4] = aux[10] + aux[31] + aux[43];

  ftmat[5] = aux[17] + aux[35] + aux[44];

  adj_m[0] = dobj_dg*matr[0];
  adj_m[1] = dobj_dg*matr[1];
  adj_m[2] = dobj_dg*matr[2];
  adj_m[3] = dobj_dg*matr[3];
  adj_m[4] = dobj_dg*matr[4];
  adj_m[5] = dobj_dg*matr[5];
  adj_m[6] = dobj_dg*matr[6];
  adj_m[7] = dobj_dg*matr[7];
  adj_m[8] = dobj_dg*matr[8];

  /* Blocks for the Hessian construction */
  loc1 = dg[0]*dobj_dfdg;
  loc2 = dg[0]*dobj_dgdg + df[0]*dobj_dfdg;
  J_A[0] = loc1*df[0] + loc2*dg[0] + dobj_df*(fmat[0] + ftmat[0] + aux[ 0]);
  J_A[1] = loc1*df[1] + loc2*dg[1] + dobj_df*(          ftmat[1] + aux[ 1]);
  J_A[2] = loc1*df[2] + loc2*dg[2] + dobj_df*(          ftmat[2] + aux[ 2]);
  J_B[0] = loc1*df[3] + loc2*dg[3] + dobj_df*(fmat[1] + aux[3]);
  J_B[1] = loc1*df[4] + loc2*dg[4] + dobj_df*aux[11] + adj_m[8];
  J_B[2] = loc1*df[5] + loc2*dg[5] + dobj_df*aux[18] - adj_m[7];
  J_C[0] = loc1*df[6] + loc2*dg[6] + dobj_df*(fmat[2] + aux[6]);
  J_C[1] = loc1*df[7] + loc2*dg[7] + dobj_df*aux[14] - adj_m[5];
  J_C[2] = loc1*df[8] + loc2*dg[8] + dobj_df*aux[21] + adj_m[4];

  loc1 = dg[1]*dobj_dfdg;
  loc2 = dg[1]*dobj_dgdg + df[1]*dobj_dfdg;
  J_A[3] = loc1*df[1] + loc2*dg[1] + dobj_df*(fmat[0] + ftmat[3] + aux[ 9]);
  J_A[4] = loc1*df[2] + loc2*dg[2] + dobj_df*(          ftmat[4] + aux[10]);
  J_B[3] = loc1*df[3] + loc2*dg[3] + dobj_df*aux[ 4] - adj_m[8];
  J_B[4] = loc1*df[4] + loc2*dg[4] + dobj_df*(fmat[1] + aux[12]);
  J_B[5] = loc1*df[5] + loc2*dg[5] + dobj_df*aux[19] + adj_m[6];
  J_C[3] = loc1*df[6] + loc2*dg[6] + dobj_df*aux[ 7] + adj_m[5];
  J_C[4] = loc1*df[7] + loc2*dg[7] + dobj_df*(fmat[2] + aux[15]);
  J_C[5] = loc1*df[8] + loc2*dg[8] + dobj_df*aux[22] - adj_m[3];

  loc1 = dg[2]*dobj_dfdg;
  loc2 = dg[2]*dobj_dgdg + df[2]*dobj_dfdg;
  J_A[5] = loc1*df[2] + loc2*dg[2] + dobj_df*(fmat[0] + ftmat[5] + aux[17]);
  J_B[6] = loc1*df[3] + loc2*dg[3] + dobj_df*aux[ 5] + adj_m[7];
  J_B[7] = loc1*df[4] + loc2*dg[4] + dobj_df*aux[13] - adj_m[6];
  J_B[8] = loc1*df[5] + loc2*dg[5] + dobj_df*(fmat[1] + aux[20]);
  J_C[6] = loc1*df[6] + loc2*dg[6] + dobj_df*aux[ 8] - adj_m[4];
  J_C[7] = loc1*df[7] + loc2*dg[7] + dobj_df*aux[16] + adj_m[3];
  J_C[8] = loc1*df[8] + loc2*dg[8] + dobj_df*(fmat[2] + aux[23]);

  loc1 = dg[3]*dobj_dfdg;
  loc2 = dg[3]*dobj_dgdg + df[3]*dobj_dfdg;
  J_D[0] = loc1*df[3] + loc2*dg[3] + dobj_df*(fmat[3] + ftmat[0] + aux[24]);
  J_D[1] = loc1*df[4] + loc2*dg[4] + dobj_df*(          ftmat[1] + aux[25]);
  J_D[2] = loc1*df[5] + loc2*dg[5] + dobj_df*(          ftmat[2] + aux[26]);
  J_E[0] = loc1*df[6] + loc2*dg[6] + dobj_df*(fmat[4] + aux[27]);
  J_E[1] = loc1*df[7] + loc2*dg[7] + dobj_df*aux[32] + adj_m[2];
  J_E[2] = loc1*df[8] + loc2*dg[8] + dobj_df*aux[36] - adj_m[1];

  loc1 = dg[4]*dobj_dfdg; 
  loc2 = dg[4]*dobj_dgdg + df[4]*dobj_dfdg;
  J_D[3] = loc1*df[4] + loc2*dg[4] + dobj_df*(fmat[3] + ftmat[3] + aux[30]);  
  J_D[4] = loc1*df[5] + loc2*dg[5] + dobj_df*(          ftmat[4] + aux[31]);
  J_E[3] = loc1*df[6] + loc2*dg[6] + dobj_df*aux[28] - adj_m[2];
  J_E[4] = loc1*df[7] + loc2*dg[7] + dobj_df*(fmat[4] + aux[33]);
  J_E[5] = loc1*df[8] + loc2*dg[8] + dobj_df*aux[37] + adj_m[0];

  loc1 = dg[5]*dobj_dfdg;
  loc2 = dg[5]*dobj_dgdg + df[5]*dobj_dfdg;
  J_D[5] = loc1*df[5] + loc2*dg[5] + dobj_df*(fmat[3] + ftmat[5] + aux[35]);
  J_E[6] = loc1*df[6] + loc2*dg[6] + dobj_df*aux[29] + adj_m[1];
  J_E[7] = loc1*df[7] + loc2*dg[7] + dobj_df*aux[34] - adj_m[0];
  J_E[8] = loc1*df[8] + loc2*dg[8] + dobj_df*(fmat[4] + aux[38]);

  loc1 = dg[6]*dobj_dfdg;
  loc2 = dg[6]*dobj_dgdg + df[6]*dobj_dfdg;
  J_F[0] = loc1*df[6] + loc2*dg[6] + dobj_df*(fmat[5] + ftmat[0] + aux[39]);
  J_F[1] = loc1*df[7] + loc2*dg[7] + dobj_df*(          ftmat[1] + aux[40]);
  J_F[2] = loc1*df[8] + loc2*dg[8] + dobj_df*(          ftmat[2] + aux[41]);

  loc1 = dg[7]*dobj_dfdg;
  loc2 = dg[7]*dobj_dgdg + df[7]*dobj_dfdg;
  J_F[3] = loc1*df[7] + loc2*dg[7] + dobj_df*(fmat[5] + ftmat[3] + aux[42]);
  J_F[4] = loc1*df[8] + loc2*dg[8] + dobj_df*(          ftmat[4] + aux[43]);

  loc1 = dg[8]*dobj_dfdg;
  loc2 = dg[8]*dobj_dgdg + df[8]*dobj_dfdg;
  J_F[5] = loc1*df[8] + loc2*dg[8] + dobj_df*(fmat[5] + ftmat[5] + aux[44]);

  /* Put everything together */
  assemble_diag(h_obj, J_A);
  assemble_diag(h_obj + 42, J_D);
  assemble_diag(h_obj + 68, J_F);

  assemble_offdiag(h_obj + 10, J_B);
  assemble_offdiag(h_obj + 26, J_C);
  assemble_offdiag(h_obj + 52, J_E);
  return 0;
}

void h_only(double h_obj[78], const double x[12])
{
  static double obj, g_obj[12];

  h_fcn(&obj, g_obj, h_obj, x);
  return;
}

#if 0
int main()
{
  double obj, g[12], h[78], x[12];
  int i;

  srand48(1003);

  for (i = 0; i < 12; ++i) {
    x[i] = 5*(drand48() - 0.5);
  }

  h_fcn(&obj, g, h, x);

  printf("o     = %10.9e\n", obj);

  for (i = 0; i < 12; ++i) {
    printf("g(%2d) = %10.9e\n", i, g[i]);
  }

  for (i = 0; i < 78; ++i) {
    printf("h(%2d) = %10.9e\n", i, h[i]);
  }
  return -1;
}
#endif

