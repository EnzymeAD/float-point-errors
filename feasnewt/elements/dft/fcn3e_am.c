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

#define tsqrt3  1.15470053837925159591885903972e+00        /*  2.0/sqrt(3.0) */
#define tsqrt6  1.22474487139158915752946794252e+00        /*  3.0/sqrt(6.0) */

#define e1a 1.00000000000	/* 1.0 / 1.0             */
#define e1b 0.57735026919	/* 1.0 / sqrt(3.0)       */
#define e1c 0.40824829046	/* 1.0 / sqrt(6.0)       */
#define e2b 0.33333333333	/* 1.0 / 3.0             */
#define e2c 0.23570226040	/* 1.0 / (3.0*sqrt(2.0)) */
#define e3c 0.16666666667	/* 1.0 / 6.0             */

int o_fcn(double *obj, const double x[12])
{
  /* 68 flops */

  static double matr[9], f, t1, t2;
  static double matd[3], g;

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
  t1 = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) +
       matr[1]*(matr[5]*matr[6] - matr[3]*matr[8]) +
       matr[2]*(matr[3]*matr[7] - matr[4]*matr[6]);
  t2 = sqrt(t1*t1 + fd2);
  g = t1 + t2;

  matd[0] = matr[0] - 1.0;
  matd[1] = matr[4] - 1.0;
  matd[2] = matr[8] - 1.0;

  /* Calculate norm(M). */
  f = matd[0]*matd[0] + matr[1]*matr[1] + matr[2]*matr[2] +
      matr[3]*matr[3] + matd[1]*matd[1] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matd[2]*matd[2];

  /* Calculate objective function. */
  (*obj) = a * f * pow(g, b);
  return 0;
}

int g_fcn(double *obj, double g_obj[12], const double x[12])
{
  /* 137 flops */

  static double matr[9], f, t1, t2;
  static double matd[3], g;
  static double adj_m[9], loc1, loc2, loc3, loc4;

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
  t1 = matr[0]*loc1 + matr[1]*loc2 + matr[2]*loc3;
  t2 = sqrt(t1*t1 + fd2);
  g = t1 + t2;

  matd[0] = matr[0] - 1.0;
  matd[1] = matr[4] - 1.0;
  matd[2] = matr[8] - 1.0;

  /* Calculate norm(M). */
  f = matd[0]*matd[0] + matr[1]*matr[1] + matr[2]*matr[2] + 
      matr[3]*matr[3] + matd[1]*matd[1] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matd[2]*matd[2];
 
  /* Calculate objective function. */
  loc4 = a * pow(g, b);
  *obj = f * loc4;

  /* Calculate the derivative of the objective function. */
  f = 2.0 * loc4;
  g = b * (*obj) / t2;

  adj_m[0] = f*matd[0] + g*loc1;
  adj_m[1] = -sqrt3*(f*matr[1] + g*loc2);
  adj_m[2] = -sqrt6*(f*matr[2] + g*loc3);

  loc1 = g*matr[0];
  loc2 = g*matr[1];
  loc3 = g*matr[2];

  adj_m[3] = f*matr[3] + loc3*matr[7] - loc2*matr[8];
  adj_m[4] = -sqrt3*(f*matd[1] + loc1*matr[8] - loc3*matr[6]);
  adj_m[5] = -sqrt6*(f*matr[5] + loc2*matr[6] - loc1*matr[7]);

  adj_m[6] = f*matr[6] + loc2*matr[5] - loc3*matr[4];
  adj_m[7] = -sqrt3*(f*matr[7] + loc3*matr[3] - loc1*matr[5]);
  adj_m[8] = -sqrt6*(f*matd[2] + loc1*matr[4] - loc2*matr[3]);

  loc1 = adj_m[1] + adj_m[2];
  g_obj[0] = loc1 - adj_m[0];
  g_obj[1] = loc1 + adj_m[0];
  g_obj[2] = adj_m[2] - 2.0*adj_m[1];
  g_obj[3] = -3.0*adj_m[2];

  loc1 = adj_m[4] + adj_m[5];
  g_obj[4] = loc1 - adj_m[3];
  g_obj[5] = loc1 + adj_m[3];
  g_obj[6] = adj_m[5] - 2.0*adj_m[4];
  g_obj[7] = -3.0*adj_m[5];

  loc1 = adj_m[7] + adj_m[8];
  g_obj[8] = loc1 - adj_m[6];
  g_obj[9] = loc1 + adj_m[6];
  g_obj[10] = adj_m[8] - 2.0*adj_m[7];
  g_obj[11] = -3.0*adj_m[8];
  return 0;
}

/*****************************************************************************/
/* The Hessian calculation is done by blocks.  Only the upper triangular     */
/* blocks are stored.  The results in the data is in the following order:    */
/*    [d1 b1 b2 d2 b3 d3 ]                                                   */
/* The matrices on the diagonal (d1,d2,d3) each contain 10 elements, while   */
/* the off-diagonal elements (b1,b2,b3) each contain 16 elements.            */
/*****************************************************************************/

static void assemble_diag(double h_obj[10], const double J_A[6])
{
  const double A0 = J_A[1] + J_A[2];
  const double A1 = J_A[3] + J_A[4];
  const double A2 = J_A[4] + J_A[5];
  static double loc1, loc2, loc3, loc4;

  loc1 = A0 + J_A[0];
  loc2 = A1 + J_A[1];
  loc3 = A2 + J_A[2];
  loc4 = loc2 + loc3;

  h_obj[0] = loc4 + loc1;
  h_obj[1] = loc4 - loc1;
  h_obj[2] = loc3 - 2.0*loc2;
  h_obj[3] = -3.0*loc3;

  loc1 = A0 - J_A[0];
  loc2 = A1 - J_A[1];
  loc3 = A2 - J_A[2];
  loc4 = loc2 + loc3;

  h_obj[4] = loc4 - loc1;
  h_obj[5] = loc3 - 2.0*loc2;
  h_obj[6] = -3.0*loc3;

  loc2 = J_A[4] - 2.0*J_A[3];
  loc3 = J_A[5] - 2.0*J_A[4];

  h_obj[7] = loc3 - 2.0*loc2;
  h_obj[8] = -3.0*loc3;

  h_obj[9] = 9.0*J_A[5];
  return;
}

static void assemble_offdiag(double h_obj[16], const double J_B[9])
{
  const double A0 = J_B[3] + J_B[6];
  const double A1 = J_B[4] + J_B[7];
  const double A2 = J_B[5] + J_B[8];
  static double loc1, loc2, loc3, loc4;

  loc1 = A0 + J_B[0];
  loc2 = A1 + J_B[1];
  loc3 = A2 + J_B[2];
  loc4 = loc2 + loc3;

  h_obj[0] = loc4 + loc1;
  h_obj[1] = loc4 - loc1;
  h_obj[2] = loc3 - 2.0*loc2;
  h_obj[3] = -3.0*loc3;

  loc1 = A0 - J_B[0];
  loc2 = A1 - J_B[1];
  loc3 = A2 - J_B[2];
  loc4 = loc2 + loc3;

  h_obj[4] = loc4 + loc1;
  h_obj[5] = loc4 - loc1;
  h_obj[6] = loc3 - 2.0*loc2;
  h_obj[7] = -3.0*loc3;

  loc1 = J_B[6] - 2.0*J_B[3];
  loc2 = J_B[7] - 2.0*J_B[4];
  loc3 = J_B[8] - 2.0*J_B[5];
  loc4 = loc2 + loc3;

  h_obj[8] = loc4 + loc1;
  h_obj[9] = loc4 - loc1;
  h_obj[10] = loc3 - 2.0*loc2;
  h_obj[11] = -3.0*loc3;

  loc1 = -3.0*J_B[6];
  loc2 = -3.0*J_B[7];
  loc3 = -3.0*J_B[8];
  loc4 = loc2 + loc3;

  h_obj[12] = loc4 + loc1;
  h_obj[13] = loc4 - loc1;
  h_obj[14] = loc3 - 2.0*loc2;
  h_obj[15] = -3.0*loc3;
  return;
}

int h_fcn(double *obj, double g_obj[12], double h_obj[78], const double x[12])
{
  /* 585 flops */

  static double matr[9], f, t1, t2;
  static double matd[3], g, t3, loc1;
  static double adj_m[9], dg[9], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;
  static double J_A[6], J_B[10], J_C[10], J_D[6], J_E[10], J_F[6];

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

  matd[0] = matr[0] - 1.0;
  matd[1] = matr[4] - 1.0;
  matd[2] = matr[8] - 1.0;

  /* Calculate norm(M). */
  f = matd[0]*matd[0] + matr[1]*matr[1] + matr[2]*matr[2] + 
      matr[3]*matr[3] + matd[1]*matd[1] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matd[2]*matd[2];

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
  adj_m[1] = -sqrt3*(dobj_df*matr[1] + dobj_dg*dg[1]);
  adj_m[2] = -sqrt6*(dobj_df*matr[2] + dobj_dg*dg[2]);
  adj_m[3] = dobj_df*matr[3] + dobj_dg*dg[3];
  adj_m[4] = -sqrt3*(dobj_df*matd[1] + dobj_dg*dg[4]);
  adj_m[5] = -sqrt6*(dobj_df*matr[5] + dobj_dg*dg[5]);
  adj_m[6] = dobj_df*matr[6] + dobj_dg*dg[6];
  adj_m[7] = -sqrt3*(dobj_df*matr[7] + dobj_dg*dg[7]);
  adj_m[8] = -sqrt6*(dobj_df*matd[2] + dobj_dg*dg[8]);

  loc1 = adj_m[1] + adj_m[2];
  g_obj[0] = loc1 - adj_m[0];
  g_obj[1] = loc1 + adj_m[0];
  g_obj[2] = adj_m[2] - 2.0*adj_m[1];
  g_obj[3] = -3.0*adj_m[2];

  loc1 = adj_m[4] + adj_m[5];
  g_obj[4] = loc1 - adj_m[3];
  g_obj[5] = loc1 + adj_m[3];
  g_obj[6] = adj_m[5] - 2.0*adj_m[4];
  g_obj[7] = -3.0*adj_m[5];

  loc1 = adj_m[7] + adj_m[8];
  g_obj[8] = loc1 - adj_m[6];
  g_obj[9] = loc1 + adj_m[6];
  g_obj[10] = adj_m[8] - 2.0*adj_m[7];
  g_obj[11] = -3.0*adj_m[8];

  /* Start of Hessian evaluation */
  adj_m[0] = dobj_dg*matr[0]; matd[0] *= dobj_dfdg;
  adj_m[1] = dobj_dg*matr[1]; matr[1] *= dobj_dfdg;
  adj_m[2] = dobj_dg*matr[2]; matr[2] *= dobj_dfdg;
  adj_m[3] = dobj_dg*matr[3]; matr[3] *= dobj_dfdg;
  adj_m[4] = dobj_dg*matr[4]; matd[1] *= dobj_dfdg;
  adj_m[5] = dobj_dg*matr[5]; matr[5] *= dobj_dfdg;
  adj_m[6] = dobj_dg*matr[6]; matr[6] *= dobj_dfdg;
  adj_m[7] = dobj_dg*matr[7]; matr[7] *= dobj_dfdg;
  adj_m[8] = dobj_dg*matr[8]; matd[2] *= dobj_dfdg;

  /* Blocks for the Hessian construction */
  loc1 = dobj_dgdg*dg[0] + matd[0];
  J_A[0] = dobj_df + dg[0]*(matd[0] + loc1);
  J_A[1] = e1b*(dg[0]*matr[1] + loc1*dg[1]);
  J_A[2] = e1c*(dg[0]*matr[2] + loc1*dg[2]);
  J_B[0] =      dg[0]*matr[3] + loc1*dg[3];
  J_B[1] = e1b*(dg[0]*matd[1] + loc1*dg[4] + adj_m[8]);
  J_B[2] = e1c*(dg[0]*matr[5] + loc1*dg[5] - adj_m[7]);
  J_C[0] =      dg[0]*matr[6] + loc1*dg[6];
  J_C[1] = e1b*(dg[0]*matr[7] + loc1*dg[7] - adj_m[5]);
  J_C[2] = e1c*(dg[0]*matd[2] + loc1*dg[8] + adj_m[4]);

  loc1 = dobj_dgdg*dg[1] + matr[1];
  J_A[3] = e2b*(dobj_df + dg[1]*(matr[1] + loc1));
  J_A[4] = e2c*(dg[1]*matr[2] + loc1*dg[2]);
  J_B[3] = e1b*(dg[1]*matr[3] + loc1*dg[3] - adj_m[8]);
  J_B[4] = e2b*(dg[1]*matd[1] + loc1*dg[4]);
  J_B[5] = e2c*(dg[1]*matr[5] + loc1*dg[5] + adj_m[6]);
  J_C[3] = e1b*(dg[1]*matr[6] + loc1*dg[6] + adj_m[5]);
  J_C[4] = e2b*(dg[1]*matr[7] + loc1*dg[7]);
  J_C[5] = e2c*(dg[1]*matd[2] + loc1*dg[8] - adj_m[3]);

  loc1 = dobj_dgdg*dg[2] + matr[2];
  J_A[5] = e3c*(dobj_df + dg[2]*(matr[2] + loc1));
  J_B[6] = e1c*(dg[2]*matr[3] + loc1*dg[3] + adj_m[7]);
  J_B[7] = e2c*(dg[2]*matd[1] + loc1*dg[4] - adj_m[6]);
  J_B[8] = e3c*(dg[2]*matr[5] + loc1*dg[5]);
  J_C[6] = e1c*(dg[2]*matr[6] + loc1*dg[6] - adj_m[4]);
  J_C[7] = e2c*(dg[2]*matr[7] + loc1*dg[7] + adj_m[3]);
  J_C[8] = e3c*(dg[2]*matd[2] + loc1*dg[8]);

  loc1 = dobj_dgdg*dg[3] + matr[3];
  J_D[0] = dobj_df + dg[3]*(matr[3] + loc1);
  J_D[1] = e1b*(dg[3]*matd[1] + loc1*dg[4]);
  J_D[2] = e1c*(dg[3]*matr[5] + loc1*dg[5]);
  J_E[0] =      dg[3]*matr[6] + loc1*dg[6];
  J_E[1] = e1b*(dg[3]*matr[7] + loc1*dg[7] + adj_m[2]);
  J_E[2] = e1c*(dg[3]*matd[2] + loc1*dg[8] - adj_m[1]);

  loc1 = dobj_dgdg*dg[4] + matd[1];
  J_D[3] = e2b*(dobj_df + dg[4]*(matd[1] + loc1));
  J_D[4] = e2c*(dg[4]*matr[5] + loc1*dg[5]);
  J_E[3] = e1b*(dg[4]*matr[6] + loc1*dg[6] - adj_m[2]);
  J_E[4] = e2b*(dg[4]*matr[7] + loc1*dg[7]);
  J_E[5] = e2c*(dg[4]*matd[2] + loc1*dg[8] + adj_m[0]);

  loc1 = dobj_dgdg*dg[5] + matr[5];
  J_D[5] = e3c*(dobj_df + dg[5]*(matr[5] + loc1));
  J_E[6] = e1c*(dg[5]*matr[6] + loc1*dg[6] + adj_m[1]);
  J_E[7] = e2c*(dg[5]*matr[7] + loc1*dg[7] - adj_m[0]);
  J_E[8] = e3c*(dg[5]*matd[2] + loc1*dg[8]);

  loc1 = dobj_dgdg*dg[6] + matr[6];
  J_F[0] = dobj_df + dg[6]*(matr[6] + loc1);
  J_F[1] = e1b*(dg[6]*matr[7] + loc1*dg[7]);
  J_F[2] = e1c*(dg[6]*matd[2] + loc1*dg[8]);

  loc1 = dobj_dgdg*dg[7] + matr[7];
  J_F[3] = e2b*(dobj_df + dg[7]*(matr[7] + loc1));
  J_F[4] = e2c*(dg[7]*matd[2] + loc1*dg[8]);

  J_F[5] = e3c*(dobj_df + dg[8]*(2.0*matd[2] + dobj_dgdg*dg[8]));

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
  /* 534 flops */

  static double matr[9], f, t1, t2;
  static double matd[3], g, t3, loc1;
  static double adj_m[9], dg[9], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;
  static double J_A[6], J_B[10], J_C[10], J_D[6], J_E[10], J_F[6];

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

  matd[0] = matr[0] - 1.0;
  matd[1] = matr[4] - 1.0;
  matd[2] = matr[8] - 1.0;

  /* Calculate norm(M). */
  f = matd[0]*matd[0] + matr[1]*matr[1] + matr[2]*matr[2] + 
      matr[3]*matr[3] + matd[1]*matd[1] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matd[2]*matd[2];

  /* Calculate objective function. */
  loc1 = a * pow(g, b);

  /* Calculate constants required */

  /* Calculate the derivative of the objective function. */
  t3 = 1.0 / t3;
  dobj_df = 2.0 * loc1;
  dobj_dg = b * f * loc1 * t3; 
  dobj_dfdg = b * dobj_df * t3;
  dobj_dgdg = dobj_dg * (bm1*t3 + fd2/(t2*g));

  /* Start of Hessian evaluation */
  adj_m[0] = dobj_dg*matr[0]; matd[0] *= dobj_dfdg;
  adj_m[1] = dobj_dg*matr[1]; matr[1] *= dobj_dfdg;
  adj_m[2] = dobj_dg*matr[2]; matr[2] *= dobj_dfdg;
  adj_m[3] = dobj_dg*matr[3]; matr[3] *= dobj_dfdg;
  adj_m[4] = dobj_dg*matr[4]; matd[1] *= dobj_dfdg;
  adj_m[5] = dobj_dg*matr[5]; matr[5] *= dobj_dfdg;
  adj_m[6] = dobj_dg*matr[6]; matr[6] *= dobj_dfdg;
  adj_m[7] = dobj_dg*matr[7]; matr[7] *= dobj_dfdg;
  adj_m[8] = dobj_dg*matr[8]; matd[2] *= dobj_dfdg;

  /* Blocks for the Hessian construction */
  loc1 = dobj_dgdg*dg[0] + matd[0];
  J_A[0] = dobj_df + dg[0]*(matd[0] + loc1);
  J_A[1] = e1b*(dg[0]*matr[1] + loc1*dg[1]);
  J_A[2] = e1c*(dg[0]*matr[2] + loc1*dg[2]);
  J_B[0] =      dg[0]*matr[3] + loc1*dg[3];
  J_B[1] = e1b*(dg[0]*matd[1] + loc1*dg[4] + adj_m[8]);
  J_B[2] = e1c*(dg[0]*matr[5] + loc1*dg[5] - adj_m[7]);
  J_C[0] =      dg[0]*matr[6] + loc1*dg[6];
  J_C[1] = e1b*(dg[0]*matr[7] + loc1*dg[7] - adj_m[5]);
  J_C[2] = e1c*(dg[0]*matd[2] + loc1*dg[8] + adj_m[4]);

  loc1 = dobj_dgdg*dg[1] + matr[1];
  J_A[3] = e2b*(dobj_df + dg[1]*(matr[1] + loc1));
  J_A[4] = e2c*(dg[1]*matr[2] + loc1*dg[2]);
  J_B[3] = e1b*(dg[1]*matr[3] + loc1*dg[3] - adj_m[8]);
  J_B[4] = e2b*(dg[1]*matd[1] + loc1*dg[4]);
  J_B[5] = e2c*(dg[1]*matr[5] + loc1*dg[5] + adj_m[6]);
  J_C[3] = e1b*(dg[1]*matr[6] + loc1*dg[6] + adj_m[5]);
  J_C[4] = e2b*(dg[1]*matr[7] + loc1*dg[7]);
  J_C[5] = e2c*(dg[1]*matd[2] + loc1*dg[8] - adj_m[3]);

  loc1 = dobj_dgdg*dg[2] + matr[2];
  J_A[5] = e3c*(dobj_df + dg[2]*(matr[2] + loc1));
  J_B[6] = e1c*(dg[2]*matr[3] + loc1*dg[3] + adj_m[7]);
  J_B[7] = e2c*(dg[2]*matd[1] + loc1*dg[4] - adj_m[6]);
  J_B[8] = e3c*(dg[2]*matr[5] + loc1*dg[5]);
  J_C[6] = e1c*(dg[2]*matr[6] + loc1*dg[6] - adj_m[4]);
  J_C[7] = e2c*(dg[2]*matr[7] + loc1*dg[7] + adj_m[3]);
  J_C[8] = e3c*(dg[2]*matd[2] + loc1*dg[8]);

  loc1 = dobj_dgdg*dg[3] + matr[3];
  J_D[0] = dobj_df + dg[3]*(matr[3] + loc1);
  J_D[1] = e1b*(dg[3]*matd[1] + loc1*dg[4]);
  J_D[2] = e1c*(dg[3]*matr[5] + loc1*dg[5]);
  J_E[0] =      dg[3]*matr[6] + loc1*dg[6];
  J_E[1] = e1b*(dg[3]*matr[7] + loc1*dg[7] + adj_m[2]);
  J_E[2] = e1c*(dg[3]*matd[2] + loc1*dg[8] - adj_m[1]);

  loc1 = dobj_dgdg*dg[4] + matd[1];
  J_D[3] = e2b*(dobj_df + dg[4]*(matd[1] + loc1));
  J_D[4] = e2c*(dg[4]*matr[5] + loc1*dg[5]);
  J_E[3] = e1b*(dg[4]*matr[6] + loc1*dg[6] - adj_m[2]);
  J_E[4] = e2b*(dg[4]*matr[7] + loc1*dg[7]);
  J_E[5] = e2c*(dg[4]*matd[2] + loc1*dg[8] + adj_m[0]);

  loc1 = dobj_dgdg*dg[5] + matr[5];
  J_D[5] = e3c*(dobj_df + dg[5]*(matr[5] + loc1));
  J_E[6] = e1c*(dg[5]*matr[6] + loc1*dg[6] + adj_m[1]);
  J_E[7] = e2c*(dg[5]*matr[7] + loc1*dg[7] - adj_m[0]);
  J_E[8] = e3c*(dg[5]*matd[2] + loc1*dg[8]);

  loc1 = dobj_dgdg*dg[6] + matr[6];
  J_F[0] = dobj_df + dg[6]*(matr[6] + loc1);
  J_F[1] = e1b*(dg[6]*matr[7] + loc1*dg[7]);
  J_F[2] = e1c*(dg[6]*matd[2] + loc1*dg[8]);

  loc1 = dobj_dgdg*dg[7] + matr[7];
  J_F[3] = e2b*(dobj_df + dg[7]*(matr[7] + loc1));
  J_F[4] = e2c*(dg[7]*matd[2] + loc1*dg[8]);

  J_F[5] = e3c*(dobj_df + dg[8]*(2.0*matd[2] + dobj_dgdg*dg[8]));

  assemble_diag(h_obj, J_A);
  assemble_diag(h_obj + 42, J_D);
  assemble_diag(h_obj + 68, J_F);

  assemble_offdiag(h_obj + 10, J_B);
  assemble_offdiag(h_obj + 26, J_C);
  assemble_offdiag(h_obj + 52, J_E);
  return;
}

#undef a
#undef b
#undef bm1

#undef d
#undef fd2 

#undef sqrt3
#undef sqrt6

#undef tsqrt3
#undef tsqrt6

#undef e1a
#undef e1b
#undef e1c
#undef e2b
#undef e2c
#undef e3c

