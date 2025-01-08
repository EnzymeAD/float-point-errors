#include <math.h>
#include "fcn.h"

#define a       5.00000000000000000000000000000e-01        /*  1.0/2.0       */
#define b      -1.00000000000000000000000000000e-00        /* -1.0/1.0       */
#define bm1    -2.00000000000000000000000000000e-00        /* -2.0/1.0       */

#define d       1.0e-4					   /*  delta         */
#define fd2     4.0e-8					   /*  4.0*delta^2   */

#define sqrt3   5.77350269189625797959429519858e-01        /*  1.0/sqrt(3.0) */

#define tsqrt3  1.15470053837925159591885903972e+00        /*  2.0/sqrt(3.0) */

#define e1a 1.00000000000	/* 1.0 / 1.0             */
#define e1b 0.57735026919	/* 1.0 / sqrt(3.0)       */
#define e2b 0.33333333333	/* 1.0 / 3.0             */

int o_fcn(double *obj, const double x[6])
{
  /* 29 flops */

  static double matr[4], f, t1, t2;
  static double matd[2], g;

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = (2.0*x[2] - x[1] - x[0])*sqrt3;

  matr[2] = x[3] - x[4];	/* Negate to obtain original */
  matr[3] = (2.0*x[5] - x[4] - x[3])*sqrt3;

  /* Calculate det(M). */
  t1 = matr[0]*matr[3] + matr[1]*matr[2];
  t2 = sqrt(t1*t1 + fd2);
  g = t1 + t2;

  matd[0] = matr[0] - 1.0;
  matd[1] = matr[3] - 1.0;

  /* Calculate norm(M). */
  f = matd[0]*matd[0] + matr[1]*matr[1] + 
      matr[2]*matr[2] + matd[1]*matd[1];

  /* Calculate objective function. */
  (*obj) = a * f * pow(g, b);
  return 0;
}

int g_fcn(double *obj, double g_obj[6], const double x[6])
{
  /* 137 flops */

  static double matr[4], f, t1, t2;
  static double matd[2], g;
  static double adj_m[4], loc1;

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = (2.0*x[2] - x[1] - x[0])*sqrt3;

  matr[2] = x[3] - x[4];	/* Negate to obtain original */
  matr[3] = (2.0*x[5] - x[4] - x[3])*sqrt3;

  /* Calculate det(M). */
  t1 = matr[0]*matr[3] + matr[1]*matr[2];
  t2 = sqrt(t1*t1 + fd2);
  g = t1 + t2;

  matd[0] = matr[0] - 1.0;
  matd[1] = matr[3] - 1.0;

  /* Calculate norm(M). */
  f = matd[0]*matd[0] + matr[1]*matr[1] + 
      matr[2]*matr[2] + matd[1]*matd[1];
 
  /* Calculate objective function. */
  loc1 = a * pow(g, b);
  *obj = f * loc1;

  /* Calculate the derivative of the objective function. */
  f = 2.0 * loc1;
  g = b * (*obj) / t2; 

  adj_m[0] = f*matd[0] + g*matr[3];
  adj_m[1] = -sqrt3*(f*matr[1] + g*matr[2]);
  adj_m[2] = f*matr[2] + g*matr[1];
  adj_m[3] = -sqrt3*(f*matd[1] + g*matr[0]);

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
/*****************************************************************************/

int h_fcn(double *obj, double g_obj[6], double h_obj[21], const double x[6])
{
  static double matr[4], f, t1, t2;
  static double matd[2], g, t3, loc1, loc2;
  static double adj_m[4], dg[4], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;
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

  t1 = matr[0]*dg[0] + matr[1]*dg[1];
  t2 = t1*t1 + fd2;
  t3 = sqrt(t2);
  g = t1 + t3;

  matd[0] = matr[0] - 1.0;
  matd[1] = matr[3] - 1.0;

  /* Calculate norm(M). */
  f = matd[0]*matd[0] + matr[1]*matr[1] + 
      matr[2]*matr[2] + matd[1]*matd[1];

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

  /* Gradient evaluation */
  adj_m[0] = dobj_df*matd[0] + dobj_dg*matr[3];
  adj_m[1] = -sqrt3*(dobj_df*matr[1] + dobj_dg*matr[2]);
  adj_m[2] = dobj_df*matr[2] + dobj_dg*matr[1];
  adj_m[3] = -sqrt3*(dobj_df*matd[1] + dobj_dg*matr[0]);

  g_obj[0] = adj_m[1] - adj_m[0];
  g_obj[1] = adj_m[1] + adj_m[0];
  g_obj[2] = -2.0*adj_m[1];

  g_obj[3] = adj_m[3] + adj_m[2];
  g_obj[4] = adj_m[3] - adj_m[2];
  g_obj[5] = -2.0*adj_m[3];

  /* Start of Hessian evaluation */
  matd[0] *= dobj_dfdg;
  matr[1] *= dobj_dfdg;
  matr[2] *= dobj_dfdg;
  matd[1] *= dobj_dfdg;

  /* First block of rows */
  loc1 = dobj_dgdg*dg[0] + matd[0];
  J_A[0] = dobj_df + dg[0]*(matd[0] + loc1);
  J_A[1] = e1b*(dg[0]*matr[1] + loc1*dg[1]);
  J_B[0] =      dg[0]*matr[2] + loc1*dg[2];
  J_B[1] = e1b*(dg[0]*matd[1] + loc1*dg[3] + dobj_dg);

  loc1 = dobj_dgdg*dg[1] + matr[1];
  J_A[2] = e2b*(dobj_df + dg[1]*(matr[1] + loc1));
  J_B[2] = e1b*(dg[1]*matr[2] + loc1*dg[2] + dobj_dg);
  J_B[3] = e2b*(dg[1]*matd[1] + loc1*dg[3]);

  loc1 = dobj_dgdg*dg[2] + matr[2];
  J_C[0] = dobj_df + dg[2]*(matr[2] + loc1);
  J_C[1] = e1b*(dg[2]*matd[1] + loc1*dg[3]);

  J_C[2] = e2b*(dobj_df + dg[3]*(2.0*matd[1] + dobj_dgdg*dg[3]));

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
  return 0;
}

void h_only(double h_obj[21], const double x[6])
{
  static double matr[4], f, t1, t2;
  static double matd[2], g, t3, loc1, loc2;
  static double dg[4], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;
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

  t1 = matr[0]*dg[0] + matr[1]*dg[1];
  t2 = t1*t1 + fd2;
  t3 = sqrt(t2);
  g = t1 + t3;

  matd[0] = matr[0] - 1.0;
  matd[1] = matr[3] - 1.0;

  /* Calculate norm(M). */
  f = matd[0]*matd[0] + matr[1]*matr[1] + 
      matr[2]*matr[2] + matd[1]*matd[1];

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
  matd[0] *= dobj_dfdg;
  matr[1] *= dobj_dfdg;
  matr[2] *= dobj_dfdg;
  matd[1] *= dobj_dfdg;

  /* First block of rows */
  loc1 = dobj_dgdg*dg[0] + matd[0];
  J_A[0] = dobj_df + dg[0]*(matd[0] + loc1);
  J_A[1] = e1b*(dg[0]*matr[1] + loc1*dg[1]);
  J_B[0] =      dg[0]*matr[2] + loc1*dg[2];
  J_B[1] = e1b*(dg[0]*matd[1] + loc1*dg[3] + dobj_dg);

  loc1 = dobj_dgdg*dg[1] + matr[1];
  J_A[2] = e2b*(dobj_df + dg[1]*(matr[1] + loc1));
  J_B[2] = e1b*(dg[1]*matr[2] + loc1*dg[2] + dobj_dg);
  J_B[3] = e2b*(dg[1]*matd[1] + loc1*dg[3]);

  loc1 = dobj_dgdg*dg[2] + matr[2];
  J_C[0] = dobj_df + dg[2]*(matr[2] + loc1);
  J_C[1] = e1b*(dg[2]*matd[1] + loc1*dg[3]);

  J_C[2] = e2b*(dobj_df + dg[3]*(2.0*matd[1] + dobj_dgdg*dg[3]));

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

#if 0
  printf("h(0): %5.4e\n", h_obj[0]);
  printf("h(1): %5.4e\n", h_obj[1]);
  printf("h(2): %5.4e\n", h_obj[2]);
  printf("h(3): %5.4e\n", h_obj[3]);
  printf("h(4): %5.4e\n", h_obj[4]);
  printf("h(5): %5.4e\n", h_obj[5]);
  printf("h(6): %5.4e\n", h_obj[6]);
  printf("h(7): %5.4e\n", h_obj[7]);
  printf("h(8): %5.4e\n", h_obj[8]);
  printf("h(9): %5.4e\n", h_obj[9]);
  printf("h(10): %5.4e\n", h_obj[10]);
  printf("h(11): %5.4e\n", h_obj[11]);
  printf("h(12): %5.4e\n", h_obj[12]);
  printf("h(13): %5.4e\n", h_obj[13]);
  printf("h(14): %5.4e\n", h_obj[14]);
  printf("h(15): %5.4e\n", h_obj[15]);
  printf("h(16): %5.4e\n", h_obj[16]);
  printf("h(17): %5.4e\n", h_obj[17]);
  printf("h(18): %5.4e\n", h_obj[18]);
  printf("h(19): %5.4e\n", h_obj[19]);
  exit(0);
#endif

  return;
}

