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

/* Assume that w is upper triangular. The t variable is the multiplication   */
/* of the rotation matrices required to obtain this matrix.                  */

/* const double w[4] = {1, -sqrt3, 0, 2*sqrt3}; */
const double w[4] = {1, 0, 0, 1};
const double t[4] = {1, 0, 0, 1};

int o_fcn(double *obj, const double x[6])
{
  static double matr[4], f, t1, t2;
  static double fmat[3], g;

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

  /* Calculate norm(M). */
  fmat[0] = matr[0]*matr[0] + matr[1]*matr[1] - 1.0;
  fmat[1] = matr[0]*matr[2] + matr[1]*matr[3];
  fmat[2] = matr[2]*matr[2] + matr[3]*matr[3] - 1.0;

  f = fmat[0]*fmat[0] + 2.0*fmat[1]*fmat[1] +
                            fmat[2]*fmat[2];

  /* Calculate objective function. */
  (*obj) = a * f * pow(g, b);
  return 0;
}

/*****************************************************************************/
/* Derivative calculation.                                                   */
/*****************************************************************************/

int g_fcn(double *obj, double g_obj[6], const double x[6])
{
  static double matr[4], f, t1, t2;
  static double fmat[3], g;
  static double adj_m[4], df[4], loc1;

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

  /* Calculate norm(M). */
  fmat[0] = matr[0]*matr[0] + matr[1]*matr[1] - 1.0;
  fmat[1] = matr[0]*matr[2] + matr[1]*matr[3];
  fmat[2] = matr[2]*matr[2] + matr[3]*matr[3] - 1.0;

  f = fmat[0]*fmat[0] + 2.0*fmat[1]*fmat[1] +
                            fmat[2]*fmat[2];

  /* Calculate objective function. */
  loc1 = a * pow(g, b);
  *obj = f * loc1;

  /* Calculate the derivative of the objective function. */
  f = 4.0 * loc1;                        /* Constant on nabla f */
  g = b * (*obj) / t2;                   /* Constant on nabla g */

  /* Compute d fmat by d mat */
  df[0] = fmat[0]*matr[0] + fmat[1]*matr[2];
  df[1] = fmat[0]*matr[1] + fmat[1]*matr[3];
  df[2] = fmat[1]*matr[0] + fmat[2]*matr[2];
  df[3] = fmat[1]*matr[1] + fmat[2]*matr[3];

  adj_m[0] = f*df[0] + g*matr[3];
  adj_m[1] = f*df[1] - g*matr[2];
  adj_m[2] = f*df[2] - g*matr[1];
  adj_m[3] = f*df[3] + g*matr[0];

  g_obj[1] = w[0]*adj_m[0] + w[1]*adj_m[1];
  g_obj[2] =                 w[3]*adj_m[1];
  g_obj[0] = -g_obj[1] - g_obj[2];

  g_obj[4] = w[0]*adj_m[2] + w[1]*adj_m[3];
  g_obj[5] =                 w[3]*adj_m[3];
  g_obj[3] = -g_obj[4] - g_obj[5] ;
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
  /* 988 operations (227 for gradient, 104 for function) */

  static double matr[4], f, t1, t2;
  static double fmat[3], ftmat[3], g, t3;
  static double adj_m[4], df[4], dg[4], loc1, loc2;
  static double dobj_df, dobj_dg, dobj_dfdg, dobj_dgdg;
  static double J_A[3], J_B[4], J_D[3];
  static double aux[10];

  /* Calculate M = A*inv(W). */
  f       = x[1] - x[0];
  matr[0] = f*w[0];
  matr[1] = f*w[1] + (x[2] - x[0])*w[3];
  
  f       = x[4] - x[3];
  matr[2] = f*w[0];
  matr[3] = f*w[1] + (x[5] - x[3])*w[3];

  /* Calculate products for M*M' */
  aux[0] = matr[0]*matr[0];
  aux[1] = matr[0]*matr[1];
  aux[2] = matr[0]*matr[2];
  aux[3] = matr[0]*matr[3];
  aux[4] = matr[1]*matr[1];
  aux[5] = matr[1]*matr[2];
  aux[6] = matr[1]*matr[3];
  aux[7] = matr[2]*matr[2];
  aux[8] = matr[2]*matr[3];
  aux[9] = matr[3]*matr[3];

  /* Calculate det(M). */
  dg[0] = matr[3];
  dg[1] = -matr[2];
  dg[2] = -matr[1];
  dg[3] = matr[0];

  t1 = matr[0]*matr[3] - matr[1]*matr[2];
  t2 = t1*t1 + fd2;
  t3 = sqrt(t2);
  g = t1 + t3;

  fmat[0] = aux[0] + aux[4] - 1.0;
  fmat[1] = aux[2] + aux[6];
  fmat[2] = aux[7] + aux[9] - 1.0;

  f = fmat[0]*fmat[0] + 2.0*fmat[1]*fmat[1] +
                            fmat[2]*fmat[2];

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

  df[0] = fmat[0]*matr[0] + fmat[1]*matr[2];
  df[1] = fmat[0]*matr[1] + fmat[1]*matr[3];
  df[2] = fmat[1]*matr[0] + fmat[2]*matr[2];
  df[3] = fmat[1]*matr[1] + fmat[2]*matr[3];

  adj_m[0] = dobj_df*df[0] + dobj_dg*dg[0];
  adj_m[1] = dobj_df*df[1] + dobj_dg*dg[1];
  adj_m[2] = dobj_df*df[2] + dobj_dg*dg[2];
  adj_m[3] = dobj_df*df[3] + dobj_dg*dg[3];

  g_obj[1] = w[0]*adj_m[0] + w[1]*adj_m[1];
  g_obj[2] =                 w[3]*adj_m[1];
  g_obj[0] = -g_obj[1] - g_obj[2];

  g_obj[4] = w[0]*adj_m[2] + w[1]*adj_m[3];
  g_obj[5] =                 w[3]*adj_m[3];
  g_obj[3] = -g_obj[4] - g_obj[5];

  /* Start of the Hessian evaluation */
  ftmat[0] = aux[0] + aux[7];
  ftmat[1] = aux[1] + aux[8];
  ftmat[2] = aux[4] + aux[9];

  adj_m[0] = dobj_dg*matr[0];
  adj_m[1] = dobj_dg*matr[1];
  adj_m[2] = dobj_dg*matr[2];
  adj_m[3] = dobj_dg*matr[3];

  /* Blocks for the Hessian construction */
  loc1 = dg[0]*dobj_dfdg;
  loc2 = dg[0]*dobj_dgdg + df[0]*dobj_dfdg;
  J_A[0] = loc1*df[0] + loc2*dg[0] + dobj_df*(fmat[0] + ftmat[0] + aux[0]);
  J_A[1] = loc1*df[1] + loc2*dg[1] + dobj_df*(          ftmat[1] + aux[1]);
  J_B[0] = loc1*df[2] + loc2*dg[2] + dobj_df*(fmat[1] + aux[2]);
  J_B[1] = loc1*df[3] + loc2*dg[3] + dobj_df*aux[5] + dobj_dg;

  loc1 = dg[1]*dobj_dfdg;
  loc2 = dg[1]*dobj_dgdg + df[1]*dobj_dfdg;
  J_A[2] = loc1*df[1] + loc2*dg[1] + dobj_df*(fmat[0] + ftmat[2] + aux[4]);
  J_B[2] = loc1*df[2] + loc2*dg[2] + dobj_df*aux[3] - dobj_dg;
  J_B[3] = loc1*df[3] + loc2*dg[3] + dobj_df*(fmat[1] + aux[6]);

  loc1 = dg[2]*dobj_dfdg;
  loc2 = dg[2]*dobj_dgdg + df[2]*dobj_dfdg;
  J_D[0] = loc1*df[2] + loc2*dg[2] + dobj_df*(fmat[2] + ftmat[0] + aux[7]);
  J_D[1] = loc1*df[3] + loc2*dg[3] + dobj_df*(          ftmat[1] + aux[8]);

  loc1 = dg[3]*dobj_dfdg; 
  loc2 = dg[3]*dobj_dgdg + df[3]*dobj_dfdg;
  J_D[2] = loc1*df[3] + loc2*dg[3] + dobj_df*(fmat[2] + ftmat[2] + aux[9]);  

  /* Put everything together */
  assemble_diag(h_obj, J_A);
  assemble_diag(h_obj + 15, J_D);
  assemble_offdiag(h_obj + 6, J_B);
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
  double obj, g[6], h[21], x[6];
  int i;

  srand48(1003);

  for (i = 0; i < 6; ++i) {
    x[i] = 5*(drand48() - 0.5);
    printf("x(%2d) := %10.9e\n", i, x[i]);
  }

  h_fcn(&obj, g, h, x);

  printf("o     = %10.9e\n", obj);

  for (i = 0; i < 6; ++i) {
    printf("g(%2d) = %10.9e\n", i, g[i]);
  }

  for (i = 0; i < 21; ++i) {
    printf("h(%2d) = %10.9e\n", i, h[i]);
  }
  return -1;
}
#endif

