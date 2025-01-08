#include <math.h>
#include <stdlib.h>
#include "fcn.h"

#define a       5.00000000000000000000000000000e-01        /*  1.0/2.0       */
#define b      -1.33333333333333333333333333333e-00        /* -4.0/3.0       */
#define bm1    -2.33333333333333333333333333333e-00        /* -7.0/3.0       */

#define d       1.0e-4					   /*  delta         */
#define fd2     4.0e-8					   /*  4.0*delta^2   */

#if !defined(NAIVE)

#include "fcn2il_am.h"

static int g_fcnl(double *obj, double g_obj[2], 
		  const double x[6], const double w[4])
{
  static double matr[4], f, t1, t2;
  static double fmat[3], g;
  static double df[2], loc1;

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
  df[0] = fmat[0]*matr[1] + fmat[1]*matr[3];
  df[1] = fmat[1]*matr[1] + fmat[2]*matr[3];

  g_obj[0] = w[3]*(f*df[0] - g*matr[2]);
  g_obj[1] = w[3]*(f*df[1] + g*matr[0]);
  return 0;
}

static int h_fcnl(double *obj, double g_obj[2], double h_obj[3],
		  const double x[6], const double w[4])
{
  static double matr[4], f, t1, t2;
  static double fmat[3], g, t3;
  static double df[2], dg[2], loc1, loc2;
  static double dobj_df, dobj_dg, dobj_dfdg, dobj_dgdg;
  static double aux[3];

  /* Calculate M = A*inv(W). */
  f       = x[1] - x[0];
  matr[0] = f*w[0];
  matr[1] = f*w[1] + (x[2] - x[0])*w[3];
  
  f       = x[4] - x[3];
  matr[2] = f*w[0];
  matr[3] = f*w[1] + (x[5] - x[3])*w[3];

  /* Calculate products for M*M' */
  aux[0] = matr[1]*matr[1];
  aux[1] = matr[1]*matr[3];
  aux[2] = matr[3]*matr[3];

  /* Calculate det(M). */
  dg[0] = -matr[2];
  dg[1] = matr[0];

  t1 = matr[0]*matr[3] - matr[1]*matr[2];
  t2 = t1*t1 + fd2;
  t3 = sqrt(t2);
  g = t1 + t3;

  fmat[0] = matr[0]*matr[0] + aux[0] - 1.0;
  fmat[1] = matr[0]*matr[2] + aux[1];
  fmat[2] = matr[2]*matr[2] + aux[2] - 1.0;

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

  df[0] = fmat[0]*matr[1] + fmat[1]*matr[3];
  df[1] = fmat[1]*matr[1] + fmat[2]*matr[3];

  g_obj[0] = w[3]*(dobj_df*df[0] - dobj_dg*matr[2]);
  g_obj[1] = w[3]*(dobj_df*df[1] + dobj_dg*matr[0]);

  /* Start of the Hessian evaluation */
  t1 = w[3]*w[3];
  t2 = aux[0] + aux[2];

  /* Blocks for the Hessian construction */
  loc1 = dg[0]*dobj_dfdg;
  loc2 = dg[0]*dobj_dgdg + df[0]*dobj_dfdg;
  h_obj[0] = t1*(loc1*df[0] + loc2*dg[0] + dobj_df*(fmat[0] + t2 + aux[0]));
  h_obj[1] = t1*(loc1*df[1] + loc2*dg[1] + dobj_df*(fmat[1] + aux[1]));

  loc1 = dg[1]*dobj_dfdg; 
  loc2 = dg[1]*dobj_dgdg + df[1]*dobj_dfdg;
  h_obj[2] = t1*(loc1*df[1] + loc2*dg[1] + dobj_df*(fmat[2] + t2 + aux[2]));
  return 0;
}
#endif

int g_fcnl_0(double *obj, double g_obj[2], const double x[6])
{
#if defined(NAIVE)
  static double g_full[6];

  if (g_fcn(obj, g_full, x)) return 1;
  g_obj[0] = g_full[0];
  g_obj[1] = g_full[3];
  return 0;
#else
  static double my_x[6];

  my_x[0] = x[1];
  my_x[1] = x[2];
  my_x[2] = x[0];

  my_x[3] = x[4];
  my_x[4] = x[5];
  my_x[5] = x[3];
  return g_fcnl(obj, g_obj, my_x, w0);
#endif
}

int g_fcnl_1(double *obj, double g_obj[2], const double x[6])
{
#if defined(NAIVE)
  static double g_full[6];

  if (g_fcn(obj, g_full, x)) return 1;
  g_obj[0] = g_full[1];
  g_obj[1] = g_full[4];
  return 0;
#else
  static double my_x[6];

  my_x[0] = x[2];
  my_x[1] = x[0];
  my_x[2] = x[1];

  my_x[3] = x[5];
  my_x[4] = x[3];
  my_x[5] = x[4];
  return g_fcnl(obj, g_obj, my_x, w1);
#endif
}

int g_fcnl_2(double *obj, double g_obj[2], const double x[6])
{
#if defined(NAIVE)
  static double g_full[6];

  if (g_fcn(obj, g_full, x)) return 1;
  g_obj[0] = g_full[2];
  g_obj[1] = g_full[5];
  return 0;
#else
  return g_fcnl(obj, g_obj, x, w2);
#endif
}

int h_fcnl_0(double *obj, double g_obj[2], double h_obj[3], const double x[6])
{
#if defined(NAIVE)
  static double g_full[6];
  static double h_full[21];

  if (h_fcn(obj, g_full, h_full, x)) return 1;
  g_obj[0] = g_full[0];
  g_obj[1] = g_full[3];

  h_obj[0] = h_full[0];
  h_obj[1] = h_full[6];
  h_obj[2] = h_full[15];
  return 0;
#else
  static double my_x[6];

  my_x[0] = x[1];
  my_x[1] = x[2];
  my_x[2] = x[0];

  my_x[3] = x[4];
  my_x[4] = x[5];
  my_x[5] = x[3];
  return h_fcnl(obj, g_obj, h_obj, my_x, w0);
#endif
}

int h_fcnl_1(double *obj, double g_obj[2], double h_obj[3], const double x[6])
{
#if defined(NAIVE)
  static double g_full[6];
  static double h_full[21];

  if (h_fcn(obj, g_full, h_full, x)) return 1;
  g_obj[0] = g_full[1];
  g_obj[1] = g_full[4];

  h_obj[0] = h_full[3];
  h_obj[1] = h_full[10];
  h_obj[2] = h_full[18];
  return 0;
#else
  static double my_x[6];

  my_x[0] = x[2];
  my_x[1] = x[0];
  my_x[2] = x[1];

  my_x[3] = x[5];
  my_x[4] = x[3];
  my_x[5] = x[4];
  return h_fcnl(obj, g_obj, h_obj, my_x, w1);
#endif
}

int h_fcnl_2(double *obj, double g_obj[2], double h_obj[3], const double x[6])
{
#if defined(NAIVE)
  static double g_full[6];
  static double h_full[21];

  if (h_fcn(obj, g_full, h_full, x)) return 1;
  g_obj[0] = g_full[2];
  g_obj[1] = g_full[5];

  h_obj[0] = h_full[5];
  h_obj[1] = h_full[14];
  h_obj[2] = h_full[20];
  return 0;
#else
  return h_fcnl(obj, g_obj, h_obj, x, w2);
#endif
}

#if 0
int main()
{ 
  double x[6];
  double o, g[6], h[21];
  double o1, g1[2], h1[3];
  int i;
  
  int hmap[3][3] = {
    {0, 6, 15},
    {3, 10, 18},
    {5, 14, 20}
  };

  srand48(1003);

  for (i = 0; i < 6; ++i) {
    x[i] = 5*(drand48() - 0.5);
  }

  h_fcn(&o, g, h, x);
  
  g_fcnl_0(&o1, g1, x);
  if (fabs(o-o1) > 1e-10) {
    printf("o1: %5.4e %5.4e\n", o, o1);
  }

  for (i = 0; i < 2; ++i) {
    if (fabs(g[3*i + 0] - g1[i]) > 1e-10) {
      printf("g2(%2d) = %5.4e\n", i, g[3*i + 0] - g1[i]);
    }
  }

  g_fcnl_1(&o1, g1, x);
  if (fabs(o-o1) > 1e-10) {
    printf("o1: %5.4e\n", o-o1);
  }

  for (i = 0; i < 2; ++i) {
    if (fabs(g[3*i + 1] - g1[i]) > 1e-10) {
      printf("g2(%2d) = %5.4e\n", i, g[3*i + 1] - g1[i]);
    }
  }

  g_fcnl_2(&o1, g1, x);
  if (fabs(o-o1) > 1e-10) {
    printf("o1: %5.4e\n", o-o1);
  }

  for (i = 0; i < 2; ++i) {
    if (fabs(g[3*i + 2] - g1[i]) > 1e-10) {
      printf("g2(%2d) = %5.4e\n", i, g[3*i + 2] - g1[i]);
    }
  }

  h_fcnl_0(&o1, g1, h1, x);
  if (fabs(o-o1) > 1e-10) {
    printf("o1: %5.4e\n", o-o1);
  }

  for (i = 0; i < 2; ++i) {
    if (fabs(g[3*i + 0] - g1[i]) > 1e-10) {
      printf("g2(%2d) = %5.4e\n", i, g[3*i + 0] - g1[i]);
    }
  }

  for (i = 0; i < 3; ++i) {
    if (fabs(h[hmap[0][i]] - h1[i]) > 1e-10) {
      printf("h2(%2d) = %5.4e\n", i, h[hmap[0][i]] - h1[i]);
    }
  }

  h_fcnl_1(&o1, g1, h1, x);
  if (fabs(o-o1) > 1e-10) {
    printf("o1: %5.4e\n", o-o1);
  }

  for (i = 0; i < 2; ++i) {
    if (fabs(g[3*i + 1] - g1[i]) > 1e-10) {
      printf("g2(%2d) = %5.4e\n", i, g[3*i + 1] - g1[i]);
    }
  }

  for (i = 0; i < 3; ++i) {
    if (fabs(h[hmap[1][i]] - h1[i]) > 1e-10) {
      printf("h2(%2d) = %5.4e\n", i, h[hmap[1][i]] - h1[i]);
    }
  }

  h_fcnl_2(&o1, g1, h1, x);
  if (fabs(o-o1) > 1e-10) {
    printf("o1: %5.4e\n", o-o1);
  }

  for (i = 0; i < 2; ++i) {
    if (fabs(g[3*i + 2] - g1[i]) > 1e-10) {
      printf("g2(%2d) = %5.4e\n", i, g[3*i + 2] - g1[i]);
    }
  }

  for (i = 0; i < 3; ++i) {
    if (fabs(h[hmap[2][i]] - h1[i]) > 1e-10) {
      printf("h2(%2d) = %5.4e\n", i, h[hmap[2][i]] - h1[i]);
    }
  }
  return -1;
}
#endif
