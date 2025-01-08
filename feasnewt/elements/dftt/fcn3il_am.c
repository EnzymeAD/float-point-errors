#include <math.h>
#include <stdlib.h>
#include "fcn.h"

#define a       5.00000000000000000000000000000e-01        /*  1.0/2.0       */
#define b      -1.33333333333333333333333333333e-00        /* -4.0/3.0       */
#define bm1    -2.33333333333333333333333333333e-00        /* -7.0/3.0       */

#define d       1.0e-4					   /*  delta         */
#define fd2     4.0e-8					   /*  4.0*delta^2   */

#if !defined(NAIVE)

#include "fcn3il_am.h"

static int g_fcnl(double *obj, double g_obj[3], 
		  const double x[12], const double w[9])
{
  /* 141 operations (104 for function) */

  static double matr[9], f, t1, t2;
  static double fmat[6], g;
  static double df[3], loc1, loc2, loc3, loc4;

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
  df[0] = fmat[0]*matr[2] + fmat[1]*matr[5] + fmat[2]*matr[8];
  df[1] = fmat[1]*matr[2] + fmat[3]*matr[5] + fmat[4]*matr[8];
  df[2] = fmat[2]*matr[2] + fmat[4]*matr[5] + fmat[5]*matr[8];

  g_obj[0] = w[8]*(f*df[0] + g*loc3);

  loc1 = g*matr[0];
  loc2 = g*matr[1];
  loc3 = g*matr[2];

  g_obj[1] = w[8]*(f*df[1] + loc2*matr[6] - loc1*matr[7]);
  g_obj[2] = w[8]*(f*df[2] + loc1*matr[4] - loc2*matr[3]);
  return 0;
}

static int h_fcnl(double *obj, double g_obj[3], double h_obj[6],
		  const double x[12], const double w[9])
{
  /* 988 operations (227 for gradient, 104 for function) */

  static double matr[9], f, t1, t2;
  static double fmat[6], g, t3;
  static double df[3], dg[3], loc1, loc2;
  static double dobj_df, dobj_dg, dobj_dfdg, dobj_dgdg;
  static double aux[6];

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
  aux[0] = matr[2]*matr[2];
  aux[1] = matr[2]*matr[5];
  aux[2] = matr[2]*matr[8];
  aux[3] = matr[5]*matr[5];
  aux[4] = matr[5]*matr[8];
  aux[5] = matr[8]*matr[8];

  /* Calculate det(M). */
  dg[0] = matr[3]*matr[7] - matr[4]*matr[6];
  dg[1] = matr[1]*matr[6] - matr[0]*matr[7];
  dg[2] = matr[0]*matr[4] - matr[1]*matr[3];

  t1 = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) +
       matr[1]*(matr[5]*matr[6] - matr[3]*matr[8]) +
       matr[2]*dg[0];
  t2 = t1*t1 + fd2;
  t3 = sqrt(t2);
  g = t1 + t3;

  fmat[0] = matr[0]*matr[0] + matr[1]*matr[1] + aux[0] - 1.0;
  fmat[1] = matr[0]*matr[3] + matr[1]*matr[4] + aux[1];
  fmat[2] = matr[0]*matr[6] + matr[1]*matr[7] + aux[2];

  fmat[3] = matr[3]*matr[3] + matr[4]*matr[4] + aux[3] - 1.0;
  fmat[4] = matr[3]*matr[6] + matr[4]*matr[7] + aux[4];

  fmat[5] = matr[6]*matr[6] + matr[7]*matr[7] + aux[5] - 1.0;

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

  df[0] = fmat[0]*matr[2] + fmat[1]*matr[5] + fmat[2]*matr[8];
  df[1] = fmat[1]*matr[2] + fmat[3]*matr[5] + fmat[4]*matr[8];
  df[2] = fmat[2]*matr[2] + fmat[4]*matr[5] + fmat[5]*matr[8];

  g_obj[0] = w[8]*(dobj_df*df[0] + dobj_dg*dg[0]);
  g_obj[1] = w[8]*(dobj_df*df[1] + dobj_dg*dg[1]);
  g_obj[2] = w[8]*(dobj_df*df[2] + dobj_dg*dg[2]);

  /* Start of the Hessian evaluation */
  t1 = w[8]*w[8];
  t2 = aux[0] + aux[3] + aux[5];

  /* Blocks for the Hessian construction */
  loc1 = dg[0]*dobj_dfdg;
  loc2 = dg[0]*dobj_dgdg + df[0]*dobj_dfdg;
  h_obj[0] = t1*(loc1*df[0] + loc2*dg[0] + dobj_df*(fmat[0] + t2 + aux[0]));
  h_obj[1] = t1*(loc1*df[1] + loc2*dg[1] + dobj_df*(fmat[1] + aux[1]));
  h_obj[2] = t1*(loc1*df[2] + loc2*dg[2] + dobj_df*(fmat[2] + aux[2]));

  loc1 = dg[1]*dobj_dfdg;
  loc2 = dg[1]*dobj_dgdg + df[1]*dobj_dfdg;
  h_obj[3] = t1*(loc1*df[1] + loc2*dg[1] + dobj_df*(fmat[3] + t2 + aux[3]));
  h_obj[4] = t1*(loc1*df[2] + loc2*dg[2] + dobj_df*(fmat[4] + aux[4]));

  loc1 = dg[2]*dobj_dfdg;
  loc2 = dg[2]*dobj_dgdg + df[2]*dobj_dfdg;
  h_obj[5] = t1*(loc1*df[2] + loc2*dg[2] + dobj_df*(fmat[5] + t2 + aux[5]));
  return 0;
}
#endif

int g_fcnl_0(double *obj, double g_obj[3], const double x[12])
{
#if defined(NAIVE)
  static double g_full[12];

  if (g_fcn(obj, g_full, x)) return 1;
  g_obj[0] = g_full[0];
  g_obj[1] = g_full[4];
  g_obj[2] = g_full[8];

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
  return g_fcnl(obj, g_obj, my_x, w0);
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
  return g_fcnl(obj, g_obj, my_x, w1);
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
  return g_fcnl(obj, g_obj, my_x, w2);
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
  return g_fcnl(obj, g_obj, x, w3);
#endif
}

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
  return h_fcnl(obj, g_obj, h_obj, my_x, w0);
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
  return h_fcnl(obj, g_obj, h_obj, my_x, w1);
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
  return h_fcnl(obj, g_obj, h_obj, my_x, w2);
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
  return h_fcnl(obj, g_obj, h_obj, x, w3);
#endif
}

#if 0
int main()
{
  double x[12];
  double o, g[12], h[78];
  double o1, g1[3], h1[6];
  int i;

  int hmap[4][6] = {
    {0, 10, 26, 42, 52, 68},
    {4, 15, 31, 46, 57, 72},
    {7, 20, 36, 49, 62, 75},
    {9, 25, 41, 51, 67, 77}
  };

  srand48(1003);

  for (i = 0; i < 12; ++i) {
    x[i] = 5*(drand48() - 0.5);
  }
  
  h_fcn(&o, g, h, x);

  g_fcnl_0(&o1, g1, x);
  if (fabs(o-o1) > 1e-10) {
    printf("o1: %5.4e\n", o-o1);
  }

  for (i = 0; i < 3; ++i) {
    if (fabs(g[4*i + 0] - g1[i]) > 1e-10) {
      printf("g2(%2d) = %5.4e\n", i, g[4*i + 0] - g1[i]);
    }
  }

  g_fcnl_1(&o1, g1, x);
  if (fabs(o-o1) > 1e-10) {
    printf("o1: %5.4e\n", o-o1);
  }

  for (i = 0; i < 3; ++i) {
    if (fabs(g[4*i + 1] - g1[i]) > 1e-10) {
      printf("g2(%2d) = %5.4e\n", i, g[4*i + 1] - g1[i]);
    }
  }

  g_fcnl_2(&o1, g1, x);
  if (fabs(o-o1) > 1e-10) {
    printf("o1: %5.4e\n", o-o1);
  }

  for (i = 0; i < 3; ++i) {
    if (fabs(g[4*i + 2] - g1[i]) > 1e-10) {
      printf("g2(%2d) = %5.4e\n", i, g[4*i + 2] - g1[i]);
    }
  }

  g_fcnl_3(&o1, g1, x);
  if (fabs(o-o1) > 1e-10) {
    printf("o1: %5.4e\n", o-o1);
  }

  for (i = 0; i < 3; ++i) {
    if (fabs(g[4*i + 3] - g1[i]) > 1e-10) {
      printf("g2(%2d) = %5.4e\n", i, g[4*i + 3] - g1[i]);
    }
  }

  h_fcnl_0(&o1, g1, h1, x);
  if (fabs(o-o1) > 1e-10) {
    printf("o1: %5.4e\n", o-o1);
  }

  for (i = 0; i < 3; ++i) {
    if (fabs(g[4*i + 0] - g1[i]) > 1e-10) {
      printf("g2(%2d) = %5.4e\n", i, g[4*i + 0] - g1[i]);
    }
  }

  for (i = 0; i < 6; ++i) {
    if (fabs(h[hmap[0][i]] - h1[i]) > 1e-10) {
      printf("h2(%2d) = %5.4e\n", i, h[hmap[0][i]] - h1[i]);
    }
  }

  h_fcnl_1(&o1, g1, h1, x);
  if (fabs(o-o1) > 1e-10) {
    printf("o1: %5.4e\n", o-o1);
  }

  for (i = 0; i < 3; ++i) {
    if (fabs(g[4*i + 1] - g1[i]) > 1e-10) {
      printf("g2(%2d) = %5.4e\n", i, g[4*i + 1] - g1[i]);
    }
  }

  for (i = 0; i < 6; ++i) {
    if (fabs(h[hmap[1][i]] - h1[i]) > 1e-10) {
      printf("h2(%2d) = %5.4e\n", i, h[hmap[1][i]] - h1[i]);
    }
  }

  h_fcnl_2(&o1, g1, h1, x);
  if (fabs(o-o1) > 1e-10) {
    printf("o1: %5.4e\n", o-o1);
  }

  for (i = 0; i < 3; ++i) {
    if (fabs(g[4*i + 2] - g1[i]) > 1e-10) {
      printf("g2(%2d) = %5.4e\n", i, g[4*i + 2] - g1[i]);
    }
  }

  for (i = 0; i < 6; ++i) {
    if (fabs(h[hmap[2][i]] - h1[i]) > 1e-10) {
      printf("h2(%2d) = %5.4e\n", i, h[hmap[2][i]] - h1[i]);
    }
  }

  h_fcnl_3(&o1, g1, h1, x);
  if (fabs(o-o1) > 1e-10) {
    printf("o1: %5.4e\n", o-o1);
  }

  for (i = 0; i < 3; ++i) {
    if (fabs(g[4*i + 3] - g1[i]) > 1e-10) {
      printf("g2(%2d) = %5.4e\n", i, g[4*i + 3] - g1[i]);
    }
  }

  for (i = 0; i < 6; ++i) {
    if (fabs(h[hmap[3][i]] - h1[i]) > 1e-10) {
      printf("h2(%2d) = %5.4e\n", i, h[hmap[3][i]] - h1[i]);
    }
  }
  return -1;
} 
#endif
