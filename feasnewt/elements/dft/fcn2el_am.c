#include <math.h>
#include <stdlib.h>
#include "fcn.h"

#define a       5.00000000000000000000000000000e-01        /*  1.0/2.0       */
#define b      -1.00000000000000000000000000000e-00        /* -1.0/1.0       */
#define bm1    -2.00000000000000000000000000000e-00        /* -2.0/1.0       */

#define d       1.0e-4					   /*  delta         */
#define fd2     4.0e-8					   /*  4.0*delta^2   */

#define sqrt3   5.77350269189625797959429519858e-01        /*  1.0/sqrt(3.0) */
#define tsqrt3  1.15470053837925159591885903972e+00        /*  2.0/sqrt(3.0) */
#define fthirds 1.33333333333333333333333333333e+00        /*  4.0/3.0       */

#if !defined(NAIVE)

#include "fcn2el_am.h"

static int g_fcnl(double *obj, double g_obj[2], 
		  const double x[6], const double t[4])
{
  static double matr[4], f, t1, t2;
  static double matd[4], g, loc1;

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = (2.0*x[2] - x[1] - x[0])*sqrt3;

  matr[2] = x[4] - x[3];
  matr[3] = (2.0*x[5] - x[4] - x[3])*sqrt3;

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

  g_obj[0] = tsqrt3*(f*matd[1] - g*matr[2]);
  g_obj[1] = tsqrt3*(f*matd[3] + g*matr[0]);
  return 0;
}

static int h_fcnl(double *obj, double g_obj[2], double h_obj[3],
		  const double x[6], const double t[4])
{
  static double matr[4], f, t1, t2;
  static double matd[4], g, t3, loc1;
  static double dg[2], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = (2.0*x[2] - x[1] - x[0])*sqrt3;

  matr[2] = x[4] - x[3];
  matr[3] = (2.0*x[5] - x[4] - x[3])*sqrt3;

  /* Calculate det(M). */
  dg[0] = -matr[2];
  dg[1] = matr[0];

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

  g_obj[0] = tsqrt3*(dobj_df*matd[1] - dobj_dg*matr[2]);
  g_obj[1] = tsqrt3*(dobj_df*matd[3] + dobj_dg*matr[0]);

  /* Start of the Hessian evaluation */
  matd[1] *= dobj_dfdg;
  matd[3] *= dobj_dfdg;

  loc1 = dobj_dgdg*dg[0] + matd[1];
  h_obj[0] = fthirds*(dobj_df + dg[0]*(matd[1] + loc1));
  h_obj[1] = fthirds*(dg[0]*matd[3] + loc1*dg[1]);

  h_obj[2] = fthirds*(dobj_df + dg[1]*(2.0*matd[3] + dobj_dgdg*dg[1]));
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
  return g_fcnl(obj, g_obj, my_x, t0);
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
  return g_fcnl(obj, g_obj, my_x, t1);
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
  return g_fcnl(obj, g_obj, x, t2);
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
  return h_fcnl(obj, g_obj, h_obj, my_x, t0);
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
  return h_fcnl(obj, g_obj, h_obj, my_x, t1);
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
  return h_fcnl(obj, g_obj, h_obj, x, t2);
#endif
}
