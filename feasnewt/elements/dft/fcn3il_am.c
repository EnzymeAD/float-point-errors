#include <math.h>
#include <stdlib.h>
#include "fcn.h"

#define a       5.00000000000000000000000000000e-01        /*  1.0/2.0       */
#define b      -6.66666666666666666666666666667e-01        /* -2.0/3.0       */
#define bm1    -1.66666666666666666666666666667e-00        /* -5.0/3.0       */

#define d       1.0e-4					   /*  delta         */
#define fd2     4.0e-8					   /*  4.0*delta^2   */

#if !defined(NAIVE)

#include "fcn3il_am.h"

static int g_fcnl(double *obj, double g_obj[3], 
		  const double x[12], const double w[9], const double t[9])
{
  /* 104 operations (83 for function) */

  static double matr[9], f, t1, t2;
  static double matd[9], g, loc1, loc2;

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
  loc1 = matr[3]*matr[7] - matr[4]*matr[6];
  t1 = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) + 
       matr[1]*(matr[5]*matr[6] - matr[3]*matr[8]) +
       matr[2]*loc1;
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
  loc2 = a * pow(g, b);
  *obj = f * loc2;

  /* Calculate the derivative of the objective function. */
  f = 2.0 * loc2;
  g = b * (*obj) / t2;

  g_obj[0] = w[8]*(f*matd[2] + g*loc1);

  loc1 = g*matr[0];
  loc2 = g*matr[1];

  g_obj[1] = w[8]*(f*matd[5] + loc2*matr[6] - loc1*matr[7]);
  g_obj[2] = w[8]*(f*matd[8] + loc1*matr[4] - loc2*matr[3]);
  return 0;
}

static int h_fcnl(double *obj, double g_obj[3], double h_obj[6],
		  const double x[12], const double w[9], const double t[9])
{
  /* 146 operations (104 for gradient, 83 for function) */

  static double matr[9], f, t1, t2;
  static double matd[9], g, t3, loc1;
  static double dg[3], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;

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
  dg[0] = matr[3]*matr[7] - matr[4]*matr[6];
  dg[1] = matr[1]*matr[6] - matr[0]*matr[7];
  dg[2] = matr[0]*matr[4] - matr[1]*matr[3];

  t1 = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) + 
       matr[1]*(matr[5]*matr[6] - matr[3]*matr[8]) +
       matr[2]*dg[0];
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

  g_obj[0] = w[8]*(dobj_df*matd[2] + dobj_dg*dg[0]);
  g_obj[1] = w[8]*(dobj_df*matd[5] + dobj_dg*dg[1]);
  g_obj[2] = w[8]*(dobj_df*matd[8] + dobj_dg*dg[2]);

  /* Start of the Hessian evaluation */
  t1 = w[8]*w[8];

  matd[2] *= dobj_dfdg;
  matd[5] *= dobj_dfdg;
  matd[8] *= dobj_dfdg;

  loc1 = dobj_dgdg*dg[0] + matd[2];
  h_obj[0] = t1*(dobj_df + dg[0]*(matd[2] + loc1));
  h_obj[1] = t1*(dg[0]*matd[5] + loc1*dg[1]);
  h_obj[2] = t1*(dg[0]*matd[8] + loc1*dg[2]);

  loc1 = dobj_dgdg*dg[1] + matd[5];
  h_obj[3] = t1*(dobj_df + dg[1]*(matd[5] + loc1));
  h_obj[4] = t1*(dg[1]*matd[8] + loc1*dg[2]);

  h_obj[5] = t1*(dobj_df + dg[2]*(2.0*matd[8] + dobj_dgdg*dg[2]));
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
  return g_fcnl(obj, g_obj, my_x, w0, t0);
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
  return g_fcnl(obj, g_obj, my_x, w1, t1);
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
  return g_fcnl(obj, g_obj, my_x, w2, t2);
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
  return g_fcnl(obj, g_obj, x, w3, t3);
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
  return h_fcnl(obj, g_obj, h_obj, my_x, w0, t0);
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
  return h_fcnl(obj, g_obj, h_obj, my_x, w1, t1);
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
  return h_fcnl(obj, g_obj, h_obj, my_x, w2, t2);
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
  return h_fcnl(obj, g_obj, h_obj, x, w3, t3);
#endif
}
