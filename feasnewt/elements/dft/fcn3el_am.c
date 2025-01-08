#include <math.h>
#include <stdlib.h>
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

#if !defined(NAIVE)

#include "fcn3el_am.h"

static int g_fcnl(double *obj, double g_obj[3], 
		  const double x[12], const double t[9])
{
  /* 95 flops */

  static double matr[9], f, t1, t2;
  static double matd[9], g, loc1, loc2;

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

  g_obj[0] = tsqrt6*(f*matd[2] + g*loc1);

  loc1 = g*matr[0];
  loc2 = g*matr[1];

  g_obj[1] = tsqrt6*(f*matd[5] + loc2*matr[6] - loc1*matr[7]);
  g_obj[2] = tsqrt6*(f*matd[8] + loc1*matr[4] - loc2*matr[3]);
  return 0;
}

static int h_fcnl(double *obj, double g_obj[3], double h_obj[6],
		  const double x[12], const double t[9])
{
  /* 136 flops */

  static double matr[9], f, t1, t2;
  static double matd[9], g, t3, loc1;
  static double dg[3], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;

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

  g_obj[0] = tsqrt6*(dobj_df*matd[2] + dobj_dg*dg[0]);
  g_obj[1] = tsqrt6*(dobj_df*matd[5] + dobj_dg*dg[1]);
  g_obj[2] = tsqrt6*(dobj_df*matd[8] + dobj_dg*dg[2]);

  /* Start of Hessian evaluation */
  matd[2] *= dobj_dfdg;
  matd[5] *= dobj_dfdg;
  matd[8] *= dobj_dfdg;

  /* Blocks for the Hessian construction */
  loc1 = dobj_dgdg*dg[0] + matd[2];
  h_obj[0] = 9.0*e3c*(dobj_df + dg[0]*(matd[2] + loc1));
  h_obj[1] = 9.0*e3c*(dg[0]*matd[5] + loc1*dg[1]);
  h_obj[2] = 9.0*e3c*(dg[0]*matd[8] + loc1*dg[2]);

  loc1 = dobj_dgdg*dg[1] + matd[5];
  h_obj[3] = 9.0*e3c*(dobj_df + dg[1]*(matd[5] + loc1));
  h_obj[4] = 9.0*e3c*(dg[1]*matd[8] + loc1*dg[2]);

  h_obj[5] = 9.0*e3c*(dobj_df + dg[2]*(2.0*matd[8] + dobj_dgdg*dg[2]));
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
  return g_fcnl(obj, g_obj, my_x, t0);
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
  return g_fcnl(obj, g_obj, my_x, t1);
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
  return g_fcnl(obj, g_obj, my_x, t2);
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
  return g_fcnl(obj, g_obj, x, t3);
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
  return h_fcnl(obj, g_obj, h_obj, my_x, t0);
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
  return h_fcnl(obj, g_obj, h_obj, my_x, t1);
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
  return h_fcnl(obj, g_obj, h_obj, my_x, t2);
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
  return h_fcnl(obj, g_obj, h_obj, x, t3);
#endif
}

