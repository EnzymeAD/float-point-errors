#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "fcn.h"

#define a       5.00000000000000000000000000000e-01        /*  1.0/2.0       */
#define b      -6.66666666666666666666666666667e-01        /* -2.0/3.0       */

#define d       1.0e-4					   /*  delta         */
#define fd2     4.0e-8					   /*  4.0*delta^2   */

#define sqrt3   5.77350269189625797959429519858e-01        /*  1.0/sqrt(3.0) */
#define sqrt6   4.08248290463863052509822647505e-01        /*  1.0/sqrt(6.0) */

int o_fcn_hm(double *obj, const double x[12])
{
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
  (*obj) = 1.0 / (a * f * pow(g, b));
  return 0;
}

