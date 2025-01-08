#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "fcn.h"

#define a       5.00000000000000000000000000000e-01        /*  1.0/2.0       */
#define b      -6.66666666666666666666666666667e-01        /* -2.0/3.0       */

#define d       1.0e-4					   /*  delta         */
#define fd2     4.0e-8					   /*  4.0*delta^2   */

int o_fcn_hm(double *obj, const double x[12])
{
  /* 50 operations */

  static double matr[9], f, t1, t2;
  static double matd[9], g;

  /* Calculate M = A*inv(W). */
  matr[0] = x[1] - x[0];
  matr[1] = x[2] - x[0];
  matr[2] = x[3] - x[0];

  matr[3] = x[5] - x[4];
  matr[4] = x[6] - x[4];
  matr[5] = x[7] - x[4];

  matr[6] = x[9] - x[8];
  matr[7] = x[10] - x[8];
  matr[8] = x[11] - x[8];

  /* Calculate det(M). */
  t1 = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) +
       matr[1]*(matr[5]*matr[6] - matr[3]*matr[8]) +
       matr[2]*(matr[3]*matr[7] - matr[4]*matr[6]);
  t2 = sqrt(t1*t1 + fd2);
  g = t1 + t2;

  matd[0] = matr[0] - 1.0;
  matd[1] = matr[1];
  matd[2] = matr[2];
  matd[3] = matr[3];
  matd[4] = matr[4] - 1.0;
  matd[5] = matr[5];
  matd[6] = matr[6];
  matd[7] = matr[7];
  matd[8] = matr[8] - 1.0;

  /* Calculate norm(M). */
  f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
      matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
      matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

  /* Calculate objective function. */
  (*obj) = 1.0 / (a * f * pow(g, b));
  return 0;
}

