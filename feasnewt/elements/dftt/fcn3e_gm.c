#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "fcn.h"

#define a       5.00000000000000000000000000000e-01        /*  1.0/2.0       */
#define b      -1.33333333333333333333333333333e-00        /* -4.0/3.0       */

#define d	1.0e-4					   /* delta          */
#define fd2     4.0e-8                                     /*  4.0*delta^2   */

#define sqrt3   5.77350269189625797959429519858e-01        /*  1.0/sqrt(3.0) */
#define sqrt6   4.08248290463863052509822647505e-01        /*  1.0/sqrt(6.0) */

int o_fcn_gm(double *obj, const double x[12])
{
  /* 104 operations */

  static double matr[9], f, t1, t2;
  static double fmat[6], g;

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
  (*obj) = log(a) + log(f) + b*log(g);
  return 0;
}

