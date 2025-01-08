#include <math.h>
#include "fcn.h"

#define a       5.00000000000000000000000000000e-01        /*  1.0/2.0       */
#define b      -1.00000000000000000000000000000e-00        /* -1.0/1.0       */

#define d       1.0e-4					   /*  delta         */
#define fd2     4.0e-8					   /*  4.0*delta^2   */

#define sqrt3   5.77350269189625797959429519858e-01        /*  1.0/sqrt(3.0) */

int o_fcn_gm(double *obj, const double x[6])
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
  (*obj) = log(a) + log(f) + b*log(g);
  return 0;
}

