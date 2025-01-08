#include <math.h>
#include "fcn.h"

/*****************************************************************************/
/* This set of functions reference tetrahedral elements to an regular        */
/* tetrahedron.  The input are the coordinates in the following order:       */
/*      [x1 x2 x3 x4 y1 y2 y3 y4 z1 z2 z3 z4]                                */
/* A zero return value indicates success, while a nonzero value indicates    */
/* failure.                                                                  */
/*****************************************************************************/
/* Not all compilers substitute out constants (especially the square root).  */
/* Therefore, they are substituted out manually.  The values below were      */
/* calculated on a solaris machine using long doubles. I believe they are    */
/* accurate.                                                                 */
/*****************************************************************************/

#define sqrt3   5.77350269189625797959429519858e-01        /*  1.0/sqrt(3.0) */
#define tsqrt3  1.15470053837925159591885903972e+00        /*  2.0/sqrt(3.0) */
#define sqrt6   4.08248290463863052509822647505e-01        /*  1.0/sqrt(6.0) */
#define tsqrt6  1.22474487139158915752946794252e+00        /*  3.0/sqrt(6.0) */
#define a       1.00000000000000000000000000000e-00        /*  1.0/1.0       */
#define b      -1.00000000000000000000000000000e-00        /* -1.0/1.0       */
#define bm1    -2.00000000000000000000000000000e-00        /* -2.0/1.0       */

int o_cond(double *obj, const double x[12])
{
  static double matr[9], cftr[9], f;
  static double g;

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
  g = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) +
      matr[1]*(matr[5]*matr[6] - matr[3]*matr[8]) +
      matr[2]*(matr[3]*matr[7] - matr[4]*matr[6]);
  if (g <= epsilon) { *obj = g; return 1; }

  /* Calculate cofactor(M). */

  cftr[0] = matr[4]*matr[8] - matr[5]*matr[7];
  cftr[1] = matr[3]*matr[8] - matr[5]*matr[6];
  cftr[2] = matr[3]*matr[7] - matr[4]*matr[6];

  cftr[3] = matr[1]*matr[8] - matr[2]*matr[7];
  cftr[4] = matr[0]*matr[8] - matr[2]*matr[6];
  cftr[5] = matr[0]*matr[7] - matr[1]*matr[6];

  cftr[6] = matr[1]*matr[5] - matr[2]*matr[4];
  cftr[7] = matr[0]*matr[5] - matr[2]*matr[3];
  cftr[8] = matr[0]*matr[4] - matr[1]*matr[3];

  /* Calculate norm(M). */
  f = sqrt(matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] +
           matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] +
           matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8]) *
      sqrt(cftr[0]*cftr[0] + cftr[1]*cftr[1] + cftr[2]*cftr[2] +
	   cftr[3]*cftr[3] + cftr[4]*cftr[4] + cftr[5]*cftr[5] +
           cftr[6]*cftr[6] + cftr[7]*cftr[7] + cftr[8]*cftr[8]);

  /* Calculate objective function. */
  (*obj) = a * f * pow(g, b);
  return 0;
}

