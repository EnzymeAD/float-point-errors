#include <math.h>
#include <stdio.h>

void getmat(double Wmat[3][3], 
            const double W0[3], const double W1[3],
	    const double W2[3], const double W3[3])
{
  Wmat[0][0] = W1[0] - W0[0];
  Wmat[1][0] = W1[1] - W0[1];
  Wmat[2][0] = W1[2] - W0[2];

  Wmat[0][1] = W2[0] - W0[0];
  Wmat[1][1] = W2[1] - W0[1];
  Wmat[2][1] = W2[2] - W0[2];

  Wmat[0][2] = W3[0] - W0[0];
  Wmat[1][2] = W3[1] - W0[1];
  Wmat[2][2] = W3[2] - W0[2];
  return;
}

void getQR(double Q[3][3], double R[3][3], double W[3][3]) 
{
  R[0][0] = sqrt(W[0][0]*W[0][0] + W[1][0]*W[1][0] + W[2][0]*W[2][0]);
  R[1][0] = 0.0;
  R[2][0] = 0.0;
  Q[0][0] = W[0][0] / R[0][0];
  Q[1][0] = W[1][0] / R[0][0];
  Q[2][0] = W[2][0] / R[0][0];

  R[0][1] = Q[0][0]*W[0][1] + Q[1][0]*W[1][1] + Q[2][0]*W[2][1];
  W[0][1] -= Q[0][0]*R[0][1];
  W[1][1] -= Q[1][0]*R[0][1];
  W[2][1] -= Q[2][0]*R[0][1];

  R[0][2] = Q[0][0]*W[0][2] + Q[1][0]*W[1][2] + Q[2][0]*W[2][2];
  W[0][2] -= Q[0][0]*R[0][2];
  W[1][2] -= Q[1][0]*R[0][2];
  W[2][2] -= Q[2][0]*R[0][2];

  R[1][1] = sqrt(W[0][1]*W[0][1] + W[1][1]*W[1][1] + W[2][1]*W[2][1]);
  R[2][1] = 0.0;
  Q[0][1] = W[0][1] / R[1][1];
  Q[1][1] = W[1][1] / R[1][1];
  Q[2][1] = W[2][1] / R[1][1];

  R[1][2] = Q[0][1]*W[0][2] + Q[1][1]*W[1][2] + Q[2][1]*W[2][2];
  W[0][2] -= Q[0][1]*R[1][2];
  W[1][2] -= Q[1][1]*R[1][2];
  W[2][2] -= Q[2][1]*R[1][2];

  R[2][2] = sqrt(W[0][2]*W[0][2] + W[1][2]*W[1][2] + W[2][2]*W[2][2]);
  Q[0][2] = W[0][2] / R[2][2];
  Q[1][2] = W[1][2] / R[2][2];
  Q[2][2] = W[2][2] / R[2][2];
  return;
}

static void divide(double W[3][3], int i, double piv)
{
  W[i][0] /= piv;
  W[i][1] /= piv;
  W[i][2] /= piv;
  return;
}

static void add(double W[3][3], int i, int j, double piv)
{
  W[i][0] -= piv*W[j][0];
  W[i][1] -= piv*W[j][1];
  W[i][2] -= piv*W[j][2];
  return;
}

static void copy(double dest[3][3],
                 double src[3][3])
{
  dest[0][0] = src[0][0];
  dest[0][1] = src[0][1];
  dest[0][2] = src[0][2];
  dest[1][0] = src[1][0];
  dest[1][1] = src[1][1];
  dest[1][2] = src[1][2];
  dest[2][0] = src[2][0];
  dest[2][1] = src[2][1];
  dest[2][2] = src[2][2];
  return;
}

static void seteye(double W[3][3])
{
  W[0][0] = 1.0;
  W[0][1] = 0.0;
  W[0][2] = 0.0;
  W[1][0] = 0.0;
  W[1][1] = 1.0;
  W[1][2] = 0.0;
  W[2][0] = 0.0;
  W[2][1] = 0.0;
  W[2][2] = 1.0;
  return;
}

void getinv(double Winv[3][3], double Wmat[3][3])
{
  double Wtmp[3][3];
  double piv;

  copy(Wtmp, Wmat);
  seteye(Winv);

  piv = Wtmp[0][0];
  divide(Winv, 0, piv);
  divide(Wtmp, 0, piv);

  piv = Wtmp[1][1];
  divide(Winv, 1, piv);
  divide(Wtmp, 1, piv);

  piv = Wtmp[2][2];
  divide(Winv, 2, piv);
  divide(Wtmp, 2, piv);

  piv = Wtmp[1][2];
  add(Winv, 1, 2, piv);
  add(Wtmp, 1, 2, piv);

  piv = Wtmp[0][2];
  add(Winv, 0, 2, piv);
  add(Wtmp, 0, 2, piv);

  piv = Wtmp[0][1];
  add(Winv, 0, 1, piv);
  add(Wtmp, 0, 1, piv);
  return;
}

void getQT(double Wmat[3][3], int idx) 
{
  double Q[3][3];
  double R[3][3];
  double Rinv[3][3];

  getQR(Q, R, Wmat);
  getinv(Rinv, R);

  printf("const double w%d[9] = {\n", idx);
  printf("  %+20.19e,\n", (fabs(Rinv[0][0]) > 1E-10) ? Rinv[0][0] : 0.0);
  printf("  %+20.19e,\n", (fabs(Rinv[0][1]) > 1E-10) ? Rinv[0][1] : 0.0);
  printf("  %+20.19e,\n", (fabs(Rinv[0][2]) > 1E-10) ? Rinv[0][2] : 0.0);
  printf("  %+20.19e,\n", (fabs(Rinv[1][0]) > 1E-10) ? Rinv[1][0] : 0.0);
  printf("  %+20.19e,\n", (fabs(Rinv[1][1]) > 1E-10) ? Rinv[1][1] : 0.0);
  printf("  %+20.19e,\n", (fabs(Rinv[1][2]) > 1E-10) ? Rinv[1][2] : 0.0);
  printf("  %+20.19e,\n", (fabs(Rinv[2][0]) > 1E-10) ? Rinv[2][0] : 0.0);
  printf("  %+20.19e,\n", (fabs(Rinv[2][1]) > 1E-10) ? Rinv[2][1] : 0.0);
  printf("  %+20.19e \n", (fabs(Rinv[2][2]) > 1E-10) ? Rinv[2][2] : 0.0);
  printf("};\n\n");

  printf("const double t%d[9] = {\n", idx);
  printf("  %+20.19e,\n", (fabs(Q[0][0]) > 1E-10) ? Q[0][0] : 0.0);
  printf("  %+20.19e,\n", (fabs(Q[0][1]) > 1E-10) ? Q[0][1] : 0.0);
  printf("  %+20.19e,\n", (fabs(Q[0][2]) > 1E-10) ? Q[0][2] : 0.0);
  printf("  %+20.19e,\n", (fabs(Q[1][0]) > 1E-10) ? Q[1][0] : 0.0);
  printf("  %+20.19e,\n", (fabs(Q[1][1]) > 1E-10) ? Q[1][1] : 0.0);
  printf("  %+20.19e,\n", (fabs(Q[1][2]) > 1E-10) ? Q[1][2] : 0.0);
  printf("  %+20.19e,\n", (fabs(Q[2][0]) > 1E-10) ? Q[2][0] : 0.0);
  printf("  %+20.19e,\n", (fabs(Q[2][1]) > 1E-10) ? Q[2][1] : 0.0);
  printf("  %+20.19e \n", (fabs(Q[2][2]) > 1E-10) ? Q[2][2] : 0.0);
  printf("};\n\n");
  return;
}

#ifdef MAIN
int main() 
{
  double W0[3];
  double W1[3];
  double W2[3];
  double W3[3];

  double Wmat[3][3];

  W0[0] = 0.0;
  W0[1] = 0.0;
  W0[2] = 0.0;

  W1[0] = 1.0;
  W1[1] = 0.0;
  W1[2] = 0.0;

  W2[0] = 1.0 / 2.0;
  W2[1] = sqrt(3.0) / 2.0;
  W2[2] = 0.0;

  W3[0] = 1.0 / 2.0;
  W3[1] = sqrt(3.0) / 6.0;
  W3[2] = sqrt(6.0) / 3.0;

  getmat(Wmat, W1, W3, W2, W0);
  getQT(Wmat, 0);

  getmat(Wmat, W0, W2, W3, W1);
  getQT(Wmat, 1);

  getmat(Wmat, W1, W0, W3, W2);
  getQT(Wmat, 2);

  getmat(Wmat, W0, W1, W2, W3);
  getQT(Wmat, 3);
  return 0;
}
#endif

