#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "opt.h"
#include "fcn.h"
#include "pre.h"

#ifndef MICROSEC
#  define MICROSEC 1000000
#endif

#undef epsilon
#define epsilon 1.0e-10

#define QNVEC	5

/*****************************************************************************/
/* These functions are for obtaining the compacted vertex representation     */
/* of the mesh and for putting the compacted vectex representation back into */
/* the mesh structure.                                                       */
/*****************************************************************************/

static void scatterMesh(Mesh *m, const double *src) 
{
  /***************************************************************************/
  /* Take the src vertices and apply the inverse permutation to scatter the  */
  /* nonfixed coordinates into the mesh structure.                           */
  /***************************************************************************/

  double *v = m->v;
  int    *ip = m->i;

  const int nn = m->nn;
  int     i, loc;

  for (i = 0; i < nn; ++i) {
    loc = *ip++;
    v[loc  ] = src[0]; 
    v[loc+1] = src[1];
    v[loc+2] = src[2];
    src += 3;
  }
  return;
}

static void gatherMesh(double *dest, const Mesh *m)
{
  /***************************************************************************/
  /* Take the src vertices and apply the permutation to gather the nonfixed  */
  /* coordinates into the dest.                                              */
  /***************************************************************************/

  double *v = m->v;
  int    *ip = m->i;

  const int nn = m->nn;
  int     i, loc;

  for (i = 0; i < nn; ++i) {
    loc = *ip++;
    dest[0] = v[loc  ]; 
    dest[1] = v[loc+1];
    dest[2] = v[loc+2];
    dest += 3;
  }
  return;
}

/*****************************************************************************/
/* We now get into the code to the algorithm.  The first is a matrix vector  */
/* multiplication using the block structure.  This is in the reduced space.  */
/* We use pointer arithmetic to speed up the code generated by gcc.  Other   */
/* compilers would have done this automatically.                             */
/*****************************************************************************/

static double norm(const double *g, const int nn)
{
  double norm_r = 0;
  int    i;
     
  for (i = 0; i < nn; ++i) {
    norm_r += g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
    g += 3;
  }
  return norm_r;
}

static double inner(const double *g, const double *w, const int nn)
{
  double inner_r = 0;
  int    i;
     
  for (i = 0; i < nn; ++i) {
    inner_r += g[0]*w[0] + g[1]*w[1] + g[2]*w[2];
    g += 3;
    w += 3;
  }
  return inner_r;
}

static void axpy(double *r, const double *x, const double c, const double *y, 
		 const int nn)
{
  int i;

  for (i = 0; i < nn; ++i) {
    r[0] = x[0] + c*y[0];
    r[1] = x[1] + c*y[1];
    r[2] = x[2] + c*y[2];
    r += 3;
    x += 3;
    y += 3;
  }

  return;
}

static void factor(double *pd, const Mesh *mesh)
{
  double *d = mesh->dat;

  const int nn = mesh->nn;
  int i;

  for (i = 0; i < nn; ++i) {
    pd[0] = 1.0  / d[0];
    pd[1] = d[1] * pd[0];
    pd[2] = d[2] * pd[0];

    pd[3] = 1.0  / (d[3] - d[1]*pd[1]);
    pd[5] = d[4] - d[2]*pd[1];
    pd[4] = pd[3] * pd[5];

    pd[5] = 1.0  / (d[5] - d[2]*pd[2] - pd[4]*pd[5]);

#ifdef CHECK_PD
    if ((pd[0] <= 0) || (pd[3] <= 0) || (pd[5] <= 0)) {
      if (d[0] + d[3] + d[5] <= 0) {
        /* Switch to diagonal */
        pd[0] = 1.0 / fabs(d[0]);
        pd[1] = 0.0;
        pd[2] = 0.0;
        pd[3] = 1.0 / fabs(d[3]);
        pd[4] = 0.0;
        pd[5] = 1.0 / fabs(d[5]);
      }
      else {
        /* Diagonal preconditioner */
        pd[0] = 1.0 / (d[0] + d[3] + d[5]);
        pd[1] = 0.0;
        pd[2] = 0.0;
        pd[3] = pd[0];
        pd[4] = 0.0;
        pd[5] = pd[3];
      }
    }
#endif

    pd += 6;
    d  += 6;
  }
  return;
}

static void solve(double *z, double *v, double *pd, const Mesh *mesh)
{
  const int nn = mesh->nn;
  int i;

  for (i = 0; i < nn; ++i) {
    z[0] = v[0];
    z[1] = v[1] - pd[1]*z[0];
    z[2] = v[2] - pd[2]*z[0] - pd[4]*z[1];

    z[0] *= pd[0];
    z[1] *= pd[3];
    z[2] *= pd[5];

    z[1] -= pd[4]*z[2];
    z[0] -= pd[1]*z[1] + pd[2]*z[2];

    pd += 6;
    v  += 3;
    z  += 3;
  }
  return;
}

int optMesh(Mesh *mesh, int max_iter, double cn_tol, int precond)
{
  double **w, *wtmp;
  double **v, *vtmp;
  double *d;
  double *x;
  double *a, *b, *r;
  double *dat;

  const double sigma = 1e-4;
  const double beta0 = 0.25;
  const double beta1 = 0.80;
  const double tol1 = 1e-8;

  double norm_r, norm_g;
  double alpha, beta;
  double obj, objn;

  const int nn = mesh->nn;
  int i, iter = 0;
  int nf = 1, ng = 1, nh = 0;
  int ferr;

  double m, m1, m2;
  double s, s1, s2;
  double t1, t2;

  struct rusage r0, r1;

#ifdef DISPLAY_MAX
  double fmax;
  int    fidx;
#endif

  if (nn <= 0) {
    /* No nodes!  Just return.                                               */
    return 0;
  }

  getrusage(RUSAGE_SELF, &r0);

#ifdef USE_WEIGHT
  mesh->w = (double *)malloc(16*sizeof(double)*mesh->ne);

  {
    double *weight = mesh->w;

    double W0[3];
    double W1[3];
    double W2[3];
    double W3[3];
    double Wmat[3][3];
    double Q[3][3];
    double R[3][3];
    double Rinv[3][3];

    int i;

    srand48(1001);

    for (i = 0; i < mesh->ne; ++i) {

#ifdef RANDOM_WEIGHT
      W0[0] = 0.0;
      W0[1] = 0.0;
      W0[2] = 0.0;

      W1[0] = 5.0*drand48() + 0.5;
      W1[1] = 0.0;
      W1[2] = 0.0;

      W2[0] = 10.0*drand48() - 5.0;
      W2[1] = 5.0*drand48() + 0.5;
      W2[2] = 0.0;

      W3[0] = 10.0*drand48() - 5.0;
      W3[1] = 10.0*drand48() - 5.0;
      W3[2] = 5.0*drand48() + 0.5;
#else
      W0[0] = 0.0;
      W0[1] = 0.0;
      W0[2] = 0.0;

      W1[0] = 1.0;
      W1[1] = 0.0;
      W1[2] = 0.0;

      W2[0] = 0.5;
      W2[1] = 0.5 * sqrt(3.0);
      W2[2] = 0.0;

      W3[0] = 0.5;
      W3[1] = 0.5 / sqrt(3.0);
      W3[2] = sqrt(6.0) / 3.0;
#endif

      getmat(Wmat, W0, W1, W2, W3);
      getQR(Q, R, Wmat);
      getinv(Rinv, R);

      weight[0] = (fabs(Rinv[0][0]) > 1E-10) ? Rinv[0][0] : 0.0;
      weight[1] = (fabs(Rinv[0][1]) > 1E-10) ? Rinv[0][1] : 0.0;
      weight[2] = (fabs(Rinv[0][2]) > 1E-10) ? Rinv[0][2] : 0.0;
      weight[3] = (fabs(Rinv[1][1]) > 1E-10) ? Rinv[1][1] : 0.0;
      weight[4] = (fabs(Rinv[1][2]) > 1E-10) ? Rinv[1][2] : 0.0;
      weight[5] = (fabs(Rinv[2][2]) > 1E-10) ? Rinv[2][2] : 0.0;

      weight[6] = (fabs(Q[0][0]) > 1E-10) ? Q[0][0] : 0.0;
      weight[7] = (fabs(Q[0][1]) > 1E-10) ? Q[0][1] : 0.0;
      weight[8] = (fabs(Q[0][2]) > 1E-10) ? Q[0][2] : 0.0;
      weight[9] = (fabs(Q[1][0]) > 1E-10) ? Q[1][0] : 0.0;
      weight[10] = (fabs(Q[1][1]) > 1E-10) ? Q[1][1] : 0.0;
      weight[11] = (fabs(Q[1][2]) > 1E-10) ? Q[1][2] : 0.0;
      weight[12] = (fabs(Q[2][0]) > 1E-10) ? Q[2][0] : 0.0;
      weight[13] = (fabs(Q[2][1]) > 1E-10) ? Q[2][1] : 0.0;
      weight[14] = (fabs(Q[2][2]) > 1E-10) ? Q[2][2] : 0.0;

      weight += 16;
    }
  }
#endif
 
  hMesh(mesh);
  if (gFcn(&obj, mesh)) {
    fprintf(stderr, "Invalid starting point.\n");
    exit(-1);
  }

  a = (double *)calloc(QNVEC, sizeof(double));
  b = (double *)calloc(QNVEC, sizeof(double));
  r = (double *)calloc(QNVEC, sizeof(double));

  w = (double **)malloc((QNVEC+1)*sizeof(double *));
  v = (double **)malloc((QNVEC+1)*sizeof(double *));
  for (i = 0; i <= QNVEC; ++i) {
    w[i] = (double *)calloc(3*nn, sizeof(double));
    v[i] = (double *)calloc(3*nn, sizeof(double));
  }

  d = (double *)malloc(3*sizeof(double)*nn);
  x = (double *)malloc(3*sizeof(double)*nn);
  dat = (double *)malloc(6*sizeof(double)*nn);

  gatherMesh(x, mesh);
  norm_r = norm(mesh->g, nn);
  norm_g = sqrt(norm_r);

  getrusage(RUSAGE_SELF, &r1);

  m1 = (double) r1.ru_utime.tv_usec;
  m2 = (double) r0.ru_utime.tv_usec;
  m = m1 - m2;
    
  s1 = (double) r1.ru_utime.tv_sec;
  s2 = (double) r0.ru_utime.tv_sec;
  s = s1 - s2;

  t1 = s + m / MICROSEC;

  m1 = (double) r1.ru_stime.tv_usec;
  m2 = (double) r0.ru_stime.tv_usec;
  m = m1 - m2;

  s1 = (double) r1.ru_stime.tv_sec;
  s2 = (double) r0.ru_stime.tv_sec;
  s = s1 - s2;

  t2 = s + m / MICROSEC;

#ifdef DISPLAY_MAX
  oMax(&fmax, &fidx, mesh);
  printf("%4d I %10.9e %5.4e            %5.4e %5d "
         "usr: %5.4e sys: %5.4e tot: %5.4e\n",
         iter, obj, norm_g, fmax, fidx, t1, t2, t1+t2);
#else
  printf("%4d I %10.9e %5.4e            "
         "usr: %5.4e sys: %5.4e tot: %5.4e\n",
         iter, obj, norm_g, t1, t2, t1+t2);
#endif


  while ((norm_g > cn_tol) && (iter < max_iter)) {
    getrusage(RUSAGE_SELF, &r0);
    ++iter;

    ++nh;
    hOnly(mesh);

    /* Store the current iterate and gradient */
    memcpy(w[QNVEC],       x, 3*sizeof(double)*nn);
    memcpy(v[QNVEC], mesh->g, 3*sizeof(double)*nn);

    memcpy(x, mesh->g, 3*sizeof(double)*nn);
    for (i = QNVEC - 1; i >= 0; --i) {
      a[i] = r[i] * inner(w[i], x, nn);
      axpy(x, x, -a[i], v[i], nn);
    }

    factor(dat, mesh);
    solve(d, x, dat, mesh);
    // dat = mesh->dat;
    // wtmp = d;
    // vtmp = x;
    // for (i = 0; i < nn; ++i) {
    //   wtmp[0] = vtmp[0] / dat[0];
    //   wtmp[1] = vtmp[1] / dat[3];
    //   wtmp[2] = vtmp[2] / dat[5];
    //   wtmp += 3; vtmp += 3; dat += 6;
    // }
    // memcpy(d, x, 3*sizeof(double)*nn);

    for (i = 0; i < QNVEC; ++i) {
      b[i] = r[i] * inner(v[i], d, nn);
      axpy(d, d, a[i] - b[i], w[i], nn);
    }

    alpha = -inner(mesh->g, d, nn);  /* direction is negated */
    if (alpha > 0.0) {
      printf("Not descent.\n");
      exit(-1);
    }
    else {
      alpha *= sigma;
    }
    beta = 1.0;

    /* Unrolled for better performance.  Only do a gradient evaluation when  */
    /* beta = 1.0.  Do not do this at other times.                           */

    axpy(x, w[QNVEC], -beta, d, nn);
    scatterMesh(mesh, x);

    ++nf;
    ++ng;
    ferr = gFcn(&objn, mesh);
    if ((!ferr && (obj - objn >= -alpha*beta - epsilon)) ||
	(!ferr && (sqrt(norm(mesh->g, nn)) < 100*cn_tol))) {
      /* Iterate is acceptable */
    }
    else {
      if (ferr) {
	/* Function not defined at trial point */
	beta *= beta0;
      }
      else {
	/* Iterate not acceptable */
	beta *= beta1;
      }

      while (beta >= tol1) {
        axpy(x, w[QNVEC], -beta, d, nn);
        scatterMesh(mesh, x);

        ++nf;
        if (oFcn(&objn, mesh)) {
	  /* Function not defined at trial point */
	  beta *= beta0;
        }
        else if (obj - objn >= -alpha*beta - epsilon) {
	  /* Iterate is acceptable */
	  break;
        }
        else {
	  /* Iterate not acceptable */
	  beta *= beta1;
        }
      }

      if (beta < tol1) {
        printf("Newton step not good\n");
        exit(-1);
      }

      ++ng;
      gOnly(mesh);
    }

    wtmp = w[0];
    vtmp = v[0];
    for (i = 0; i < QNVEC-1; ++i) {
      r[i] = r[i+1];
      w[i] = w[i+1];
      v[i] = v[i+1];
    }
    w[i] = wtmp;
    v[i] = vtmp;

    axpy(w[QNVEC-1],       x, -1.0, w[QNVEC], nn);
    axpy(v[QNVEC-1], mesh->g, -1.0, v[QNVEC], nn);
    r[QNVEC-1] = 1.0 / inner(w[QNVEC-1], v[QNVEC-1], nn);

    /* Print convergence statistics */
    obj = objn;
    norm_r = norm(mesh->g, nn);
    norm_g = sqrt(norm_r);

    getrusage(RUSAGE_SELF, &r1);

    m1 = (double) r1.ru_utime.tv_usec;
    m2 = (double) r0.ru_utime.tv_usec;
    m = m1 - m2;
    
    s1 = (double) r1.ru_utime.tv_sec;
    s2 = (double) r0.ru_utime.tv_sec;
    s = s1 - s2;
    
    t1 = s + m / MICROSEC;
    
    m1 = (double) r1.ru_stime.tv_usec;
    m2 = (double) r0.ru_stime.tv_usec;
    m = m1 - m2;
    
    s1 = (double) r1.ru_stime.tv_sec;
    s2 = (double) r0.ru_stime.tv_sec;
    s = s1 - s2;
    
    t2 = s + m / MICROSEC;
    
#ifdef DISPLAY_MAX
    oMax(&fmax, &fidx, mesh);
    printf("%4d N %10.9e %5.4e %5.4e %5.4e %5d "
           "usr: %5.4e sys: %5.4e tot: %5.4e\n",
           iter, obj, norm_g, beta, fmax, fidx, t1, t2, t1+t2);
#else
    printf("%4d N %10.9e %5.4e %5.4e usr: %5.4e sys: %5.4e tot: %5.4e\n",
           iter, obj, norm_g, beta, t1, t2, t1+t2);
#endif
  }

  printf("Function evals: %4d\nGradient evals: %4d\nHessian  evals: %4d\n", 
	 nf, ng, nh);

  for (i = 0; i <= QNVEC; ++i) {
    free(w[i]);
    free(v[i]);
  }
  free(w);
  free(v);

  free(x);
  free(d);
  free(dat);

  free(a);
  free(b);
  free(r);
  return 0;
}
