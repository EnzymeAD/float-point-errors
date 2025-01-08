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

static void matmul(double *w, const Mesh *mesh, const double *p)
{
  int    *len = mesh->len;
  int    *col = mesh->col;
  double *dat = mesh->dat;
  double *n;
  double *o;

  double x[3];
  double y[3];
  double m[3];

  const int nn = mesh->nn;
  int c;
  int i, j, l;

  for (i = 0; i < nn; ++i) {
    /* Get x components and modification for diagonal block */
    l = *len++;
    c = *col++;
    o = w + c;

    x[0] = p[c];
    x[1] = p[c+1];
    x[2] = p[c+2];

    m[0] = dat[0]*x[0] + dat[1]*x[1] + dat[2]*x[2];
    m[1] = dat[1]*x[0] + dat[3]*x[1] + dat[4]*x[2];
    m[2] = dat[2]*x[0] + dat[4]*x[1] + dat[5]*x[2];

    dat += 6;

    /* Calculate the off diagonal blocks */
    for (j = 1; j < l; ++j) {
      c = *col++;
      n = w + c;
     
      y[0] = p[c];
      y[1] = p[c+1];
      y[2] = p[c+2];

      m[0] += dat[0]*y[0] + dat[1]*y[1] + dat[2]*y[2];
      m[1] += dat[3]*y[0] + dat[4]*y[1] + dat[5]*y[2];
      m[2] += dat[6]*y[0] + dat[7]*y[1] + dat[8]*y[2];

      n[0] += dat[0]*x[0] + dat[3]*x[1] + dat[6]*x[2];
      n[1] += dat[1]*x[0] + dat[4]*x[1] + dat[7]*x[2];
      n[2] += dat[2]*x[0] + dat[5]*x[1] + dat[8]*x[2];
      dat += 9;
    }

    /* Add modifiation */
    o[0] += m[0];
    o[1] += m[1];
    o[2] += m[2];
  }
  return;
}

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

static void negate(double *r, const double *g, const int nn)
{
  int i;

  for (i = 0; i < nn; ++i) {
    r[0] = -g[0];
    r[1] = -g[1];
    r[2] = -g[2];
    r += 3;
    g += 3;
  }
  return;
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

static void maxpy(double *r, const double *x, const double c, const double *y, 
		  const int nn)
{
  int i;

  for (i = 0; i < nn; ++i) {
    r[0] = c*y[0] - x[0];
    r[1] = c*y[1] - x[1];
    r[2] = c*y[2] - x[2];
    r += 3;
    x += 3;
    y += 3;
  }
  return;
}

int optMesh(Mesh *mesh, int max_iter, double cn_tol, int precond)
{
  const double cg_tol = 1e-2;
  const double sigma = 1e-4;
  const double beta0 = 0.25;
  const double beta1 = 0.80;
  const double beta_tol = 1e-8;
  const double eta_1 = 0.01;
  const double eta_2 = 0.90;
  const double tr_incr = 10;
  const double tr_decr_def = 0.25;
  const double tr_decr_undef = 0.25;
  const double tr_num_tol = 1e-6;
  const int max_cg_iter = 200;
  
  double radius = 1000;		/* delta*delta */

  Precond  *prec;

  double *v;
  double *w;
  double *z;
  double *d;
  double *p;
  double *r;
  double *t;

  double norm_r, norm_g;
  double alpha, beta, kappa;
  double rz, rzm1;
  double dMp, norm_d, norm_dp1, norm_p;
  double obj, objn;

  const int nn = mesh->nn;
  int iter = 0, cg_iter;
  int nf = 1, ng = 1, nh = 0;
  int needLS;
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

  prec = preCreate(precond, nn, mesh->nz);

  v  = (double *)malloc(3*sizeof(double)*nn);
  w  = (double *)malloc(3*sizeof(double)*nn);
  z  = (double *)malloc(3*sizeof(double)*nn);
  d  = (double *)malloc(3*sizeof(double)*nn);
  p  = (double *)malloc(3*sizeof(double)*nn);
  r  = (double *)malloc(3*sizeof(double)*nn);

  gatherMesh(v, mesh);
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
  printf("%4d I      %10.9e %5.4e            %5.4e %5d "
         "usr: %5.4e sys: %5.4e tot: %5.4e\n",
         iter, obj, norm_g, fmax, fidx, t1, t2, t1+t2);
#else
  printf("%4d I      %10.9e %5.4e            "
         "usr: %5.4e sys: %5.4e tot: %5.4e\n",
         iter, obj, norm_g, t1, t2, t1+t2);
#endif

  while ((norm_g > cn_tol) && (iter < max_iter) && (radius > 1e-20)) {
    getrusage(RUSAGE_SELF, &r0);
    ++iter;

    ++nh;
    hOnly(mesh);
    prec->calc(prec, mesh);

    memset(d, 0, 3*sizeof(double)*nn);
    memcpy(r, mesh->g, 3*sizeof(double)*nn);
    norm_g *= cg_tol;

    prec->apply(z, r, prec, mesh);
    negate(p, z, nn);
    rz = inner(r, z, nn);

    dMp    = 0;
    norm_p = rz;
    norm_d = 0;

    cg_iter = 0;
    while ((sqrt(norm_r) > norm_g) && (cg_iter < max_cg_iter)) {
      ++cg_iter;

      memset(w, 0, 3*sizeof(double)*nn);
      matmul(w, mesh, p);

      kappa = inner(p, w, nn);
      if (kappa <= 0.0) {
        printf("Negative curvature.\n");
        alpha = (sqrt(dMp*dMp+norm_p*(radius-norm_d))-dMp)/norm_p;
        axpy(d, d, alpha, p, nn);
	break;
      }

      alpha = rz / kappa;

      norm_dp1 = norm_d + 2.0*alpha*dMp + alpha*alpha*norm_p;
      if (norm_dp1 >= radius) {
        printf("Radius Boundary.\n");
        alpha = (sqrt(dMp*dMp+norm_p*(radius-norm_d))-dMp)/norm_p;
        axpy(d, d, alpha, p, nn);
	break;
      }

      axpy(d, d, alpha, p, nn);
      axpy(r, r, alpha, w, nn);
      norm_r = norm(r, nn);

      prec->apply(z, r, prec, mesh);

      rzm1 = rz;
      rz = inner(r, z, nn);
      beta = rz / rzm1;
      maxpy(p, z, beta, p, nn);

      dMp = beta*(dMp + alpha*norm_p);
      norm_p = rz + beta*beta*norm_p;
      norm_d = norm_dp1;
    }

    alpha = inner(mesh->g, d, nn);

    memset(p, 0, 3*sizeof(double)*nn);
    matmul(p, mesh, d);
    beta = 0.5*inner(p, d, nn);
    kappa = alpha + beta;

    /* Put the new point into the locations */
    axpy(w, v, 1.0, d, nn);
    scatterMesh(mesh, w);

    ++nf;
    ++ng;
    ferr = gFcn(&objn, mesh);

#ifdef OUTPUT_STATS
    printf("alpha = %20.19e\n", alpha);
    printf("beta  = %20.19e\n", beta);
    printf("kappa = %20.19e\n", kappa);

    if (!ferr) {
      printf("new   = %20.19e\n", objn);
      printf("old   = %20.19e\n", obj);
      printf("ared  = %20.19e\n", objn - obj);
      printf("pred  = %20.19e\n", kappa);
      printf("ratio = %20.19e\n", (objn - obj) / kappa);
    }
#endif

    needLS = 0;
    beta = 1.0;

    if (!ferr) {
      if ((fabs(kappa) <= tr_num_tol) && (fabs(objn - obj) <= tr_num_tol)) {
	kappa = 1;
      }
      else {
	kappa = (objn - obj) / kappa;
      }
      
      if (kappa >= eta_1) {
	/* Iterate is acceptable */

        if (kappa >= eta_2) {
	  /* Iterate is a very good step, increase radius */
	  radius *= tr_incr;
	  if (radius > 1e20) {
	    radius = 1e20;
	  }
	}
      }
      else {
	/* Iterate is unacceptable */
	radius *= tr_decr_def;
	alpha *= sigma;
	beta *= beta1;
	needLS = 1;
      }
    }
    else {
      /* Function not defined at trial point */
      radius *= tr_decr_undef;
      alpha *= sigma;
      beta *= beta0;
      needLS = 1;
    }

    if (needLS) {
      /* Iterate is unacceptable, linesearch and reset radius */
      while (beta >= beta_tol) {
        axpy(w, v, beta, d, nn);
        scatterMesh(mesh, w);

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

      if (beta < beta_tol) {
        printf("Newton step not good\n");
        exit(-1);
      }

      /* Need a gradient evaluation for termination test */
      ++ng;
      gOnly(mesh);
    }

    /* Update the iterate (v = current point, w = new point) so swap */
    t = v;
    v = w;
    w = t;
    
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
    printf("%4d N %4d %10.9e %5.4e %5.4e %5.4e %5.4e %5d "
           "usr: %5.4e sys: %5.4e tot: %5.4e\n",
           iter, cg_iter, obj, norm_g, beta, radius, fmax, fidx, t1, t2, t1+t2);
#else
    printf("%4d N %4d %10.9e %5.4e %5.4e %5.4e usr: %5.4e sys: %5.4e tot: %5.4e\n",
           iter, cg_iter, obj, norm_g, beta, radius, t1, t2, t1+t2);
#endif

  }

  printf("Function evals: %4d\nGradient evals: %4d\nHessian  evals: %4d\n", 
	 nf, ng, nh);

  prec->destroy(prec);

  free(v);
  free(w);
  free(z);
  free(d);
  free(p);
  free(r);
  return 0;
}
