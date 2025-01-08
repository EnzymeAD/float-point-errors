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
/* Preconditioner calculations for a single element.                         */
/*   Cannot use full space preconditioner.                                   */
/*****************************************************************************/

#if 0
static void preCalc0(double p[3], const double h[3])
{
  return;
}

static void preApply0(double z[2], const double v[2], const double p[3])
{
  z[0] = v[0];
  z[1] = v[1];
  return;
}

static void preCalc1(double p[3], const double h[3])
{
  p[0] = 1.0 / (h[0] + h[2]);
  return;
}

static void preApply1(double z[2], const double v[2], const double p[3])
{
  z[0] = v[0]*p[0];
  z[1] = v[1]*p[0];
  return;
}
#endif

static void preCalc2(double p[3], const double h[3])
{
  p[0] = 1.0  / h[0];
  p[1] = h[1] * p[0];

  p[2] = 1.0  / (h[2]-h[1]*p[1]);
  return;
}

static void preApply2(double z[2], const double v[2], const double p[3])
{
  z[0] = v[0];
  z[1] = v[1] - p[1]*z[0];

  z[0] *= p[0];
  z[1] *= p[2];

  z[0] -= p[1]*z[1];
  return;
}

static double norm(const double *g, const int nn)
{
  double norm_r = 0;
  int    i;
     
  for (i = 0; i < nn; ++i) {
    norm_r += g[0]*g[0] + g[1]*g[1];
    g += 2;
  }
  return norm_r;
}

int optMesh(Mesh *mesh, int max_iter, double cn_tol, int precond)
{
  const double sigma  = 1e-4;
  const double beta0  = 0.25;
  const double beta1  = 0.80;
  const double tol1   = 1e-8;
  const int max_min_iter = 1;

  const int n = mesh->nv;
  const int n1 = n + 1;
  const int e = mesh->ne;
  const int nn = mesh->nn;

  double *mv;
  int *sta;
  int *tri;
  int *t;

  double g_obj[2];
  double h_obj[3];
  double p_obj[3];

  double v[2];
  double d[2];
  double r[2];

  double conv;
  double norm_r, norm_g;
  double alpha, beta = 1.0;
  double obj, objn;

  int iter = 0, min_iter, cg_iter = 0;
  int nf = 1, ng = 1, nh = 0;
  int ferr;
  int i, j, k;

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

  /* Start off by initializing the gradient for the entire mesh.             */
  gMesh(mesh);
  if (gFcn(&obj, mesh)) {
    fprintf(stderr, "Invalid starting point.\n");
    exit(-1);
  }

  /* Now calculate the node->elements listing.                               */
  sta = (int *)calloc(n1+1, sizeof(int)); ++sta;
  tri = (int *)malloc(sizeof(int)*3*e);

  t = mesh->e;
  for (i = 0; i < e; ++i) {
    ++sta[t[0]];
    ++sta[t[1]];
    ++sta[t[2]];
    t += 3;
  }

  /* Calculate the starts */
  k = sta[0];
  sta[0] = 0;
  for (i = 1; i < n1; ++i) {
    j = sta[i];
    sta[i] = k;
    k += j;
  }

  t = mesh->e;
  for (i = 0; i < 3*e; i += 3) {
    tri[sta[t[0]]++] = i;
    tri[sta[t[1]]++] = i;
    tri[sta[t[2]]++] = i;
    t += 3;
  }

  /* Recover the starts */
  --sta;

  norm_r = norm(mesh->g, nn);
  conv   = sqrt(norm_r);

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
         iter, obj, conv, fmax, fidx, t1, t2, t1+t2);
#else
  printf("%4d I      %10.9e %5.4e            "
         "usr: %5.4e sys: %5.4e tot: %5.4e\n",
         iter, obj, conv, t1, t2, t1+t2);
#endif

  while ((conv > cn_tol) && (iter < max_iter)) {
    getrusage(RUSAGE_SELF, &r0);

    ++iter;

    ++nf;
    ++ng;
    ++nh;

    mv = mesh->v;
    for (i = 0; i < n; ++i) {
      if (mesh->p[i] < 0) {
	mv += 2;
	continue;
      }

      hFcnl(&obj, g_obj, h_obj, i, tri + sta[i], sta[i+1] - sta[i], mesh);

      norm_r = g_obj[0]*g_obj[0] + g_obj[1]*g_obj[1];
      norm_g = sqrt(norm_r);

      min_iter = 0;
      while ((norm_g > conv / 1000.0) && (min_iter < max_min_iter)) {
        ++min_iter;

	v[0] = mv[0];
	v[1] = mv[1];
	
        preCalc2(p_obj, h_obj);
	r[0] = -g_obj[0];
	r[1] = -g_obj[1];
        preApply2(d, r, p_obj);

	alpha = g_obj[0]*d[0] + g_obj[1]*d[1];
	if (alpha > 0.0) {
	  printf("Not descent.\n");
	  exit(-1);
	}
	else {
	  alpha *= sigma;
	}
	beta = 1.0;
	
	/* Unrolled for better performance.  Only do a gradient evaluation   */
	/* when beta = 1.0.  Do not do this at other times.                  */

	mv[0] = v[0] + beta*d[0];
	mv[1] = v[1] + beta*d[1];

	ferr = gFcnl(&objn, g_obj, i, tri + sta[i], sta[i+1] - sta[i], mesh);
	if ((!ferr && (obj - objn >= -alpha*beta - epsilon)) ||
	    (!ferr && (sqrt(g_obj[0]*g_obj[0] + 
			    g_obj[1]*g_obj[1]) < 100*cn_tol))) {
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
	    mv[0] = v[0] + beta*d[0];
	    mv[1] = v[1] + beta*d[1];

	    if (oFcnl(&objn, i, tri + sta[i], sta[i+1] - sta[i], mesh)) {
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

	  gFcnl(&objn, g_obj, i, tri + sta[i], sta[i+1] - sta[i], mesh);
	}

	obj = objn;
        norm_r = g_obj[0]*g_obj[0] + g_obj[1]*g_obj[1];
        norm_g = sqrt(norm_r);

        if ((norm_g > conv / 1000.0) && (min_iter < max_min_iter)) {
          hFcnl(&objn, g_obj, h_obj, i, tri + sta[i], sta[i+1] - sta[i], mesh);
        }
      }
      mv += 2;
    }

    ++nf;
    ++ng;
    if (gFcn(&obj, mesh)) {
      fprintf(stderr, "Invalid intermediate point.\n");
      exit(-1);
    }

    norm_r = norm(mesh->g, nn);
    conv   = sqrt(norm_r);

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
    printf("%4d N %4d %10.9e %5.4e %5.4e %5.4e %5d "
           "usr: %5.4e sys: %5.4e tot: %5.4e\n",
           iter, cg_iter, obj, conv, beta, fmax, fidx, t1, t2, t1+t2);
#else
    printf("%4d N %4d %10.9e %5.4e %5.4e usr: %5.4e sys: %5.4e tot: %5.4e\n",
           iter, cg_iter, obj, conv, beta, t1, t2, t1+t2);
#endif
  }

  printf("Function evals: %4d\nGradient evals: %4d\nHessian  evals: %4d\n", 
	 nf, ng, nh);

  free(sta);
  free(tri);
  return 0;
}
