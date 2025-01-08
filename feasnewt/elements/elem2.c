#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fcn.h"

#ifdef CHECK
int cFcn(const Mesh *m)
{
  const int e = m->ne;

  double *v = m->v;
  double *w;

  int    *t = m->e;
  int    *p = m->p;

  double  x[6];
  double  o_am = 0.0, o_gm = 0.0, o_hm = 0.0, o_co = 0.0, f;
  double  o_am_min1 = 0.0, o_am_min2 = 0.0;
  double  o_am_max1 = 0.0, o_am_max2 = 0.0;

  int     v1, v2, v3;
  int     i;

  int     o_am_min1_idx = -1, o_am_min2_idx = -1;
  int     o_am_max1_idx = -1, o_am_max2_idx = -1;
  int     small = 0, medium = 0, large = 0;

  for (i = 0; i < e; ++i) {
    v1 = t[0];
    v2 = t[1];
    v3 = t[2];
    t += 3;

    w = v + 2*v1;
    x[0] = w[0];
    x[3] = w[1];

    w = v + 2*v2;
    x[1] = w[0];
    x[4] = w[1];

    w = v + 2*v3;
    x[2] = w[0];
    x[5] = w[1];

    v1 = p[v1];
    v2 = p[v2];
    v3 = p[v3];

    if (o_fcn(&f, x)) return 1;
    o_am += f;

    if ((v1 >= 0) && (v2 >= 0) && (v3 >= 0)) {
      if ((o_am_min1_idx < 0) || (f < o_am_min1)) {
        o_am_min1 = f;
        o_am_min1_idx = i;
      }

      if ((o_am_max1_idx < 0) || (f > o_am_max1)) {
        o_am_max1 = f;
        o_am_max1_idx = i;
      }
    }

    if ((v1 >= 0) || (v2 >= 0) || (v3 >= 0)) {
      if ((o_am_min2_idx < 0) || (f < o_am_min2)) {
        o_am_min2 = f;
        o_am_min2_idx = i;
      }

      if ((o_am_max2_idx < 0) || (f > o_am_max2)) {
        o_am_max2 = f;
        o_am_max2_idx = i;
      }
    }

    if (f <= 1.25) {
      ++small;
    }
    else if (f <= 3.00) {
      ++medium;
    }
    else {
      ++large;
    }

    if (o_fcn_gm(&f, x)) return 1;
    o_gm += f;

    if (o_fcn_hm(&f, x)) return 1;
    o_hm += f;

    if (o_cond(&f, x)) return 1;
    o_co += f;
  }

  o_am /= e;
  o_gm  = exp(o_gm / e);
  o_hm /= e;
  o_co /= e;

  printf("Minimum AF: %5.4e (%d)\n", o_am_min1, o_am_min1_idx);
  printf("Maximum AF: %5.4e (%d)\n\n", o_am_max1, o_am_max1_idx);

  printf("Minimum OF: %5.4e (%d)\n", o_am_min2, o_am_min2_idx);
  printf("Maximum OF: %5.4e (%d)\n\n", o_am_max2, o_am_max2_idx);

  printf("Small     : %d\n", small);
  printf("Medium    : %d\n", medium);
  printf("Large     : %d\n\n", large);

  printf("Arithmetic: %5.4e\n", o_am);
  printf("Geometric : %5.4e\n", o_gm);
  printf("Harmonic  : %5.4e\n", o_hm);
  printf("Condition : %5.4e\n", o_co);
  return 0;
}
#endif

/*****************************************************************************/
/* The next set of functions calculate the function, gradient, and Hessian   */
/* for the entire mesh.  The coordinates are stores in the input mesh, which */
/* contains the values for the fixed vertices.  The output is in the         */
/* compacted representation.                                                 */
/*****************************************************************************/

int oFcn(double *obj, const Mesh *m)
{
  const int e = m->ne;

  double *v = m->v;
  double *w;
  int    *t = m->e;

  double  x[6];
  double  o, f;
  int     v1, v2, v3;
  int     i;
  
#ifdef USE_WEIGHT
  double *weight = m->w;
#endif

  *obj = 0.0;

  o = 0.0;
  for (i = 0; i < e; ++i) {
    v1 = t[0];
    v2 = t[1];
    v3 = t[2];
    t += 3;

    w = v + 2*v1;
    x[0] = w[0];
    x[3] = w[1];

    w = v + 2*v2;
    x[1] = w[0];
    x[4] = w[1];

    w = v + 2*v3;
    x[2] = w[0];
    x[5] = w[1];

#ifndef USE_WEIGHT
    if (o_fcn(&f, x)) return 1;
#else
    if (o_fcn(&f, x, weight, weight+3)) return 1;
    weight += 9;
#endif

    o += f;
  }

  *obj = o;
  return 0;
}

int oMax(double *max, int *idx, const Mesh *m)
{
  const int e = m->ne;

  double *v = m->v;
  double *w;

  int    *t = m->e;
  int    *p = m->p;

  double  x[6];
  double  f;
  int     v1, v2, v3;
  int     i;
  int     cnt;

#ifdef USE_WEIGHT
  double *weight = m->w;
#endif

  *max =  0.0;
  *idx = -1;

  for (i = 0; i < e; ++i) {
    cnt = 0;

    v1 = t[0];
    v2 = t[1];
    v3 = t[2];
    t += 3;

    w = v + 2*v1;
    x[0] = w[0];
    x[3] = w[1];

    w = v + 2*v2;
    x[1] = w[0];
    x[4] = w[1];

    w = v + 2*v3;
    x[2] = w[0];
    x[5] = w[1];

#ifndef USE_WEIGHT
    if (o_fcn(&f, x)) return 1;
#else
    if (o_fcn(&f, x, weight, weight+3)) return 1;
    weight += 9;
#endif

    v1 = p[v1];
    v2 = p[v2];
    v3 = p[v3];

    if (v1 >= 0) {
      ++cnt;
    }

    if (v2 >= 0) {
      ++cnt;
    }

    if (v3 >= 0) {
      ++cnt;
    }

    if (cnt && (f > *max)) {
      *max = f;
      *idx = i;
    }
  }
  return 0;
}

int gFcn(double *obj, const Mesh *m)
{
  const int e = m->ne;

  double *v = m->v;
  double *g_obj = m->g;
  double *w;

  int    *t = m->e;
  int    *p = m->p;

  double  x[6];
  double  d[6];
  double  o, f;

  int     v1, v2, v3;
  int     i;

#ifdef USE_WEIGHT
  double *weight = m->w;
#endif

  *obj = 0.0;
  memset(g_obj, 0, 2*sizeof(double)*m->nn);

  o = 0.0;
  for (i = 0; i < e; ++i) {
    v1 = t[0];
    v2 = t[1];
    v3 = t[2];
    t += 3;

    w = v + 2*v1;
    x[0] = w[0];
    x[3] = w[1];

    w = v + 2*v2;
    x[1] = w[0];
    x[4] = w[1];

    w = v + 2*v3;
    x[2] = w[0];
    x[5] = w[1];

#ifndef USE_WEIGHT
    if (g_fcn(&f, d, x)) return 1;
#else
    if (g_fcn(&f, d, x, weight, weight+3)) return 1;
    weight += 9;
#endif

    v1 = p[v1];
    v2 = p[v2];
    v3 = p[v3];

    if (v1 >= 0) {
      w = g_obj + v1;
      w[0] += d[0];
      w[1] += d[3];
    }

    if (v2 >= 0) {
      w = g_obj + v2;
      w[0] += d[1];
      w[1] += d[4];
    }

    if (v3 >= 0) {
      w = g_obj + v3;
      w[0] += d[2];
      w[1] += d[5];
    }

    o += f;
  }

  *obj = o;
  return 0;
}

void gOnly(const Mesh *m)
{
  const int e = m->ne;

  double *v = m->v;
  double *g_obj = m->g;
  double *w;

  int    *t = m->e;
  int    *p = m->p;

  double  x[6];
  double  d[6];
  double  f;
  int     v1, v2, v3;
  int     i;

#ifdef USE_WEIGHT
  double *weight = m->w;
#endif

  memset(g_obj, 0, 2*sizeof(double)*m->nn);

  for (i = 0; i < e; ++i) {
    v1 = t[0];
    v2 = t[1];
    v3 = t[2];
    t += 3;

    w = v + 2*v1;
    x[0] = w[0];
    x[3] = w[1];

    w = v + 2*v2;
    x[1] = w[0];
    x[4] = w[1];

    w = v + 2*v3;
    x[2] = w[0];
    x[5] = w[1];

#ifndef USE_WEIGHT
    g_fcn(&f, d, x);
#else
    g_fcn(&f, d, x, weight, weight+3);
    weight += 9;
#endif

    v1 = p[v1];
    v2 = p[v2];
    v3 = p[v3];

    if (v1 >= 0) {
      w = g_obj + v1;
      w[0] += d[0];
      w[1] += d[3];
    }

    if (v2 >= 0) {
      w = g_obj + v2;
      w[0] += d[1];
      w[1] += d[4];
    }

    if (v3 >= 0) {
      w = g_obj + v3;
      w[0] += d[2];
      w[1] += d[5];
    }
  }
  return;
}

void gMesh(Mesh *m)
{
  int *p  = m->p;
  int *ip = m->i;
  int  n  = m->nv;
  int  nn = m->nn;

  int  i;

  if ((m->dat != NULL) || (m->g != NULL)) {
    /* Already calculated the structure.  Just return.                       */
    return;
  }

  m->g = (double *)malloc(2*sizeof(double)*nn);

  for (i = 0; i < n; ++i) {
    *p++ *= 2;
  }

  for (i = 0; i < nn; ++i) {
    *ip++ *= 2;
  }
  return;
}

double gNorm(const Mesh *m)
{
  double *g = m->g;
  double  norm_r = 0;
  int     i, nn = m->nn;

  for (i = 0; i < nn; ++i) {
    norm_r += g[0]*g[0] + g[1]*g[1];
    g += 2;
  }
  return sqrt(norm_r);
}

/*****************************************************************************/
/* The following three functions are for dealing with blocks.  These are     */
/* used in the accumulation of the upper triangular part of the Hessian      */
/* matrix.  We differentiate between diagonal elements and off-diagonal      */
/* elements.  For the off diagonal ones, we either add the block normally    */
/* or the transpose depending on the ordering in the element.                */
/*****************************************************************************/

static void add_diag(double *res, const double *src)
{
  res[0] += src[0];
  res[1] += src[1];
  res[2] += src[2];
  return;
}

static void add_block(double *res, const double *src)
{
  res[0] += src[0];
  res[1] += src[1];
  res[2] += src[2];
  res[3] += src[3];
  return;
}

static void add_blockT(double *res, const double *src)
{
  res[0] += src[0];
  res[1] += src[2];
  res[2] += src[1];
  res[3] += src[3];
  return;
}

int hFcn(double *obj, const Mesh *m)
{
  const int e = m->ne;

  double *v = m->v;
  double *g_obj = m->g;
  double *h_obj = m->dat;
  double *w;

  int    *t = m->e;
  int    *p = m->p;
  int    *inst  = m->inst;

  double  x[6];
  double  d[6];
  double  h[21];
  double  A[4];
  double  o, f;

  int     v1, v2, v3;
  int     i;

#ifdef USE_WEIGHT
  double *weight = m->w;
#endif

  *obj = 0.0;
  memset(g_obj, 0, 2*sizeof(double)*m->nn);
  memset(h_obj, 0,   sizeof(double)*m->nz);

  o = 0.0;
  for (i = 0; i < e; ++i) {
    v1 = t[0];
    v2 = t[1];
    v3 = t[2];
    t += 3;

    w = v + 2*v1;
    x[0] = w[0];
    x[3] = w[1];

    w = v + 2*v2;
    x[1] = w[0];
    x[4] = w[1];

    w = v + 2*v3;
    x[2] = w[0];
    x[5] = w[1];

#ifndef USE_WEIGHT
    if (h_fcn(&f, d, h, x)) return 1;
#else
    if (h_fcn(&f, d, h, x, weight, weight+3)) return 1;
    weight += 9;
#endif

    v1 = p[v1];
    v2 = p[v2];
    v3 = p[v3];

    if (v1 >= 0) {
      /* Add the gradient information */
      w = g_obj + v1;
      w[0] += d[0];
      w[1] += d[3];

      /* Add diagonal block */

      A[0] = h[0];
      A[1] = h[6];
      A[2] = h[15];

      add_diag(h_obj + *inst++, A);

      if (v2 >= 0) {

	A[0] = h[1];
	A[1] = h[7];
  	A[2] = h[9];
	A[3] = h[16];

	/* Add first off diagonal block */
	if (v1 < v2) {
	  add_block(h_obj + *inst++, A);
        }
        else {
	  add_blockT(h_obj + *inst++, A);
        }
      }

      if (v3 >= 0) {

	A[0] = h[2];
	A[1] = h[8];
  	A[2] = h[12];
	A[3] = h[17];

	/* Add second off diagonal block */
	if (v1 < v3) {
	  add_block(h_obj + *inst++, A);
        }
        else {
	  add_blockT(h_obj + *inst++, A);
        }
      }
    }

    if (v2 >= 0) {
      /* Add the gradient information */
      w = g_obj + v2;
      w[0] += d[1];
      w[1] += d[4];

      /* Add diagonal block */

      A[0] = h[3];
      A[1] = h[10];
      A[2] = h[18];

      add_diag(h_obj + *inst++, A);

      if (v3 >= 0) {

	A[0] = h[4];
	A[1] = h[11];
  	A[2] = h[13];
	A[3] = h[19];

	/* Add first off diagonal block */
	if (v2 < v3) {
	  add_block(h_obj + *inst++, A);
        }
        else {
	  add_blockT(h_obj + *inst++, A);
        }
      }
    }

    if (v3 >= 0) {
      /* Add the gradient information */
      w = g_obj + v3;
      w[0] += d[2];
      w[1] += d[5];

      /* Add diagonal block */

      A[0] = h[5];
      A[1] = h[14];
      A[2] = h[20];

      add_diag(h_obj + *inst++, A);
    }

    o += f;
  }

  *obj = o;
  return 0;
}

void hOnly(const Mesh *m)
{
  const int e = m->ne;

  double *v = m->v;
  double *h_obj = m->dat;
  double *w;

  int    *t = m->e;
  int    *p = m->p;
  int    *inst  = m->inst;

  double  x[6];
  double  h[21];
  double  A[4];

  int     v1, v2, v3;
  int     i;

#ifdef USE_WEIGHT
  double *weight = m->w;
#endif

  memset(h_obj, 0, sizeof(double)*m->nz);

  for (i = 0; i < e; ++i) {
    v1 = t[0];
    v2 = t[1];
    v3 = t[2];
    t += 3;

    w = v + 2*v1;
    x[0] = w[0];
    x[3] = w[1];

    w = v + 2*v2;
    x[1] = w[0];
    x[4] = w[1];

    w = v + 2*v3;
    x[2] = w[0];
    x[5] = w[1];

#ifndef USE_WEIGHT
    h_only(h, x);
#else
    h_only(h, x, weight, weight+3);
    weight += 9;
#endif

    v1 = p[v1];
    v2 = p[v2];
    v3 = p[v3];

    if (v1 >= 0) {
      /* Add diagonal block */
      A[0] = h[0];
      A[1] = h[6];
      A[2] = h[15];

      add_diag(h_obj + *inst++, A);

      if (v2 >= 0) {

	A[0] = h[1];
	A[1] = h[7];
  	A[2] = h[9];
	A[3] = h[16];

	/* Add first off diagonal block */
	if (v1 < v2) {
	  add_block(h_obj + *inst++, A);
        }
        else {
	  add_blockT(h_obj + *inst++, A);
        }
      }

      if (v3 >= 0) {

	A[0] = h[2];
	A[1] = h[8];
  	A[2] = h[12];
	A[3] = h[17];

	/* Add second off diagonal block */
	if (v1 < v3) {
	  add_block(h_obj + *inst++, A);
        }
        else {
	  add_blockT(h_obj + *inst++, A);
        }
      }
    }

    if (v2 >= 0) {
      /* Add diagonal block */
      A[0] = h[3];
      A[1] = h[10];
      A[2] = h[18];

      add_diag(h_obj + *inst++, A);

      if (v3 >= 0) {

	A[0] = h[4];
	A[1] = h[11];
  	A[2] = h[13];
	A[3] = h[19];

	/* Add first off diagonal block */
	if (v2 < v3) {
	  add_block(h_obj + *inst++, A);
        }
        else {
	  add_blockT(h_obj + *inst++, A);
        }
      }
    }

    if (v3 >= 0) {
      /* Add diagonal block */
      A[0] = h[5];
      A[1] = h[14];
      A[2] = h[20];

      add_diag(h_obj + *inst++, A);
    }
  }
  return;
}

/*****************************************************************************/
/* The next function calculates the structure of the hessian.  This needs to */
/* be rewritten to be better.                                                */
/*    Currently only needs nn rows.                                          */
/*****************************************************************************/

void hMesh(Mesh *m)
{
  const int  n  = m->nv;
  const int  nn = m->nn;
  const int  nn1 = nn + 1;
  const int  e = m->ne;

  int *t  = m->e;
  int *p  = m->p;
  int *ip = m->i;

  int *cl;
  int *cr;
  int *ci;

  int *rl;
  int *rc;
  int *ri;

  int  nb = m->nb;

  int  perm[3];
  int  i, j, k;
  int  r, c;
  int  st, en;
  int  ndb = 0, nob = 0;

  if ((m->dat != NULL) || (m->g != NULL)) {
    /* Already calculated the structure.  Just return.                       */
    return;
  }

  /* Calculate the structure of the hessian matrix.  This is done by         */
  /* blocks of coordinates.                                                  */

  cl = (int *)calloc(nn1+1, sizeof(int)); ++cl; /* Length of the column */
  cr = (int *)malloc(nb*sizeof(int));     /* Row index for the column   */
  ci = (int *)malloc(nb*sizeof(int));     /* Instruction for the column */
  
  rl = (int *)calloc(nn1+1, sizeof(int)); ++rl; /* Length of the row   */
  rc = (int *)malloc(nb*sizeof(int));     /* Row index for the row     */
  ri = (int *)malloc(nb*sizeof(int));     /* Instruction for the row   */
  
  /* Start by calculating a compressed sparse column representation of the   */
  /* matrix.  Begin by counting the number of elements in each column.       */
  for (i = 0; i < e; ++i) {
    perm[0]  = p[t[0]];
    perm[1]  = p[t[1]];
    perm[2]  = p[t[2]];
    t += 3;

    /* Make the resulting hessian upper triangular (slightly better locality */
    /* of reference)                                                         */

    for (j = 0; j < 3; ++j) {
      r = perm[j]; 
      if (r >= 0) { 
        for (k = j; k < 3; ++k) {
          c = perm[k];
	  if (c >= r) {
	    ++cl[c];
	  }
	  else if (c >= 0) {
	    ++cl[r];
	  }
        }
      }
    }
  }
  
  /* Calculate the column starts                                             */
  nb = cl[0];
  cl[0] = 0;
  for (i = 1; i < nn1; ++i) {
    j = cl[i];
    cl[i] = nb;
    nb += j;
  }

  nb = 0;
  t = m->e;
  for (i = 0; i < e; ++i) {
    perm[0]  = p[t[0]];
    perm[1]  = p[t[1]];
    perm[2]  = p[t[2]];
    t += 3;

    /* Make the resulting hessian upper triangular (slightly better locality */
    /* of reference)                                                         */

    for (j = 0; j < 3; ++j) {
      r = perm[j];
      if (r >= 0) {
        for (k = j; k < 3; ++k) {
          c = perm[k];
          if (c >= r) {
	    cr[cl[c]] = r;
	    ci[cl[c]++] = nb++;
	  }
	  else if (c >= 0) {
	    cr[cl[r]] = c;
	    ci[cl[r]++] = nb++;
	  }
        }
      }
    }
  }

  /* Recover the column starts                                               */
  --cl;

  /* Sort by row index                                                       */
  /* Calculate the lengths of the rows                                       */
  memset(rl, 0, nn1*sizeof(int));
  for (i = 0; i < nb; ++i) {
    rl[cr[i]]++;
  }

  /* Calculate the starts of the rows                                        */
  nb = rl[0];
  rl[0] = 0;
  for (i = 1; i < nn1; ++i) {
    j = rl[i];
    rl[i] = nb;
    nb += j;
  }

  /* Now do the sorting into an array                                        */
  for (i = 0; i < nn; ++i) {
    st = cl[i];
    en = cl[i+1];

    while (st < en) {
      k = rl[cr[st]]++;
      rc[k] = i;
      ri[k] = ci[st++];
    }
  }

  /* Recover the row starts                                                  */
  --rl;

  /* Compact (remove fixed, put diagonal first, ...) and patch instructions  */
  /* ii is the address where the element comes from                          */

  memset(cl, 0, nn1*sizeof(int));

  j = 0;
  k = 0;

  for (r = 0; r < nn; ++r) {
    st = rl[r];
    en = rl[r+1];

    while (st < en) {
      c = rc[st];

      while ((st < en) && (c == rc[st])) {
        ci[ri[st++]] = j;
      }

      cl[r]++;
      cr[k++] = 2*c;
      if (r != c) {
        j += 4;
	++nob;
      }
      else {
        j += 3;
	++ndb;
      }
    }
  }

  /* Deallocate unnecessary stuff */
  free(rl);
  free(rc);
  free(ri);

  m->nz = 3*ndb + 4*nob;

#if 0
  printf("final     diagonal blocks: %d\n", ndb);
  printf("final off diagonal blocks: %d\n", nob);
  printf("final nonzero count      : %d\n", m->nz);
#endif

  /* Fill in return values */
  m->len  = cl;
  m->col  = (int *)realloc(cr, (ndb+nob)*sizeof(int));
  m->dat  = (double *)malloc(m->nz*sizeof(double));
  m->inst = ci;

  m->g = (double *)malloc(2*sizeof(double)*nn);

  for (i = 0; i < n; ++i) {
    p[i] *= 2;
  }

  for (i = 0; i < nn; ++i) {
    ip[i] *= 2;
  }
  return;
}

