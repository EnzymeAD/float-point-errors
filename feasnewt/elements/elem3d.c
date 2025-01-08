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

  double  x[12];
  double  o_am = 0.0, o_gm = 0.0, o_hm = 0.0, o_co = 0.0, o_sq = 0.0, f;
  double  o_am_min1 = 0.0, o_am_min2 = 0.0;
  double  o_am_max1 = 0.0, o_am_max2 = 0.0;

  int     v1, v2, v3, v4;
  int     i;
  
  int     o_am_min1_idx = -1, o_am_min2_idx = -1;
  int     o_am_max1_idx = -1, o_am_max2_idx = -1;
  int     small = 0, medium = 0, large = 0;

  for (i = 0; i < e; ++i) {
    v1 = t[0];
    v2 = t[1];
    v3 = t[2];
    v4 = t[3];
    t += 4;

    w = v + 3*v1;
    x[0] = w[0];
    x[4] = w[1];
    x[8] = w[2];

    w = v + 3*v2;
    x[1] = w[0];
    x[5] = w[1];
    x[9] = w[2];

    w = v + 3*v3;
    x[2] = w[0];
    x[6] = w[1];
    x[10]= w[2];

    w = v + 3*v4;
    x[3] = w[0];
    x[7] = w[1];
    x[11]= w[2];

    v1 = p[v1];
    v2 = p[v2];
    v3 = p[v3];
    v4 = p[v4];

    if (o_fcn(&f, x)) return 1;
    o_am += f;
    o_sq += f*f;

    if ((v1 >= 0) && (v2 >= 0) && (v3 >= 0) && (v4 >= 0)) {
      if ((o_am_min1_idx < 0) || (f < o_am_min1)) {
        o_am_min1 = f;
        o_am_min1_idx = i;
      }

      if ((o_am_max1_idx < 0) || (f > o_am_max1)) {
        o_am_max1 = f;
        o_am_max1_idx = i;
      }
    }

    if ((v1 >= 0) || (v2 >= 0) || (v3 >= 0) || (v4 >= 0)) {
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

  printf("Sum of Squ: %5.4e\n", o_sq);
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

int histo(double *obj, const Mesh *m)
{
  const int e = m->ne;

  double *v = m->v;
  double *w;
  int    *t = m->e;

  double  x[12];
  double  o, f;
  int     v1, v2, v3, v4;
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
    v4 = t[3];
    t += 4;

    w = v + 3*v1;
    x[0] = w[0];
    x[4] = w[1];
    x[8] = w[2];

    w = v + 3*v2;
    x[1] = w[0];
    x[5] = w[1];
    x[9] = w[2];

    w = v + 3*v3;
    x[2] = w[0];
    x[6] = w[1];
    x[10]= w[2];

    w = v + 3*v4;
    x[3] = w[0];
    x[7] = w[1];
    x[11]= w[2];

#ifndef USE_WEIGHT
    if (o_fcn(&f, x)) return 1;
#else
    if (o_fcn(&f, x, weight, weight+6)) return 1;
    weight += 16;
#endif

    printf("%10.9e\n", f);

    o += f;
  }

  *obj = o;
  return 0;
}

int oFcn(double *obj, const Mesh *m)
{
  const int e = m->ne;

  double *v = m->v;
  double *w;
  int    *t = m->e;

  double  x[12];
  double  o, f;
  int     v1, v2, v3, v4;
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
    v4 = t[3];
    t += 4;

    w = v + 3*v1;
    x[0] = w[0];
    x[4] = w[1];
    x[8] = w[2];

    w = v + 3*v2;
    x[1] = w[0];
    x[5] = w[1];
    x[9] = w[2];

    w = v + 3*v3;
    x[2] = w[0];
    x[6] = w[1];
    x[10]= w[2];

    w = v + 3*v4;
    x[3] = w[0];
    x[7] = w[1];
    x[11]= w[2];

#ifndef USE_WEIGHT
    if (o_fcn(&f, x)) return 1;
#else
    if (o_fcn(&f, x, weight, weight+6)) return 1;
    weight += 16;
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

  double  x[12];
  double  f;

  int     v1, v2, v3, v4;
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
    v4 = t[3];
    t += 4;

    w = v + 3*v1;
    x[0] = w[0];
    x[4] = w[1];
    x[8] = w[2];

    w = v + 3*v2;
    x[1] = w[0];
    x[5] = w[1];
    x[9] = w[2];

    w = v + 3*v3;
    x[2] = w[0];
    x[6] = w[1];
    x[10]= w[2];

    w = v + 3*v4;
    x[3] = w[0];
    x[7] = w[1];
    x[11]= w[2];

#ifndef USE_WEIGHT
    if (o_fcn(&f, x)) return 1;
#else
    if (o_fcn(&f, x, weight, weight+6)) return 1;
    weight += 16;
#endif

    v1 = p[v1];
    v2 = p[v2];
    v3 = p[v3];
    v4 = p[v4];

    if (v1 >= 0) {
      ++cnt;
    }

    if (v2 >= 0) {
      ++cnt;
    }

    if (v3 >= 0) {
      ++cnt;
    }

    if (v4 >= 0) {
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

  double  x[12];
  double  d[12];
  double  o, f;

  int     v1, v2, v3, v4;
  int     i;

#ifdef USE_WEIGHT
  double *weight = m->w;
#endif

  *obj = 0.0;
  memset(g_obj, 0, 3*sizeof(double)*m->nn);

  o = 0.0;
  for (i = 0; i < e; ++i) {
    v1 = t[0];
    v2 = t[1];
    v3 = t[2];
    v4 = t[3];
    t += 4;

    w = v + 3*v1;
    x[0] = w[0];
    x[4] = w[1];
    x[8] = w[2];

    w = v + 3*v2;
    x[1] = w[0];
    x[5] = w[1];
    x[9] = w[2];

    w = v + 3*v3;
    x[2] = w[0];
    x[6] = w[1];
    x[10]= w[2];

    w = v + 3*v4;
    x[3] = w[0];
    x[7] = w[1];
    x[11]= w[2];

#ifndef USE_WEIGHT
    if (g_fcn(&f, d, x)) return 1;
#else
    if (g_fcn(&f, d, x, weight, weight+6)) return 1;
    weight += 16;
#endif

    v1 = p[v1];
    v2 = p[v2];
    v3 = p[v3];
    v4 = p[v4];

    if (v1 >= 0) {
      w = g_obj + v1;
      w[0] += d[0];
      w[1] += d[4];
      w[2] += d[8];
    }

    if (v2 >= 0) {
      w = g_obj + v2;
      w[0] += d[1];
      w[1] += d[5];
      w[2] += d[9];
    }

    if (v3 >= 0) {
      w = g_obj + v3;
      w[0] += d[2];
      w[1] += d[6];
      w[2] += d[10];
    }

    if (v4 >= 0) {
      w = g_obj + v4;
      w[0] += d[3];
      w[1] += d[7];
      w[2] += d[11];
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

  double  x[12];
  double  d[12];
  double  f;

  int     v1, v2, v3, v4;
  int     i;

#ifdef USE_WEIGHT
  double *weight = m->w;
#endif

  memset(g_obj, 0, 3*sizeof(double)*m->nn);

  for (i = 0; i < e; ++i) {
    v1 = t[0];
    v2 = t[1];
    v3 = t[2];
    v4 = t[3];
    t += 4;

    w = v + 3*v1;
    x[0] = w[0];
    x[4] = w[1];
    x[8] = w[2];

    w = v + 3*v2;
    x[1] = w[0];
    x[5] = w[1];
    x[9] = w[2];

    w = v + 3*v3;
    x[2] = w[0];
    x[6] = w[1];
    x[10]= w[2];

    w = v + 3*v4;
    x[3] = w[0];
    x[7] = w[1];
    x[11]= w[2];

#ifndef USE_WEIGHT
    g_fcn(&f, d, x);
#else
    g_fcn(&f, d, x, weight, weight+6);
    weight += 16;
#endif

    v1 = p[v1];
    v2 = p[v2];
    v3 = p[v3];
    v4 = p[v4];

    if (v1 >= 0) {
      w = g_obj + v1;
      w[0] += d[0];
      w[1] += d[4];
      w[2] += d[8];
    }

    if (v2 >= 0) {
      w = g_obj + v2;
      w[0] += d[1];
      w[1] += d[5];
      w[2] += d[9];
    }

    if (v3 >= 0) {
      w = g_obj + v3;
      w[0] += d[2];
      w[1] += d[6];
      w[2] += d[10];
    }

    if (v4 >= 0) {
      w = g_obj + v4;
      w[0] += d[3];
      w[1] += d[7];
      w[2] += d[11];
    }
  }
  return;
}

void gMesh(Mesh *m)
{
  const int  n  = m->nv;
  const int  nn = m->nn;

  int *p  = m->p;
  int *ip = m->i;

  int  i;

  if ((m->dat != NULL) || (m->g != NULL)) {
    /* Already calculated the structure.  Just return.                       */
    return;
  }

  m->g = (double *)malloc(3*sizeof(double)*nn);

  for (i = 0; i < n; ++i) {
    p[i] *= 3;
  }

  for (i = 0; i < nn; ++i) {
    ip[i] *= 3;
  }
  return;
}

double gNorm(const Mesh *m)
{
  const int nn = m->nn;

  double *g = m->g;
  double  norm_r = 0;
  int     i;

  for (i = 0; i < nn; ++i) {
    norm_r += g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
    g += 3;
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
  res[3] += src[3];
  res[4] += src[4];
  res[5] += src[5];
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

  double  x[12];
  double  d[12];
  double  h[24];
  double  A[6];
  double  o, f;

  int     v1, v2, v3, v4;
  int     i;

#ifdef USE_WEIGHT
  double *weight = m->w;
#endif

  *obj = 0.0;
  memset(g_obj, 0, 3*sizeof(double)*m->nn);
  memset(h_obj, 0, 6*sizeof(double)*m->nn);

  o = 0.0;
  for (i = 0; i < e; ++i) {
    v1 = t[0];
    v2 = t[1];
    v3 = t[2];
    v4 = t[3];
    t += 4;

    w = v + 3*v1;
    x[0] = w[0];
    x[4] = w[1];
    x[8] = w[2];

    w = v + 3*v2;
    x[1] = w[0];
    x[5] = w[1];
    x[9] = w[2];

    w = v + 3*v3;
    x[2] = w[0];
    x[6] = w[1];
    x[10]= w[2];

    w = v + 3*v4;
    x[3] = w[0];
    x[7] = w[1];
    x[11]= w[2];

#ifndef USE_WEIGHT
    if (h_fcn(&f, d, h, x)) return 1;
#else
    if (h_fcn(&f, d, h, x, weight, weight+6)) return 1;
    weight += 16;
#endif

    v1 = p[v1];
    v2 = p[v2];
    v3 = p[v3];
    v4 = p[v4];

    if (v1 >= 0) {
      /* Add the gradient information */
      w = g_obj + v1;
      w[0] += d[0];
      w[1] += d[4];
      w[2] += d[8];

      /* Add diagonal block */

      A[0] = h[0];
      A[1] = h[4];
      A[2] = h[8];
      A[3] = h[12];
      A[4] = h[16];
      A[5] = h[20];

      add_diag(h_obj + 2*v1, A);
    }

    if (v2 >= 0) {
      /* Add the gradient information */
      w = g_obj + v2;
      w[0] += d[1];
      w[1] += d[5];
      w[2] += d[9];

      /* Add diagonal block */

      A[0] = h[1];
      A[1] = h[5];
      A[2] = h[9];
      A[3] = h[13];
      A[4] = h[17];
      A[5] = h[21];

      add_diag(h_obj + 2*v2, A);
    }

    if (v3 >= 0) {
      /* Add the gradient information */
      w = g_obj + v3;
      w[0] += d[2];
      w[1] += d[6];
      w[2] += d[10];

      /* Add diagonal block */

      A[0] = h[2];
      A[1] = h[6];
      A[2] = h[10];
      A[3] = h[14];
      A[4] = h[18];
      A[5] = h[22];

      add_diag(h_obj + 2*v3, A);
    }

    if (v4 >= 0) {
      /* Add the gradient information */
      w = g_obj + v4;
      w[0] += d[3];
      w[1] += d[7];
      w[2] += d[11];

      /* Add diagonal block */

      A[0] = h[3];
      A[1] = h[7];
      A[2] = h[11];
      A[3] = h[15];
      A[4] = h[19];
      A[5] = h[23];

      add_diag(h_obj + 2*v4, A);
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

  double  x[12];
  double  h[24];
  double  A[6];

  int     v1, v2, v3, v4;
  int     i;

#ifdef USE_WEIGHT
  double *weight = m->w;
#endif

  memset(h_obj, 0, 6*sizeof(double)*m->nn);

  for (i = 0; i < e; ++i) {
    v1 = t[0];
    v2 = t[1];
    v3 = t[2];
    v4 = t[3];
    t += 4;

    w = v + 3*v1;
    x[0] = w[0];
    x[4] = w[1];
    x[8] = w[2];

    w = v + 3*v2;
    x[1] = w[0];
    x[5] = w[1];
    x[9] = w[2];

    w = v + 3*v3;
    x[2] = w[0];
    x[6] = w[1];
    x[10]= w[2];

    w = v + 3*v4;
    x[3] = w[0];
    x[7] = w[1];
    x[11]= w[2];

#ifndef USE_WEIGHT
    h_only(h, x);
#else
    h_only(h, x, weight, weight+6);
    weight += 16;
#endif

    v1 = p[v1];
    v2 = p[v2];
    v3 = p[v3];
    v4 = p[v4];

    if (v1 >= 0) {
      /* Add diagonal block */
      A[0] = h[0];
      A[1] = h[4];
      A[2] = h[8];
      A[3] = h[12];
      A[4] = h[16];
      A[5] = h[20];

      add_diag(h_obj + 2*v1, A);
    }

    if (v2 >= 0) {
      /* Add diagonal block */
      A[0] = h[1];
      A[1] = h[5];
      A[2] = h[9];
      A[3] = h[13];
      A[4] = h[17];
      A[5] = h[21];

      add_diag(h_obj + 2*v2, A);
    }

    if (v3 >= 0) {
      /* Add diagonal block */
      A[0] = h[2];
      A[1] = h[6];
      A[2] = h[10];
      A[3] = h[14];
      A[4] = h[18];
      A[5] = h[22];

      add_diag(h_obj + 2*v3, A);
    }

    if (v4 >= 0) {
      /* Add diagonal block */
      A[0] = h[3];
      A[1] = h[7];
      A[2] = h[11];
      A[3] = h[15];
      A[4] = h[19];
      A[5] = h[23];

      add_diag(h_obj + 2*v4, A);
    }
  }
  return;
}

void hMesh(Mesh *m)
{
  const int  n  = m->nv;
  const int  nn = m->nn;

  int *p  = m->p;
  int *ip = m->i;

  int  i;

  if ((m->dat != NULL) || (m->g != NULL)) {
    /* Already calculated the structure.  Just return.                       */
    return;
  }

  m->g = (double *)malloc(3*sizeof(double)*nn);
  m->dat = (double *)malloc(6*sizeof(double)*nn);

  for (i = 0; i < n; ++i) {
    p[i] *= 3;
  }

  for (i = 0; i < nn; ++i) {
    ip[i] *= 3;
  }
  return;
}

