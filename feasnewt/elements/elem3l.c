#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fcn.h"

int oFcnl(double *obj,
          const int  vert,
	  const int *elems, const int nelems, const Mesh *m)
{
  int    *t = m->e;
  double *v = m->v;
  double *w;

  double  x[12];
  double  o, f;

  int     v1, v2, v3, v4;
  int     i, e;

#ifdef USE_WEIGHT
  double *weight;
#endif

  *obj = 0.0;

  o = 0.0;
  for (i = 0; i < nelems; ++i) {
    e = *elems++;
    v1 = t[e+0];
    v2 = t[e+1];
    v3 = t[e+2];
    v4 = t[e+3];

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
    weight = m->w + 4*e;  /* 16*e/4 */
    if (o_fcn(&f, x, weight, weight+6)) return 1;
#endif

    /* Add in objective function */
    o += f;
  }

  *obj = o;
  return 0;
}

int gFcnl(double *obj, double *g_obj,
          const int  vert,
	  const int *elems, const int nelems, const Mesh *m)
{
  int    *t = m->e;
  double *v = m->v;
  double *w;

  double  x[12];
  double  d[3];
  double  o, f;

  int     v1, v2, v3, v4;
  int     i, e;

#ifdef USE_WEIGHT
  double *weight;
#endif

  *obj = 0.0;

  g_obj[0] = 0.0;
  g_obj[1] = 0.0;
  g_obj[2] = 0.0;

  o = 0.0;
  for (i = 0; i < nelems; ++i) {
    e = *elems++;
    v1 = t[e+0];
    v2 = t[e+1];
    v3 = t[e+2];
    v4 = t[e+3];

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
    if (v1 == vert) {
      if (g_fcnl_0(&f, d, x)) return 1;
    }
    else if (v2 == vert) {
      if (g_fcnl_1(&f, d, x)) return 1;
    }
    else if (v3 == vert) {
      if (g_fcnl_2(&f, d, x)) return 1;
    }
    else if (v4 == vert) {
      if (g_fcnl_3(&f, d, x)) return 1;
    }
#else
    weight = m->w + 4*e;  /* 16*e/4 */
    if (v1 == vert) {
      if (g_fcnl_0(&f, d, x, weight, weight+6)) return 1;
    }
    else if (v2 == vert) {
      if (g_fcnl_1(&f, d, x, weight, weight+6)) return 1;
    }
    else if (v3 == vert) {
      if (g_fcnl_2(&f, d, x, weight, weight+6)) return 1;
    }
    else if (v4 == vert) {
      if (g_fcnl_3(&f, d, x, weight, weight+6)) return 1;
    }
#endif

    /* Add the gradient information */
    g_obj[0] += d[0];
    g_obj[1] += d[1];
    g_obj[2] += d[2];
    
    /* Add in objective function */
    o += f;
  }

  *obj = o;
  return 0;
}

int hFcnl(double *obj, double *g_obj, double *h_obj, 
          const int  vert,
	  const int *elems, const int nelems, const Mesh *m)
{
  int    *t = m->e;
  double *v = m->v;
  double *w;

  double  x[12];
  double  d[3];
  double  h[6];
  double  o, f;

  int     v1, v2, v3, v4;
  int     i, e;

#ifdef USE_WEIGHT
  double *weight;
#endif

  *obj = 0.0;

  g_obj[0] = 0.0;
  g_obj[1] = 0.0;
  g_obj[2] = 0.0;

  h_obj[0] = 0.0;
  h_obj[1] = 0.0;
  h_obj[2] = 0.0;
  h_obj[3] = 0.0;
  h_obj[4] = 0.0;
  h_obj[5] = 0.0;

  o = 0.0;
  for (i = 0; i < nelems; ++i) {
    e = *elems++;
    v1 = t[e+0];
    v2 = t[e+1];
    v3 = t[e+2];
    v4 = t[e+3];

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
    if (v1 == vert) {
      if (h_fcnl_0(&f, d, h, x)) return 1;
    }
    else if (v2 == vert) {
      if (h_fcnl_1(&f, d, h, x)) return 1;
    }
    else if (v3 == vert) {
      if (h_fcnl_2(&f, d, h, x)) return 1;
    }
    else if (v4 == vert) {
      if (h_fcnl_3(&f, d, h, x)) return 1;
    }
#else
    weight = m->w + 4*e;  /* 16*e/4 */
    if (v1 == vert) {
      if (h_fcnl_0(&f, d, h, x, weight, weight+6)) return 1;
    }
    else if (v2 == vert) {
      if (h_fcnl_1(&f, d, h, x, weight, weight+6)) return 1;
    }
    else if (v3 == vert) {
      if (h_fcnl_2(&f, d, h, x, weight, weight+6)) return 1;
    }
    else if (v4 == vert) {
      if (h_fcnl_3(&f, d, h, x, weight, weight+6)) return 1;
    }
#endif

    /* Add the gradient information */
    g_obj[0] += d[0];
    g_obj[1] += d[1];
    g_obj[2] += d[2];
    
    /* Add diagonal block */
    
    h_obj[0] += h[0];
    h_obj[1] += h[1];
    h_obj[2] += h[2];
    h_obj[3] += h[3];
    h_obj[4] += h[4];
    h_obj[5] += h[5];

    /* Add in objective function */
    o += f;
  }

  *obj = o;
  return 0;
}

