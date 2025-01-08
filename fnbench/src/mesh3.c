#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mesh.h"

int reorderMesh(Mesh *m) 
{
  const int n = m->nv;
  const int e = m->ne;

  double *v, *v2;
  int    *p, *p2;
  int    *t, *t2;

  int    *sta;
  int    *tet;

  int    *ord;
  int    *per;
  int    *pel;
  int    *iper;

  int    *que, *que2, *tmp;

  double  max, tst;
  int     idx, loc;
  int     ne, st, en;
  int     i, j;

  int     q, ql, ql2;

  if (n <= 0) {
    return 0;
  }

  /* Allocate memory for the list of elements containing a node and the      */
  /* permutation vector for the nodes.                                       */

  sta = (int *)malloc(sizeof(int)*(n+1));
  ord = (int *)malloc(sizeof(int)*n);
  per = (int *)malloc(sizeof(int)*n);
  tet = (int *)malloc(sizeof(int)*4*e);

  que  = (int *)malloc(sizeof(int)*n);
  que2 = (int *)malloc(sizeof(int)*n);

  memset(sta, 0, sizeof(int)*(n+1));
  memset(per, 0, sizeof(int)*n);

  /* Calculate the list of elements containing a node.                       */

  t = m->e;
  for (i = 0; i < e; ++i) {
    ++sta[t[0]];
    ++sta[t[1]];
    ++sta[t[2]];
    ++sta[t[3]];
    t += 4;
  }

  ne = sta[0];
  sta[0] = 0;
  for (i = 1; i <= n; ++i) {
    j = sta[i];
    sta[i] = ne;
    ne += j;
  }

  t = m->e;
  for (i = 0; i < e; ++i) {
    tet[sta[t[0]]++] = i;
    tet[sta[t[1]]++] = i;
    tet[sta[t[2]]++] = i;
    tet[sta[t[3]]++] = i;
    t += 4;
  }

  for (i = n; i > 0; --i) {
    sta[i] = sta[i-1];
  }
  sta[i] = 0;

  /* Find the node that is largest in magnitue that has not been marked.     */
  /* This will become the node from which we will start the Cuthill-McKee    */
  /* algorithm.                                                              */

  t = m->e;
  v = m->v;

  loc = 0;
  max =  0.0;
  idx = -1;

  for (i = 0; i < n; ++i) {
    tst = fabs(v[0]) + fabs(v[1]) + fabs(v[2]);
    if (tst >= max) {
      max = tst;
      idx = i;
    }
    v += 3;
  }

  while (idx >= 0) {

    /* We have a node that has not been visited.  Visit it now. */

    ql = 1;
    que[0] = idx;
    per[idx] = 1;

    while (ql) {
      q = 0;
      ql2 = 0;

      while (q < ql) {
        idx = que[q++];
        ord[loc++] = idx;

        st = sta[idx];
        en = sta[idx+1];
        while (st < en) {
          tmp = t + 4*tet[st++];

	  idx = tmp[0];
	  if (!per[idx]) {
	    que2[ql2++] = idx;
	    per[idx] = 1;
	  }

	  idx = tmp[1];
	  if (!per[idx]) {
	    que2[ql2++] = idx;
	    per[idx] = 1;
	  }

	  idx = tmp[2];
	  if (!per[idx]) {
	    que2[ql2++] = idx;
	    per[idx] = 1;
	  }

	  idx = tmp[3];
	  if (!per[idx]) {
	    que2[ql2++] = idx;
	    per[idx] = 1;
	  }
        }
      }

      ql   = ql2;

      tmp  = que;
      que  = que2;
      que2 = tmp;
    }

    if (loc >= n) {
      break;
    }

    v = m->v;
    max =  0.0;
    idx = -1;
    for (i = 0; i < n; ++i) {
      if (!per[i]) {
        tst = fabs(v[0]) + fabs(v[1]) + fabs(v[2]);
        if (tst >= max) {
          max = tst;
          idx = i;
        }
      }
      v += 3;
    }
  }

  free(que);
  free(que2);

  /* Calculate the permutation vector for the vertices and elements */

  pel = (int *)malloc(sizeof(int)*e);
  for (i = 0; i < e; ++i) {
    pel[i] = -1;
  }

  ne = 0;
  for (i = 0; i < n; ++i) {
    loc = ord[n-1-i];

    per[loc] = i;

    st = sta[loc];
    en = sta[loc+1];
    while (st < en) {
      loc = tet[st++];
      if (pel[loc] < 0) {
        pel[loc] = ne++;
      }
    }
  }

  free(ord);
  free(sta);
  free(tet);

  /* Reorder the nodes */

  v  = m->v;
  p  = m->p;

  v2   = (double *)malloc(sizeof(double)*3*n);
  p2   = (int    *)malloc(sizeof(int)*n);
  iper = (int    *)malloc(sizeof(int)*n);

  for (i = 0; i < n; ++i) {
    loc = per[i];
    iper[loc] = i;

    p2[loc  ] = *p++;
    
    loc = 3*per[i];
    v2[loc  ] = v[0];
    v2[loc+1] = v[1];
    v2[loc+2] = v[2];
    v += 3;
  }

  free(m->v);
  free(m->p);

  m->v = v2;
  m->p = p2;

  /* Reorder the elements */
  t  = m->e;
  t2 = (int *)malloc(sizeof(int)*4*e);

  for (i = 0; i < e; ++i) {
    loc = 4*pel[i];

    t2[loc  ] = per[t[0]];
    t2[loc+1] = per[t[1]];
    t2[loc+2] = per[t[2]];
    t2[loc+3] = per[t[3]];

    t += 4;
  }

  free(m->e);
  m->e = t2;

  /* Rerder the element categories */
  if (NULL != m->ecat) {
    t = m->ecat;
    t2 = (int *)malloc(sizeof(int)*e);

    for (i = 0; i < e; ++i) {
      t2[pel[i]] = t[0];
      ++t;
    }

    free(m->ecat);
    m->ecat = t2;
  }

  /* Rerder the element data */
  if (NULL != m->edat) {
    v = m->edat;
    v2 = (double *)malloc(sizeof(double)*e);

    for (i = 0; i < e; ++i) {
      v2[pel[i]] = v[0];
      ++v;
    }

    free(m->edat);
    m->edat = v2;
  }

  /* Track the permutations */
  m->per = per;
  m->pel = pel;
  m->iper = iper;
  return 0;
}

int finishMesh(Mesh *m)
{
  int    *e;
  int    *p;
  int    *ip;

  int     i, nv, ne;
  int     nf, nn, ndb, nob;

  nv = m->nv;
  ne = m->ne;

  e = m->e;
  p = m->p;

  nn = 0;
  nf = 0;
  for (i = 0; i < nv; ++i) {
    if (*p) {
      *p++ = -1;
      ++nf;
    }
    else {
      *p++ = nn;
      ++nn;
    }
  }

  m->nf = nf;
  m->nn = nn;

  if (nn <= 0) {
    fprintf(stderr, "All nodes are fixed.\n");

    m->i = (int *)malloc(sizeof(int));
    
    m->nb  = 0;
    m->ndb = 0;
    m->nob = 0;
    return -1;
  }

  m->i = (int *)malloc(sizeof(int)*nn);

  p  = m->p;
  ip = m->i;

  nn = 0;
  for (i = 0; i < nv; ++i) {
    if (*p++ >= 0) {
      ip[nn++] = i;
    }
  }

  p = m->p;
  nob = 0;
  ndb = 0;

  for (i = 0; i < ne; ++i) {
    nn = 0;
    if (p[e[0]] >= 0) { ++nn; }
    if (p[e[1]] >= 0) { ++nn; }
    if (p[e[2]] >= 0) { ++nn; }
    if (p[e[3]] >= 0) { ++nn; }
    
    ndb += nn;
    nob += (nn*(nn-1))/2;
    e += 4;
  }

  m->nb  = ndb + nob;
  m->ndb = ndb;
  m->nob = nob;
  return 0;
}
