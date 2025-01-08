#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mesh.h"

#define EQU_TOL 1e-10

int vert2(const void *p1, const void *p2)
{
  const double *x1 = (double *) p1;
  const double *x2 = (double *) p2;

  int i;
  
  for (i = 0; i < 2; ++i) {
    if (fabs(x1[i] - x2[i]) < EQU_TOL) {
      continue;
    }
    
    if (x1[i] < x2[i]) {
      return -1;
    }
    
    if (x1[i] > x2[i]) {
      return  1;
    }
  }
  return 0;
}

int vert3(const void *p1, const void *p2)
{
  const double *x1 = (double *) p1;
  const double *x2 = (double *) p2;

  int i;

  for (i = 0; i < 3; ++i) {
    if (fabs(x1[i] - x2[i]) < EQU_TOL) {
      continue;
    }

    if (x1[i] < x2[i]) {
      return -1;
    }

    if (x1[i] > x2[i]) {
      return  1;
    }
  }
  return 0;
}

int vert4(const void *p1, const void *p2)
{
  const double *x1 = (double *) p1;
  const double *x2 = (double *) p2;

  int i;

  for (i = 0; i < 4; ++i) {
    if (fabs(x1[i] - x2[i]) < EQU_TOL) {
      continue;
    }

    if (x1[i] < x2[i]) {
      return -1;
    }

    if (x1[i] > x2[i]) {
      return  1;
    }
  }
  return 0;
}

void face3(int *l)
{
  int min, idx;
  int i;
  
  min = l[0];
  idx = 0;
  
  for (i = 1; i < 3; ++i) {
    if (l[i] < min) {
      min = l[i];
      idx = i;
    }
  }
  
  switch(idx) {
  case 0:
    break;
    
  case 1:
    l[1] = l[2];
    l[2] = l[0];
    l[0] = min;
    break;
    
  case 2:
    l[2] = l[1];
    l[1] = l[0];
    l[0] = min;
    break;
  }
  return;
}

void face4(int *l)
{
  int min, idx;
  int i;

  min = l[0];
  idx = 0;

  for (i = 1; i < 4; ++i) {
    if (l[i] < min) {
      min = l[i];
      idx = i;
    }
  }

  switch(idx) {
  case 0:
    break;

  case 1:
    l[1] = l[2];
    l[2] = l[3];
    l[3] = l[0];
    l[0] = min;
    break;

  case 2:
    l[2] = l[3];
    l[3] = l[1];
    l[1] = l[2];
    l[2] = l[0];
    l[0] = min;
    break;

  case 3:
    l[3] = l[2];
    l[2] = l[1];
    l[1] = l[0];
    l[0] = min;
    break;
  }
  return;
}

void sort2(int *l)
{
  int t;

  if (l[0] > l[1]) {
    t = l[0];
    l[0] = l[1];
    l[1] = t;
  }
  return;
}

void sort3(int *l)
{
  int i, j, t;
  
  for (j = 2; j >= 1; --j) {
    for (i = 0; i < j; ++i) {
      if (l[i] > l[i+1]) {
        t = l[i];
        l[i] = l[i+1];
        l[i+1] = t;
      }
    }
  }
  return;
}

void sort4(int *l)
{
  int i, j, t;

  for (j = 3; j >= 1; --j) {
    for (i = 0; i < j; ++i) {
      if (l[i] > l[i+1]) {
        t = l[i];
        l[i] = l[i+1];
        l[i+1] = t;
      }
    }
  }
  return;
}

void sort8(int *l)
{
  int i, j, t;

  for (j = 7; j >= 1; --j) {
    for (i = 0; i < j; ++i) {
      if (l[i] > l[i+1]) {
        t = l[i];
        l[i] = l[i+1];
        l[i+1] = t;
      }
    }
  }
  return;
}

void radix2(int *dest, int *src, int *len,
	    const int idx, const int nv, const int ne)
{
  int i, j, n;

  memset(len, 0, sizeof(int)*nv);
  for (i = 0; i < ne; ++i) {
    ++len[src[idx]];
    src += 4;
  }
  src -= 4*ne;

  n = len[0];
  len[0] = 0;
  for (i = 1; i < nv; ++i) {
    j = len[i];
    len[i] = n;
    n += j;
  }

  for (i = 0; i < ne; ++i) {
    n = 4*len[src[idx]]++;

    dest[n  ] = src[0];
    dest[n+1] = src[1];
    dest[n+2] = src[2];
    dest[n+3] = src[3];
    src += 4;
  }
  return;
}

void radix3(int *dest, int *src, int *len,
	    const int idx, const int nv, const int ne)
{
  int i, j, n;

  memset(len, 0, sizeof(int)*nv);
  for (i = 0; i < ne; ++i) {
    ++len[src[idx]];
    src += 5;
  }
  src -= 5*ne;

  n = len[0];
  len[0] = 0;
  for (i = 1; i < nv; ++i) {
    j = len[i];
    len[i] = n;
    n += j;
  }

  for (i = 0; i < ne; ++i) {
    n = 5*len[src[idx]]++;

    dest[n  ] = src[0];
    dest[n+1] = src[1];
    dest[n+2] = src[2];
    dest[n+3] = src[3];
    dest[n+4] = src[4];
    src += 5;
  }
  return;
}

void radix4(int *dest, int *src, int *len,
	    const int idx, const int nv, const int ne)
{
  int i, j, n;

  memset(len, 0, sizeof(int)*nv);
  for (i = 0; i < ne; ++i) {
    ++len[src[idx]];
    src += 6;
  }
  src -= 6*ne;

  n = len[0];
  len[0] = 0;
  for (i = 1; i < nv; ++i) {
    j = len[i];
    len[i] = n;
    n += j;
  }

  for (i = 0; i < ne; ++i) {
    n = 6*len[src[idx]]++;

    dest[n  ] = src[0];
    dest[n+1] = src[1];
    dest[n+2] = src[2];
    dest[n+3] = src[3];
    dest[n+4] = src[4];
    dest[n+5] = src[5];
    src += 6;
  }
  return;
}

int allocMeshData(MeshData **m)
{
  if (NULL != *m) {
    return -1;
  }

  (*m) = (MeshData *)malloc(sizeof(MeshData));
  (*m)->v = NULL;
  (*m)->e = NULL;
  (*m)->f = NULL;
  (*m)->b = NULL;
  (*m)->part = NULL;
  (*m)->se = NULL;
  (*m)->si = NULL;
  (*m)->ss = NULL;
  return 0; 
}

int allocMesh(Mesh **m)
{
  if (NULL != *m) {
    return -1;
  }

  (*m) = (Mesh *)malloc(sizeof(Mesh));

  (*m)->d = NULL;

  (*m)->v = NULL;
  (*m)->e = NULL;

  (*m)->p = NULL;
  (*m)->i = NULL;
  (*m)->g = NULL;
  
  (*m)->len = NULL;
  (*m)->col = NULL;
  (*m)->dat = NULL;
  (*m)->inst = NULL;
  
  (*m)->per = NULL;
  return 0;
}

int freeMeshData(MeshData **m)
{
  if (NULL == *m) {
    return -1;
  }

  if (NULL != (*m)->v) { free((*m)->v); }
  if (NULL != (*m)->e) { free((*m)->e); }
  if (NULL != (*m)->f) { free((*m)->f); }
  if (NULL != (*m)->b) { free((*m)->b); }
  if (NULL != (*m)->part) { free((*m)->part); }
  if (NULL != (*m)->se) { free((*m)->se); }
  if (NULL != (*m)->si) { free((*m)->si); }
  if (NULL != (*m)->ss) { free((*m)->ss); }

  free(*m);
  (*m) = NULL;
  return 0;
}

int freeMesh(Mesh **m)
{
  if (NULL == *m) {
    return -1;
  }

  if (NULL != (*m)->d) { freeMeshData(&((*m)->d)); }

  if (NULL != (*m)->v) { free((*m)->v); }
  if (NULL != (*m)->e) { free((*m)->e); }

  if (NULL != (*m)->p) { free((*m)->p); }
  if (NULL != (*m)->i) { free((*m)->i); }
  if (NULL != (*m)->g) { free((*m)->g); }
  if (NULL != (*m)->len) { free((*m)->len); }
  if (NULL != (*m)->col) { free((*m)->col); }
  if (NULL != (*m)->dat) { free((*m)->dat); }
  if (NULL != (*m)->inst) { free((*m)->inst); }
  if (NULL != (*m)->per) { free((*m)->per); }

  free(*m);
  (*m) = NULL;
  return 0;
}

