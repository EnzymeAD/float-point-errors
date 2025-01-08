#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mesh.h"
#include "fcn.h"

static int readMeshData(const char *fname, MeshData *m)
{
  FILE   *fp;

  double *v;
  int    *e;
  int    *f;
  int    *part;

  int i, na, nv, ne, nf, np, no;

  int errs = 0;
  int line = 0;

  char buf[1024];

  /***************************************************************************/
  /* Open the file to read the mesh data from.                               */
  /***************************************************************************/

  fp = fopen(fname, "r");
  if (NULL == fp) {
    printf("Invalid input mesh: %s\n", fname);
    return -1;
  }

  /***************************************************************************/
  /* First line is an integer specifying the number of vertices.             */
  /* Optionally: number of partitions (default 1) and amount of overlap      */
  /* (default 1).                                                            */
  /***************************************************************************/

  buf[0] = '\0';
  fgets(buf, 1024, fp); ++line;

  na = sscanf(buf, "%d %d %d", &nv, &np, &no);
  if ((na < 1) || (nv <= 0)) {
    printf("Read error on line %d: %s\n", line, 
           "expecting positive number of vertices\n");
    return -2;
  }

  if (na >= 2) {
    if (np <= 0) {
      printf("Read error on line %d: %s\n", line, 
             "expecting positive number of partitions\n");
      return -2;
    }
  }
  else {
    np = 1;
  }

  if ((np >= 2) && (na >= 3)) {
    if (no <= 0) {
      printf("Read error on line %d: %s\n", line, 
             "expecting positive amount of overlap\n");
      return -2;
    }
  }
  else {
    no = 1;
  }

  m->nv = nv;
  m->np = np;
  m->no = no;

  m->v = (double *)malloc(3*sizeof(double)*nv);
  v = m->v;

  if (1 == np) {
    /*************************************************************************/
    /* Read three coordinates for each vertex.                               */
    /*************************************************************************/

    for (i = 0; i < nv; ++i) {
      buf[0] = '\0';
      fgets(buf, 1024, fp); ++line;
      if (3 != sscanf(buf, "%lf%lf%lf", v, v+1, v+2)) {
        printf("Read error on line %d: %s\n", line,
	       "expecting three coordinate values");
        errs = 1;
      }

      v += 3;
    }
  }
  else {
    /*************************************************************************/
    /* Read three coordinates for each vertex + partition number             */
    /*************************************************************************/

    m->part = (int *)malloc(sizeof(int)*nv);
    part = m->part;

    for (i = 0; i < nv; ++i) {
      buf[0] = '\0';
      fgets(buf, 1024, fp); ++line;
      if (4 != sscanf(buf, "%lf%lf%lf%d", v, v+1, v+2, part)) {
        printf("Read error on line %d: %s\n", line,
	       "expecting three coordinate values and parition");
        errs = 1;
      }

      if ((*part < 0) || (*part >= np)) {
        printf("Error on line %d: %s\n", line,
	       "parition not in valid range");
        errs = 1;
      }

      v += 3;
      ++part;
    }
  }

  /***************************************************************************/
  /* Next is an integer specifying the number of elements.                   */
  /***************************************************************************/

  buf[0] = '\0';
  fgets(buf, 1024, fp); ++line;
  if ((1 != sscanf(buf, "%d", &ne)) || (ne <= 0)) {
    printf("Read error on line %d: %s\n", line, 
           "expecting positive number of elements");
    fclose(fp);
    return -2;
  }

  m->ne = ne;
  m->e = (int *)malloc(4*sizeof(int)*ne);

  /***************************************************************************/
  /* Read four indices for each vertex.                                      */
  /***************************************************************************/

  e = m->e;
  for (i = 0; i < ne; ++i) {
    buf[0] = '\0';
    fgets(buf, 1024, fp); ++line;
    if (4 != sscanf(buf, "%d%d%d%d", e, e+1, e+2, e+3)) {
      printf("Read error on line %d: %s\n", line,
	     "expecting four integer vertex indices");
      errs = 1;
    }

    e += 4;
  }

  /***************************************************************************/
  /* Next is an integer specifying the number of fixed vertices (optional).  */
  /***************************************************************************/

  buf[0] = '\0';
  fgets(buf, 1024, fp); ++line;
  if ((1 != sscanf(buf, "%d", &nf)) || (nf < 0)) {
    nf = 0;
  }

  m->nf = nf;
  if (nf > 0) {
    m->f = (int *)malloc(sizeof(int)*nf);

    /*************************************************************************/
    /* Read one index for each vertex.                                       */
    /*************************************************************************/

    f = m->f;
    for (i = 0; i < nf; ++i) {
      buf[0] = '\0';
      fgets(buf, 1024, fp); ++line;
      if (1 != sscanf(buf, "%d", f)) {
        printf("Read error on line %d: %s\n", line,
	       "expecting one integer vertex index");
        errs = 1;
      }

      f += 1;
    }
  }

  fclose(fp);
  return errs;
}

static int checkMeshData(MeshData *m)
{
  int *e = m->e;
  int *f = m->f;

  int nv = m->nv;
  int ne = m->ne;
  int nf = m->nf;

  int i, j;

  int errs = 0;
  
  for (i = 0; i < ne; ++i) {
    for (j = 0; j < 4; ++j) {
      if ((e[j] < 0) || (e[j] >= nv)) {
        printf("Index error for element %d\n", i);
	errs = 1;
	break;
      }
    }
    e += 4;
  }
 
  for (i = 0; i < nf; ++i) {
    if ((f[0] < 0) || (f[0] >= nv)) {
      printf("Index error for fixed vertex %d\n", i);
      errs = 1;
    }
    f += 1;
  }

  return errs;
}

static void sort3(int *l)
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

static void radix3(int *dest, int *src, int *len,
                   const int idx, const int nv, const int ne)
{
  int i, j, n;

  memset(len, 0, sizeof(int)*nv);
  for (i = 0; i < ne; ++i) {
    ++len[src[idx]];
    src += 3;
  }
  src -= 3*ne;

  n = len[0];
  len[0] = 0;
  for (i = 1; i < nv; ++i) {
    j = len[i];
    len[i] = n;
    n += j;
  }

  for (i = 0; i < ne; ++i) {
    n = 3*len[src[idx]]++;

    dest[n  ] = src[0];
    dest[n+1] = src[1];
    dest[n+2] = src[2];
    src += 3;
  }
  return;
}

static int boundaryMeshData(MeshData *m)
{
  int *e = m->e;
  int *f = m->f;
  int *b;

  int *l1, *l2, *ls;

  int nv = m->nv;
  int ne = m->ne;
  int nf = m->nf;

  int i, j;

  m->b = (int *)malloc(sizeof(int)*nv);
  b = m->b;

  memset(b, 0, sizeof(int)*nv);

  for (i = 0; i < ne; ++i) {
    for (j = 0; j < 4; ++j) {
      ++b[e[j]];
    }
    e += 4;
  }
 
  for (i = 0; i < nf; ++i) {
    ++b[f[0]];
    f += 1;
  }

  for (i = 0; i < nv; ++i) {
    if (0 == b[0]) {
      printf("Unreferenced vertex %d\n", i);
      b[0] = 1;
    }
    else {
      b[0] = 0;
    }
    ++b;
  }

  f = m->f;
  b = m->b;
  for (i = 0; i < nf; ++i) {
    b[f[0]] = 1;
    f += 1;
  }

  ls = (int *)malloc(   sizeof(int)*nv);         /* vertex starts */
  l1 = (int *)malloc(12*sizeof(int)*(ne+1));     /* list of faces */
  l2 = (int *)malloc(12*sizeof(int)*(ne+1));     /* list of faces */

  e = m->e;
  for (i = 0; i < ne; ++i) {
    l1[0] = e[0]; l1[1] = e[1]; l1[2] = e[2];
    sort3(l1); l1 += 3;

    l1[0] = e[0]; l1[1] = e[1]; l1[2] = e[3];
    sort3(l1); l1 += 3;

    l1[0] = e[0]; l1[1] = e[2]; l1[2] = e[3];
    sort3(l1); l1 += 3;

    l1[0] = e[1]; l1[1] = e[2]; l1[2] = e[3];
    sort3(l1); l1 += 3;

    e += 4;
  }

  l1[0] = -1;
  l1[1] = -1;
  l1[2] = -1;

  l1 -= 12*ne;

  /* Now perform the radix sort */
  radix3(l2, l1, ls, 2, nv, 4*ne);
  radix3(l1, l2, ls, 1, nv, 4*ne);
  radix3(l2, l1, ls, 0, nv, 4*ne);
  memcpy(l1, l2, 12*sizeof(int)*ne);

  free(ls);
  free(l2);

  for (i = 0; i < 4*ne; ++i) {
    if ((l1[0] == l1[3]) && (l1[1] == l1[4]) && (l1[2] == l1[5])) {

      while ((i < 4*ne) &&
             (l1[0] == l1[3]) && (l1[1] == l1[4]) && (l1[2] == l1[5])) {
        l1 += 3;
        ++i;
      }
    }
    else {
      b[l1[0]] = 1;
      b[l1[1]] = 1;
      b[l1[2]] = 1;
    }

    l1 += 3;
  }

  l1 -= 12*ne;
  free(l1);
  return 0;
}

static int computeMesh(Mesh *m) 
{
  const MeshData *md = m->d;

  const int hn = md->nv;
  const int he = md->ne;

  m->nv = hn;
  m->ne = he;

  m->v = (double *)malloc(3*sizeof(double)*hn);
  m->e = (int    *)malloc(4*sizeof(int   )*he);
  m->p = (int    *)malloc(1*sizeof(int   )*hn);

  memcpy(m->v, md->v, 3*sizeof(double)*hn);
  memcpy(m->e, md->e, 4*sizeof(int   )*he);
  memcpy(m->p, md->b,   sizeof(int   )*hn);
  return 0;
}

static int checkMesh(Mesh *m) 
{
  const int ne = m->d->ne;

  double *v = m->v;
  double *w;
  int    *e = m->e;

  double  x[12];
  double  f;
  int     v1, v2, v3, v4;
  int     i, errs = 0;
  
  for (i = 0; i < ne; ++i) {
    v1 = e[0];
    v2 = e[1];
    v3 = e[2];
    v4 = e[3];
    e += 4;

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

    if (o_fcn(&f, x)) {
      if (f < -epsilon) {
	printf("Invalid element: %d: inverted.\n", i);
      } 
      else {
	printf("Invalid element: %d: zero area: %5.4e.\n", i, f);
      }
      ++errs;
    }
  }
  return errs;
}

int readMesh(const char *fname, Mesh **m)
{
  MeshData *md = NULL;

  if (allocMeshData(&md)) {
    printf("Allocation failure.\n");
    return -1;
  }

  if (readMeshData(fname, md)) {
    printf("Read failure.\n");
    freeMeshData(&md);
    return -2;
  }

  if (checkMeshData(md)) {
    printf("Index check failure.\n");
    freeMeshData(&md);
    return -2;
  }

  if (boundaryMeshData(md)) {
    printf("Boundary calculation failure.\n");
    freeMeshData(&md);
    return -2;
  }

  if (allocMesh(m)) {
    printf("Mesh allocation failure.\n");
    freeMeshData(&md);
    return -1;
  }

  (*m)->d = md;

  if (computeMesh(*m)) {
    printf("Mesh construction failure.\n");
    freeMesh(m);
    return -2;
  }

  if (checkMesh(*m)) {
    printf("Inverted element failure.\n");
    freeMesh(m);
    return -2;
  }

#ifdef REORDER
  if (reorderMesh(*m)) {
    printf("Reordering failure.\n");
    freeMesh(m);
    return -2;
  }
#endif

  if (finishMesh(*m)) {
    printf("Mesh initialization failure.\n");
    freeMesh(m);
    return -2;
  }
  return 0;
}

static int writeMeshData(const char *fname, MeshData *m)
{
  FILE   *fp;

  double *v;
  int    *e;
  int    *f;
  int     i, nv, ne, nf;

  fp = fopen(fname, "w");
  if (NULL == fp) {
    fprintf(stderr, "Could not open output mesh.\n");
    return -1;
  }

  nv = m->nv;
  ne = m->ne;
  nf = m->nf;

  fprintf(fp, "%d\n", nv);

  v = m->v;
  for (i = 0; i < nv; ++i) {
    fprintf(fp, "%15.14e %15.14e %15.14e\n", v[0], v[1], v[2]);
    v += 3;
  }

  fprintf(fp, "%d\n", ne);

  e = m->e;
  for (i = 0; i < ne; ++i) {
    fprintf(fp, "%d %d %d %d\n", e[0], e[1], e[2], e[3]);
    e += 4;
  }

  if (nf) {
    fprintf(fp, "%d\n", nf);

    f = m->f;
    for (i = 0; i < nf; ++i) {
      fprintf(fp, "%d\n", f[0]);
      f += 1;
    }
  }

  fclose(fp);
  return 0;
}

int writeMesh(const char *fname, Mesh *m)
{
  const int n = m->nv;

  double *hv = m->d->v;
  double *tv = m->v;

#ifdef REORDER
  int    *tp = m->per;
  int i, loc;

  for (i = 0; i < n; ++i) {
    loc = 3*(*tp++);

    hv[0] = tv[loc  ];
    hv[1] = tv[loc+1];
    hv[2] = tv[loc+2];

    hv += 3;
  }
#else
  memcpy(hv, tv, 3*sizeof(double)*n);
#endif

  return writeMeshData(fname, m->d);
}

static int writeMeshData_AMPL(const char *fname, MeshData *m)
{
  FILE   *fp;

  double *v;
  int    *e;
  int     i, nv, ne;

  fp = fopen(fname, "w");
  if (NULL == fp) {
    fprintf(stderr, "Could not open output mesh.\n");
    return -1;
  }

  nv = m->nv;
  ne = m->ne;

  fprintf(fp, "param V := %d;\n", nv);
  fprintf(fp, "param E := %d;\n", ne);

  fprintf(fp, "\nvar x : 1 2 3 =\n");
  v = m->v;
  for (i = 0; i < nv; ++i) {
    fprintf(fp, "%6d %15.14e %15.14e %15.14e\n", i+1, v[0], v[1], v[2]);
    v += 3;
  }

  fprintf(fp, ";\n\nparam TETS : 1 2 3 4 =\n");
  e = m->e;
  for (i = 0; i < ne; ++i) {
    fprintf(fp, "%6d %6d %6d %6d %6d\n", i+1, e[0]+1, e[1]+1, e[2]+1, e[3]+1);
    e += 4;
  }
  fprintf(fp, ";\n\n");

  for (i = 0; i < nv; ++i) {
    if (m->b[i]) {
      fprintf(fp, "fix {i in COORDS} x[%d,i];\n", i+1);
    }
  }

  fclose(fp);
  return 0;
}

static int writeMeshOpt_AMPL(const char *fname, Mesh *m)
{
  FILE   *fp;

  double *v;
  int    *e;
  int     i, nv, ne;

  fp = fopen(fname, "w");
  if (NULL == fp) {
    fprintf(stderr, "Could not open output mesh.\n");
    return -1;
  }

  nv = m->nv;
  ne = m->ne;

  fprintf(fp, "param V := %d;\n", nv);
  fprintf(fp, "param E := %d;\n", ne);

  fprintf(fp, "\nvar x : 1 2 3 =\n");
  v = m->v;
  for (i = 0; i < nv; ++i) {
    fprintf(fp, "%6d %15.14e %15.14e %15.14e\n", i+1, v[0], v[1], v[2]);
    v += 3;
  }

  fprintf(fp, ";\n\nparam TETS : 1 2 3 4 =\n");
  e = m->e;
  for (i = 0; i < ne; ++i) {
    fprintf(fp, "%6d %6d %6d %6d %6d\n", i+1, e[0]+1, e[1]+1, e[2]+1, e[3]+1);
    e += 4;
  }
  fprintf(fp, ";\n\n");

  for (i = 0; i < nv; ++i) {
    if (m->p[i] < 0) {
      fprintf(fp, "fix {i in COORDS} x[%d,i];\n", i+1);
    }
  }

  fclose(fp);
  return 0;
}

int writeMesh_AMPL(const char *fname, const char *fname_opt, Mesh *m)
{
  const int n = m->nv;

  double *hv = m->d->v;
  double *tv = m->v;

#ifdef REORDER
  int    *tp = m->per;
  int i, loc;

  for (i = 0; i < n; ++i) {
    loc = 3*(*tp++);

    hv[0] = tv[loc  ];
    hv[1] = tv[loc+1];
    hv[2] = tv[loc+2];

    hv += 3;
  }
#else
  memcpy(hv, tv, 3*sizeof(double)*n);
#endif

  writeMeshOpt_AMPL(fname_opt, m);
  return writeMeshData_AMPL(fname, m->d);
}

