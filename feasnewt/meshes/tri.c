#include <assert.h>
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

  m->v = (double *)malloc(2*sizeof(double)*nv);
  v = m->v;

  if (1 == np) {
    /*************************************************************************/
    /* Read two coordinates for each vertex.                                 */
    /*************************************************************************/

    for (i = 0; i < nv; ++i) {
      buf[0] = '\0';
      fgets(buf, 1024, fp); ++line;
      if (2 != sscanf(buf, "%lf%lf", v, v+1)) {
        printf("Read error on line %d: %s\n", line,
	       "expecting two coordinate values");
        errs = 1;
      }

      v += 2;
    }
  }
  else {
    /*************************************************************************/
    /* Read two coordinates for each vertex + partition number               */
    /*************************************************************************/

    m->part = (int *)malloc(sizeof(int)*nv);
    part = m->part;

    for (i = 0; i < nv; ++i) {
      buf[0] = '\0';
      fgets(buf, 1024, fp); ++line;
      if (3 != sscanf(buf, "%lf%lf%d", v, v+1, part)) {
        printf("Read error on line %d: %s\n", line,
	       "expecting two coordinate values and parition");
        errs = 1;
      }

      if ((*part < 0) || (*part >= np)) {
        printf("Error on line %d: %s\n", line,
	       "parition not in valid range");
        errs = 1;
      }

      v += 2;
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
  m->e = (int *)malloc(3*sizeof(int)*ne);

  /***************************************************************************/
  /* Read three indices for each vertex.                                     */
  /***************************************************************************/

  e = m->e;
  for (i = 0; i < ne; ++i) {
    buf[0] = '\0';
    fgets(buf, 1024, fp); ++line;
    if (3 != sscanf(buf, "%d%d%d", e, e+1, e+2)) {
      printf("Read error on line %d: %s\n", line,
	     "expecting three integer vertex indices");
      errs = 1;
    }

    e += 3;
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

#ifdef REFINE
#ifdef SPLINTER
  // This refinement creates splinters in the mesh and can only be
  // user a small number of times before the refined mesh becomes
  // impossible to optimize due to numerical errors
  if (!errs) {
    double *newv;
    int *newe;

    newv = (double *)malloc(2*sizeof(double)*(nv + ne));
    newe = (int *)malloc(9*sizeof(int)*ne);

    memcpy(newv, m->v, 2*sizeof(double)*nv);
    for (i = 0; i < ne; ++i) {
      newv[2*(nv+i)+0] = (m->v[2*m->e[3*i+0]+0] +
                          m->v[2*m->e[3*i+1]+0] +
                          m->v[2*m->e[3*i+2]+0]) / 3.0;

      newv[2*(nv+i)+1] = (m->v[2*m->e[3*i+0]+1] +
                          m->v[2*m->e[3*i+1]+1] +
                          m->v[2*m->e[3*i+2]+1]) / 3.0;
    }

    free(m->v);
    m->v = newv;
    m->nv = nv+ne;

    for (i = 0; i < ne; ++i) {
      newe[9*i+0] = m->e[3*i+0];
      newe[9*i+1] = m->e[3*i+1];
      newe[9*i+2] = nv + i;

      newe[9*i+3] = m->e[3*i+1];
      newe[9*i+4] = m->e[3*i+2];
      newe[9*i+5] = nv + i;

      newe[9*i+6] = m->e[3*i+2];
      newe[9*i+7] = m->e[3*i+0];
      newe[9*i+8] = nv + i;
    }

    free(m->e);
    m->e = newe;
    m->ne = 3*ne;
  }

#else
  // This version is a "uniform" refinement of the mesh.  Three vertices
  // are added per element, one vertex that bisects each edge.  The
  // triangle is then split into four pieces.  If we start with an
  // equilateral triangle, then the refined one will be four equilateral
  // triangles.

  if (!errs) {
    double *newv;
    int *perm;
    int *newe;

    int nnv, vert, j;

    // We start off by adding the vertices to the mesh.  One vertex is
    // added per edge and most edges appear in many elements.  The upper
    // bound on the number of new vertices 3*ne, but the final count will
    // be much lower once the duplicates are removed.  For the vertices,
    // we store three values, two doubles for the coordinates and an
    // integer for the initial position.  We reallocate the data once
    // we know the correct sizes.

    newv = (double *)malloc(3*sizeof(double)*(nv + 3*ne));
    perm = (int *)malloc(sizeof(int)*(nv + 3*ne));

    // There will be four elements in the new discretization per element in
    // the old discretization; each element is a set of three indices.
    // We do not need to reallocate this data.

    newe = (int *)malloc(12*sizeof(int)*ne);

    // Copy the original vertex list into the allocated vertex list and
    // add the vertex numbers.

    for (i = 0; i < nv; ++i) {
      newv[3*i + 0] = m->v[2*i + 0];
      newv[3*i + 1] = m->v[2*i + 1];
      newv[3*i + 2] = i;
      perm[i] = -1;
    }

    // Add the new vertices and elements

    for (i = 0; i < ne; ++i) {
      // First edge: a: (0,1)
      newv[3*(nv + 3*i + 0) + 0] = 0.5*(m->v[2*m->e[3*i + 0] + 0] +
                                        m->v[2*m->e[3*i + 1] + 0]);
      newv[3*(nv + 3*i + 0) + 1] = 0.5*(m->v[2*m->e[3*i + 0] + 1] +
                                        m->v[2*m->e[3*i + 1] + 1]);
      newv[3*(nv + 3*i + 0) + 2] = nv + 3*i + 0;
      perm[nv + 3*i + 0] = -1;

      // Second edge: b: (0,2)
      newv[3*(nv + 3*i + 1) + 0] = 0.5*(m->v[2*m->e[3*i + 0] + 0] +
                                        m->v[2*m->e[3*i + 2] + 0]);
      newv[3*(nv + 3*i + 1) + 1] = 0.5*(m->v[2*m->e[3*i + 0] + 1] +
                                        m->v[2*m->e[3*i + 2] + 1]);
      newv[3*(nv + 3*i + 1) + 2] = nv + 3*i + 1;
      perm[nv + 3*i + 1] = -1;

      // Third edge: c: (1,2)
      newv[3*(nv + 3*i + 2) + 0] = 0.5*(m->v[2*m->e[3*i + 1] + 0] +
                                        m->v[2*m->e[3*i + 2] + 0]);
      newv[3*(nv + 3*i + 2) + 1] = 0.5*(m->v[2*m->e[3*i + 1] + 1] +
                                        m->v[2*m->e[3*i + 2] + 1]);
      newv[3*(nv + 3*i + 2) + 2] = nv + 3*i + 2;
      perm[nv + 3*i + 2] = -1;

      // First element (0, a, b)
      newe[12*i +  0] = m->e[3*i + 0];
      newe[12*i +  1] = nv + 3*i + 0;
      newe[12*i +  2] = nv + 3*i + 1;

      // Second element (1, c, a)
      newe[12*i +  3] = m->e[3*i + 1];
      newe[12*i +  4] = nv + 3*i + 2;
      newe[12*i +  5] = nv + 3*i + 0;

      // Third element (2, b, c)
      newe[12*i +  6] = m->e[3*i + 2];
      newe[12*i +  7] = nv + 3*i + 1;
      newe[12*i +  8] = nv + 3*i + 2;

      // Fourth element (a, c, b)
      newe[12*i +  9] = nv + 3*i + 0;
      newe[12*i + 10] = nv + 3*i + 2;
      newe[12*i + 11] = nv + 3*i + 1;
    }

    // Now sort the vertices by their coordinates and index number
    qsort(newv, nv + 3*ne, 3*sizeof(double), vert3);

    // Compute the permutation
    nnv = nv;
    for (i = 0; i < nv + 3*ne; ) {
      // Get the current vertex number
      vert = (int) newv[3*i + 2];

      if (vert < nv) {
        // Current vertex is an original vertex assumed unique
        perm[vert] = vert;
        assert(vert2(newv + 3*i, newv + 3*(i+1)));
        ++i;
      }
      else {
        // Current vertex is a new vertex; permute to next location
        // Get the vertices that are the same and permute them
        for (j = i; j < nv + 3*ne; ++j) {
          if (vert2(newv + 3*i, newv + 3*j)) {
            // i and j not the same
            break;
          }
          perm[(int) newv[3*j + 2]] = nnv;
        }
        ++nnv;

        // Skip to the next distinct vertex
        i = j;
      }
    }

    // No permutation for initial vertices
    for (i = 0; i < nv; ++i) {
      assert(perm[i] == i);
    }

    // Some permutations for rest
    for (; i < nv + 3*ne; ++i) {
      assert(perm[i] > 0 && perm[i] < nnv);
    }

    // Now permute the vertices
    free(m->v);
    m->v = (double *)malloc(2*sizeof(double)*nnv);
    m->nv = nnv;

    for (i = 0; i < nv + 3*ne; ++i) {
      vert = perm[(int) newv[3*i + 2]];

      m->v[2*vert + 0] = newv[3*i + 0];
      m->v[2*vert + 1] = newv[3*i + 1];
    }

    free(newv);

    // Now permute the elements

    free(m->e);
    m->e = newe;
    m->ne = 4*ne;

    for (i = 0; i < 4*ne; ++i) {
      m->e[3*i + 0] = perm[m->e[3*i + 0]];
      m->e[3*i + 1] = perm[m->e[3*i + 1]];
      m->e[3*i + 2] = perm[m->e[3*i + 2]];
    }

    free(perm);
  }

#endif
#endif

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
    for (j = 0; j < 3; ++j) {
      if ((e[j] < 0) || (e[j] >= nv)) {
        printf("Index error for element %d\n", i);
	errs = 1;
	break;
      }
    }
    e += 3;
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

#ifdef CHECK_MESH
static int checkVertData(MeshData *m)
{
  double *v = m->v;
  int nv = m->nv;

  double *vertData, *vd;
  int i;

  /* Each record is:                 */
  /* | c1 | c2 | vert # | */

  vertData = (double *)malloc(3*sizeof(double)*nv);

  vd = vertData;
  for (i = 0; i < nv; ++i) {
    vd[0] = v[0];
    vd[1] = v[1];
    vd[2] = i;

    v += 2;
    vd += 3;
  }

  vd = vertData;
  qsort(vd, nv, 3*sizeof(double), vert2);

  for (i = 0; i < nv - 1; ++i) {
    if (!vert2(vd, vd+3)) {
      /* Two vertices match */
      printf("Same vertex: %d %d: %5.4e\n", (int) vd[2], (int) vd[5],
	     fabs(vd[0] - vd[3]) + fabs(vd[1] - vd[4]));

#ifdef HEAL_MESH
      {
        const int o = (int) vd[2];
        const int n = (int) vd[5];
	const int ne = m->ne;
        int *ed = m->e;
	int j, k;

        for (j = 0; j < ne; ++j) {
          for (k = 0; k < 3; ++k) {
            if (ed[k] == o) {
              ed[k] = n;
              printf("Fixing element %d vertex %d\n", j, k);
            }
          }
	  ed += 3;
        }
      }
#endif
    }
    vd += 3;
  }
  
  free(vertData);
  return 0;
}
#endif

static int boundaryMeshData(MeshData *m)
{
  int *e = m->e;
  int *f = m->f;
  int *b;

  int *l1, *l2, *ls;

  int nv = m->nv;
  int ne = m->ne;
  int nf = m->nf;

  int i, j, cnt;

  m->b = (int *)malloc(sizeof(int)*nv);
  b = m->b;

  memset(b, 0, sizeof(int)*nv);

  for (i = 0; i < ne; ++i) {
    for (j = 0; j < 3; ++j) {
      ++b[e[j]];
    }
    e += 3;
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

#ifdef HEAL_MESH
  {
    /* Remove unreferenced vertices */
    double *v = m->v;
    int tv, k;

    tv = 0;
    b = m->b;
    for (i = 0; i < nv; ++i) {
      if (0 == b[i]) {
        /* vertex is referenced */
        b[tv] = 0;
        v[2*tv + 0] = v[2*i + 0];
        v[2*tv + 1] = v[2*i + 1];
  
        e = m->e;
        for (j = 0; j < ne; ++j) {
          for (k = 0; k < 3; ++k) {
            if (e[k] == i) {
              e[k] = tv;
            }
          }
          e += 3;
        }
        ++tv;
      }
    }
    m->nv = tv;
    nv = tv;
  }
#endif

  f = m->f;
  b = m->b;
  for (i = 0; i < nf; ++i) {
    b[f[0]] = 1;
    f += 1;
  }

  ls = (int *)malloc(  sizeof(int)*nv);         /* vertex starts */
  l1 = (int *)malloc(12*sizeof(int)*(ne+1));     /* list of faces */
  l2 = (int *)malloc(12*sizeof(int)*(ne+1));     /* list of faces */

  e = m->e;
  for (i = 0; i < ne; ++i) {
    l1[0] = e[0]; l1[1] = e[1];
    l1[2] = i; l1[3] = 0;
    sort2(l1); l1 += 4;

    l1[0] = e[0]; l1[1] = e[2];
    l1[2] = i; l1[3] = 1;
    sort2(l1); l1 += 4;

    l1[0] = e[1]; l1[1] = e[2];
    l1[2] = i; l1[3] = 2;
    sort2(l1); l1 += 4;

    e += 3;
  }

  l1[0] = -1;
  l1[1] = -1;
  l1[2] = -1;
  l1[3] = -1;

  l1 -= 12*ne;

  /* Now perform the radix sort */
  radix2(l2, l1, ls, 1, nv, 3*ne);
  radix2(l1, l2, ls, 0, nv, 3*ne);

  free(ls);
  free(l2);

  for (i = 0; i < 3*ne; ++i) {
    if ((l1[0] == l1[4]) && (l1[1] == l1[5])) {
      /* Face must be in only two elements! */
      cnt = 1;

#ifdef CHECK_MESH
      /* Check compatibility conditions for the faces             */
      /* Make sure one is clockwise and other is counter clockwise*/
      e = m->e + 3*l1[2];

      switch(l1[3]) {
      case 0:
        l1[0] = e[0]; l1[1] = e[1];
        break;

      case 1:
        l1[0] = e[2]; l1[1] = e[0];
        break;

      case 2:
        l1[0] = e[1]; l1[1] = e[2];
        break;
      }

      e = m->e + 3*l1[6];
      switch(l1[7]) {
      case 0:
        l1[4] = e[1]; l1[5] = e[0];
        break;

      case 1:
        l1[4] = e[0]; l1[5] = e[2];
        break;

      case 2:
        l1[4] = e[2]; l1[5] = e[1];
        break;
      }

      if ((l1[0] != l1[4]) || (l1[1] != l1[5])) {
        printf("Face Compatiblity: %d.%d and %d.%d\n",
               l1[2], l1[3], l1[6], l1[7]);
      }

      sort2(l1);
      sort2(l1 + 4);
#endif

      /* Remove all replicated vertices */
      while ((i < 3*ne) && (l1[0] == l1[4]) && (l1[1] == l1[5])) {
        l1 += 4;
        ++i;
	++cnt;
      }

      if (cnt > 2) {
	printf("Face Count Problem: %d.%d = %d\n", l1[2], l1[3], cnt);
      }
    }
    else {
      b[l1[0]] = 1;
      b[l1[1]] = 1;
    }

    l1 += 4;
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

  m->v = (double *)malloc(2*sizeof(double)*hn);
  m->e = (int    *)malloc(3*sizeof(int   )*he);
  m->p = (int    *)malloc(1*sizeof(int   )*hn);

  memcpy(m->v, md->v, 2*sizeof(double)*hn);
  memcpy(m->e, md->e, 3*sizeof(int   )*he);
  memcpy(m->p, md->b,   sizeof(int   )*hn);
  return 0;
}

#ifdef CHECK_MESH
static int checkMesh(Mesh *m) 
{
  const int ne = m->d->ne;

  double *v = m->v;
  double *w;
  int    *e = m->e;

  double  x[6];
  double  f;
  int     v1, v2, v3;
  int     i, errs = 0;
  
  for (i = 0; i < ne; ++i) {
    v1 = e[0];
    v2 = e[1];
    v3 = e[2];
    e += 3;

    w = v + 2*v1;
    x[0] = w[0];
    x[3] = w[1];
  
    w = v + 2*v2;
    x[1] = w[0];
    x[4] = w[1];

    w = v + 2*v3;
    x[2] = w[0];
    x[5] = w[1];

    if (o_fcn(&f, x)) {
      if (f < -epsilon) {
	printf("Invalid element: %d: inverted.\n", i);
      }
      else {
	printf("Invalid element: %d: zero area.\n", i);
      } 
      ++errs;
    }
  }
  return errs;
}
#endif

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

#ifdef CHECK_MESH
  if (checkVertData(md)) {
    printf("Vertex check failure.\n");
    freeMeshData(&md);
    return -2;
  }
#endif

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

#ifdef CHECK_MESH
  if (checkMesh(*m)) {
    printf("Inverted element failure.\n");
    freeMesh(m);
    return -2;
  }
#endif

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
    fprintf(fp, "% 15.14e % 15.14e\n", v[0], v[1]);
    v += 2;
  }

  fprintf(fp, "%d\n", ne);

  e = m->e;
  for (i = 0; i < ne; ++i) {
    fprintf(fp, "%d %d %d\n", e[0], e[1], e[2]);
    e += 3;
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
    loc = 2*(*tp++);

    hv[0] = tv[loc  ];
    hv[1] = tv[loc+1];

    hv += 2;
  }
#else
  memcpy(hv, tv, 2*sizeof(double)*n);
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

  fprintf(fp, "\nvar x : 1 2 =\n");
  v = m->v;
  for (i = 0; i < nv; ++i) {
    fprintf(fp, "%6d % 15.14e % 15.14e\n", i+1, v[0], v[1]);
    v += 2;
  }

  fprintf(fp, ";\n\nparam TRIS : 1 2 3 =\n");
  e = m->e;
  for (i = 0; i < ne; ++i) {
    fprintf(fp, "%6d %6d %6d %6d\n", i+1, e[0]+1, e[1]+1, e[2]+1);
    e += 3;
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

  fprintf(fp, "\nvar x : 1 2 =\n");
  v = m->v;
  for (i = 0; i < nv; ++i) {
    fprintf(fp, "%6d % 15.14e % 15.14e\n", i+1, v[0], v[1]);
    v += 2;
  }

  fprintf(fp, ";\n\nparam TRIS : 1 2 3 =\n");
  e = m->e;
  for (i = 0; i < ne; ++i) {
    fprintf(fp, "%6d %6d %6d %6d\n", i+1, e[0]+1, e[1]+1, e[2]+1);
    e += 3;
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
    loc = 2*(*tp++);

    hv[0] = tv[loc  ];
    hv[1] = tv[loc+1];

    hv += 2;
  }
#else
  memcpy(hv, tv, 2*sizeof(double)*n);
#endif

  writeMeshOpt_AMPL(fname_opt, m);
  return writeMeshData_AMPL(fname, m->d);
}

int createMesh(double *verts, const int nv, 
	       int    *elems, const int ne, 
	       int    *fixed, const int nf, 
	       Mesh **m)
{
  MeshData *md = NULL;

  if (allocMeshData(&md)) {
    printf("Allocation failure.\n");
    return -1;
  }
  
  md->nv = nv;
  md->np = 1;
  md->no = 1;

  md->v = (double *)malloc(2*sizeof(double)*nv);
  memcpy(md->v, verts, 2*sizeof(double)*nv);

  md->ne = ne;
  md->e = (int *)malloc(3*sizeof(int)*ne);
  memcpy(md->e, elems, 3*sizeof(int)*ne);

  md->nf = nf;
  if (nf > 0) {
    md->f = (int *)malloc(sizeof(int)*nf);
    memcpy(md->f, fixed, sizeof(int)*nf);
  }

  if (checkMeshData(md)) {
    printf("Index check failure.\n");
    freeMeshData(&md);
    return -2;
  }

#ifdef CHECK_MESH
  if (checkVertData(md)) {
    printf("Vertex check failure.\n");
    freeMeshData(&md);
    return -2;
  }
#endif

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

#ifdef CHECK_MESH
  if (checkMesh(*m)) {
    printf("Inverted element failure.\n");
    freeMesh(m);
    return -2;
  }
#endif

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

