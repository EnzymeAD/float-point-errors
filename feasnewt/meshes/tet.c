#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mesh.h"
#include "fcn.h"

#define FACTOR_TOL 3.667e-11
#define PLANE_TOL  3.667e-11

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

#ifdef REFINE
#ifdef SPLINTER
  // This refinement creates splinters in the mesh and can only be
  // user a small number of times before the refined mesh becomes
  // impossible to optimize due to numerical errors
  if (!errs) {
    double *newv;
    int *newe;

    newv = (double *)malloc(3*sizeof(double)*(nv + ne));
    newe = (int *)malloc(16*sizeof(int)*ne);

    memcpy(newv, m->v, 3*sizeof(double)*nv);
    for (i = 0; i < ne; ++i) {
      newv[3*(nv+i)+0] = 0.25*(m->v[3*m->e[4*i+0]+0] + 
			       m->v[3*m->e[4*i+1]+0] +
			       m->v[3*m->e[4*i+2]+0] +
			       m->v[3*m->e[4*i+3]+0]);

      newv[3*(nv+i)+1] = 0.25*(m->v[3*m->e[4*i+0]+1] + 
			       m->v[3*m->e[4*i+1]+1] +
			       m->v[3*m->e[4*i+2]+1] +
			       m->v[3*m->e[4*i+3]+1]);

      newv[3*(nv+i)+2] = 0.25*(m->v[3*m->e[4*i+0]+2] + 
			       m->v[3*m->e[4*i+1]+2] +
			       m->v[3*m->e[4*i+2]+2] +
			       m->v[3*m->e[4*i+3]+2]);
    }

    free(m->v);
    m->v = newv;
    m->nv = nv+ne;

    for (i = 0; i < ne; ++i) {
      newe[16*i+ 0] = m->e[4*i+0];
      newe[16*i+ 1] = m->e[4*i+1];
      newe[16*i+ 2] = m->e[4*i+2];
      newe[16*i+ 3] = nv + i;

      newe[16*i+ 4] = m->e[4*i+0];
      newe[16*i+ 5] = m->e[4*i+3];
      newe[16*i+ 6] = m->e[4*i+1];
      newe[16*i+ 7] = nv + i;

      newe[16*i+ 8] = m->e[4*i+0];
      newe[16*i+ 9] = m->e[4*i+2];
      newe[16*i+10] = m->e[4*i+3];
      newe[16*i+11] = nv + i;

      newe[16*i+12] = m->e[4*i+1];
      newe[16*i+13] = m->e[4*i+3];
      newe[16*i+14] = m->e[4*i+2];
      newe[16*i+15] = nv + i;
    }

    free(m->e);
    m->e = newe;
    m->ne = 4*ne;
  }

#else
  // This version is a "uniform" refinement of the mesh.  Six vertices
  // are added per element, one vertex that bisects each edge.  The
  // tetrahedron is then split into eight pieces.  If we start with
  // regular tetrahedron, then the refined one will have four regular
  // ones on the corners and four right tetrahedron in the centers.
  // The latter are the reason why this is not a truly uniform
  // refinement operation.

  if (!errs) {
    double *newv;
    int *perm;
    int *newe;

    int nnv, vert, j;

    // We start off by adding the vertices to the mesh.  One vertex is
    // added per edge and most edges appear in many elements.  The upper 
    // bound on the number of new vertices 6*ne, but the final count will
    // be much lower once the duplicates are removed.  For the vertices,
    // we store four values, three doubles for the coordinates and an
    // integer for the initial position.  We reallocate the data once
    // we know the correct sizes.

    newv = (double *)malloc(4*sizeof(double)*(nv + 6*ne));
    perm = (int *)malloc(sizeof(int)*(nv + 6*ne));

    // There will be eight elements in the new discretization per element in
    // the old discretization; each element is a set of four indices.
    // We do not need to reallocate this data.

    newe = (int *)malloc(32*sizeof(int)*ne);

    // Copy the original vertex list into the allocated vertex list and
    // add the vertex numbers.

    for (i = 0; i < nv; ++i) {
      newv[4*i + 0] = m->v[3*i + 0];
      newv[4*i + 1] = m->v[3*i + 1];
      newv[4*i + 2] = m->v[3*i + 2];
      newv[4*i + 3] = i;
      perm[i] = -1;
    }

    // Add the new vertices and elements

    for (i = 0; i < ne; ++i) {
      // First edge: a: (0,1)
      newv[4*(nv + 6*i + 0) + 0] = 0.5*(m->v[3*m->e[4*i + 0] + 0] + 
			                m->v[3*m->e[4*i + 1] + 0]);
      newv[4*(nv + 6*i + 0) + 1] = 0.5*(m->v[3*m->e[4*i + 0] + 1] + 
			                m->v[3*m->e[4*i + 1] + 1]);
      newv[4*(nv + 6*i + 0) + 2] = 0.5*(m->v[3*m->e[4*i + 0] + 2] + 
			                m->v[3*m->e[4*i + 1] + 2]);
      newv[4*(nv + 6*i + 0) + 3] = nv + 6*i + 0;
      perm[nv + 6*i + 0] = -1;

      // Second edge: b: (0,2)
      newv[4*(nv + 6*i + 1) + 0] = 0.5*(m->v[3*m->e[4*i + 0] + 0] + 
			                m->v[3*m->e[4*i + 2] + 0]);
      newv[4*(nv + 6*i + 1) + 1] = 0.5*(m->v[3*m->e[4*i + 0] + 1] + 
			                m->v[3*m->e[4*i + 2] + 1]);
      newv[4*(nv + 6*i + 1) + 2] = 0.5*(m->v[3*m->e[4*i + 0] + 2] + 
			                m->v[3*m->e[4*i + 2] + 2]);
      newv[4*(nv + 6*i + 1) + 3] = nv + 6*i + 1;
      perm[nv + 6*i + 1] = -1;

      // Third edge: c: (0,3)
      newv[4*(nv + 6*i + 2) + 0] = 0.5*(m->v[3*m->e[4*i + 0] + 0] + 
			                m->v[3*m->e[4*i + 3] + 0]);
      newv[4*(nv + 6*i + 2) + 1] = 0.5*(m->v[3*m->e[4*i + 0] + 1] + 
			                m->v[3*m->e[4*i + 3] + 1]);
      newv[4*(nv + 6*i + 2) + 2] = 0.5*(m->v[3*m->e[4*i + 0] + 2] + 
			                m->v[3*m->e[4*i + 3] + 2]);
      newv[4*(nv + 6*i + 2) + 3] = nv + 6*i + 2;
      perm[nv + 6*i + 2] = -1;

      // Fourth edge: d: (1,2)
      newv[4*(nv + 6*i + 3) + 0] = 0.5*(m->v[3*m->e[4*i + 1] + 0] + 
			                m->v[3*m->e[4*i + 2] + 0]);
      newv[4*(nv + 6*i + 3) + 1] = 0.5*(m->v[3*m->e[4*i + 1] + 1] + 
			                m->v[3*m->e[4*i + 2] + 1]);
      newv[4*(nv + 6*i + 3) + 2] = 0.5*(m->v[3*m->e[4*i + 1] + 2] + 
			                m->v[3*m->e[4*i + 2] + 2]);
      newv[4*(nv + 6*i + 3) + 3] = nv + 6*i + 3;
      perm[nv + 6*i + 3] = -1;

      // Fifth edge: e: (1,3)
      newv[4*(nv + 6*i + 4) + 0] = 0.5*(m->v[3*m->e[4*i + 1] + 0] + 
			                m->v[3*m->e[4*i + 3] + 0]);
      newv[4*(nv + 6*i + 4) + 1] = 0.5*(m->v[3*m->e[4*i + 1] + 1] + 
			                m->v[3*m->e[4*i + 3] + 1]);
      newv[4*(nv + 6*i + 4) + 2] = 0.5*(m->v[3*m->e[4*i + 1] + 2] + 
			                m->v[3*m->e[4*i + 3] + 2]);
      newv[4*(nv + 6*i + 4) + 3] = nv + 6*i + 4;
      perm[nv + 6*i + 4] = -1;

      // Sixth edge: f: (2,3)
      newv[4*(nv + 6*i + 5) + 0] = 0.5*(m->v[3*m->e[4*i + 2] + 0] + 
			                m->v[3*m->e[4*i + 3] + 0]);
      newv[4*(nv + 6*i + 5) + 1] = 0.5*(m->v[3*m->e[4*i + 2] + 1] + 
			                m->v[3*m->e[4*i + 3] + 1]);
      newv[4*(nv + 6*i + 5) + 2] = 0.5*(m->v[3*m->e[4*i + 2] + 2] + 
			                m->v[3*m->e[4*i + 3] + 2]);
      newv[4*(nv + 6*i + 5) + 3] = nv + 6*i + 5;
      perm[nv + 6*i + 5] = -1;

      // First element (0, a, b, c)
      newe[32*i +  0] = m->e[4*i+0];
      newe[32*i +  1] = nv + 6*i + 0;
      newe[32*i +  2] = nv + 6*i + 1;
      newe[32*i +  3] = nv + 6*i + 2;

      // Second element (1, d, a, e)
      newe[32*i +  4] = m->e[4*i+1];
      newe[32*i +  5] = nv + 6*i + 3;
      newe[32*i +  6] = nv + 6*i + 0;
      newe[32*i +  7] = nv + 6*i + 4;

      // Third element (2, b, d, f)
      newe[32*i +  8] = m->e[4*i+2];
      newe[32*i +  9] = nv + 6*i + 1;
      newe[32*i + 10] = nv + 6*i + 3;
      newe[32*i + 11] = nv + 6*i + 5;

      // Fourth element (3, f, e, c)
      newe[32*i + 12] = m->e[4*i+3];
      newe[32*i + 13] = nv + 6*i + 5;
      newe[32*i + 14] = nv + 6*i + 4;
      newe[32*i + 15] = nv + 6*i + 2;

      // Fifth element (b, d, c, a)
      newe[32*i + 16] = nv + 6*i + 1;
      newe[32*i + 17] = nv + 6*i + 3;
      newe[32*i + 18] = nv + 6*i + 2;
      newe[32*i + 19] = nv + 6*i + 0;

      // Sixth element (c, d, e, a)
      newe[32*i + 20] = nv + 6*i + 2;
      newe[32*i + 21] = nv + 6*i + 3;
      newe[32*i + 22] = nv + 6*i + 4;
      newe[32*i + 23] = nv + 6*i + 0;

      // Seventh element (b, c, d, f)
      newe[32*i + 24] = nv + 6*i + 1;
      newe[32*i + 25] = nv + 6*i + 2;
      newe[32*i + 26] = nv + 6*i + 3;
      newe[32*i + 27] = nv + 6*i + 5;

      // Eighth element (c, e, d, f)
      newe[32*i + 28] = nv + 6*i + 2;
      newe[32*i + 29] = nv + 6*i + 4;
      newe[32*i + 30] = nv + 6*i + 3;
      newe[32*i + 31] = nv + 6*i + 5;
    }

    // Now sort the vertices by their coordinates and index number
    qsort(newv, nv + 6*ne, 4*sizeof(double), vert4);

    // Compute the permutation
    nnv = nv;
    for (i = 0; i < nv + 6*ne; ) {
      // Get the current vertex number
      vert = (int) newv[4*i + 3];

      if (vert < nv) {
	// Current vertex is an original vertex assumed unique
        perm[vert] = vert;
        assert(vert3(newv + 4*i, newv + 4*(i+1)));
        ++i;
      }
      else {
        // Current vertex is a new vertex; permute to next location
        // Get the vertices that are the same and permute them
        for (j = i; j < nv + 6*ne; ++j) {
          if (vert3(newv + 4*i, newv + 4*j)) {
	    // i and j not the same
            break;
          }
          perm[(int) newv[4*j + 3]] = nnv;
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
    for (; i < nv + 6*ne; ++i) {
      assert(perm[i] > 0 && perm[i] < nnv);
    }

    // Now permute the vertices
    free(m->v);
    m->v = (double *)malloc(3*sizeof(double)*nnv);
    m->nv = nnv;

    for (i = 0; i < nv + 6*ne; ++i) {
      vert = perm[(int) newv[4*i + 3]];

      m->v[3*vert + 0] = newv[4*i + 0];
      m->v[3*vert + 1] = newv[4*i + 1];
      m->v[3*vert + 2] = newv[4*i + 2];
    }

    free(newv);

    // Now permute the elements

    free(m->e);
    m->e = newe;
    m->ne = 8*ne;

    for (i = 0; i < 8*ne; ++i) {
      m->e[4*i + 0] = perm[m->e[4*i + 0]];
      m->e[4*i + 1] = perm[m->e[4*i + 1]];
      m->e[4*i + 2] = perm[m->e[4*i + 2]];
      m->e[4*i + 3] = perm[m->e[4*i + 3]];
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

#ifdef CHECK_MESH
static int checkVertData(MeshData *m)
{
#ifdef HEAL_MESH
  int *perm;
#endif

  double *v = m->v;
  int nv = m->nv;

  double *vertData, *vd;
  int i;

  /* Each record is:                 */
  /* | c1 | c2 | c3 | vert # | */

  vertData = (double *)malloc(4*sizeof(double)*nv);

  vd = vertData;
  for (i = 0; i < nv; ++i) {
    vd[0] = v[0];
    vd[1] = v[1];
    vd[2] = v[2];
    vd[3] = i;

    v += 3;
    vd += 4;
  }

  vd = vertData;
  qsort(vd, nv, 4*sizeof(double), vert4);

#ifdef HEAL_MESH
  perm = (int *)malloc(sizeof(int)*nv);
  for (i = 0; i < nv; ++i) {
    perm[i] = i;
  }
#endif

  for (i = 0; i < nv - 1; ++i) {
    if (!vert3(vd, vd+4)) {
      /* Two vertices match */
      
      printf("Same vertex: %d %d: %5.4e\n", (int) vd[3], (int) vd[7],
	     fabs(vd[0]-vd[4]) + fabs(vd[1]-vd[5]) + fabs(vd[2]-vd[6]));
      
#ifdef HEAL_MESH
      {
        const int o = (int) vd[3];
        const int n = (int) vd[7];

	if (o == perm[o]) {
	  perm[n] = o;
	}
	else {
	  perm[n] = perm[o];
	}
      }
#endif
    }
    vd += 4;
  }

#ifdef HEAL_MESH
  {
    const int ne = m->ne;
    int *ed = m->e;
    int j, k;

    for (j = 0; j < ne; ++j) {
      for (k = 0; k < 4; ++k) {
	if (ed[k] != perm[ed[k]]) {
	  ed[k] = perm[ed[k]];
	  printf("Fixing %d.%d\n", j, k);
	}
      }
      ed += 4;
    }
  }
  free(perm);
#endif
  
  free(vertData);
  return 0;
}
#endif

static int boundaryMeshData(MeshData *m)
{
  /* Set b[i] = 0 if vertex i unrestricted  */
  /* Set b[i] = 1 if vertex i on boundary   */
  /* Set b[i] = 2 if vertex i fixed by user */
  /* Set b[i] = 4 if vertex i unreferenced  */

#if defined(CHECK_MESH) || defined(SEGMENT_MESH)
  int *se;
#endif

  int *e = m->e;
  int *f = m->f;
  int *b;

  int *l1, *l2, *ls;

  int nv = m->nv;
  int ne = m->ne;
  int nf = m->nf;
  int nse = 0;

  int i, j, cnt = 0;
  int errs = 0;

  m->b = (int *)calloc(nv, sizeof(int));
  b = m->b;

  /* Counting argument for unreferenced vertices */
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
      b[0] = 4;
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

    /* Redefine b to be a permutation vector. */
    tv = 0;
    b = m->b;
    for (i = 0; i < nv; ++i) {
      if (0 == b[0]) {
	b[0] = tv++;
      }
      else {
	b[0] = -1;
      }
      ++b;
    }

    /* Apply the permutation vector to vertices */
    b = m->b;
    for (i = 0; i < nv; ++i) {
      if (b[0] >= 0) {
	j = 3*b[0];
	k = 3*i;

	v[j+0] = v[k+0];
	v[j+1] = v[k+1];
	v[j+2] = v[k+2];
      }
      ++b;
    }

    /* Apply the permutation vector to elements */
    b = m->b;
    e = m->e;
    for (i = 0; i < ne; ++i) {
      for (j = 0; j < 4; ++j) {
	e[j] = b[e[j]];
      }
      e += 4;
    }

    /* Apply the permutation vector to fixed list */
    f = m->f;
    for (i = 0; i < nf; ++i) {
      f[i] = b[f[i]];
    }

    /* Now all vertices are referenced */
    memset(b, 0, sizeof(int)*tv);

    /* Reset number of vertices. */
    m->nv = tv;
    nv = tv;
  }
#endif

  /* Fixed variables by the user */
  f = m->f;
  b = m->b;
  for (i = 0; i < nf; ++i) {
    b[f[0]] = 2;
    f += 1;
  }

  /* Finish with boundary vertices; compute and sort the faces */
  ls = (int *)malloc(   sizeof(int)*nv);         /* vertex starts */
  l1 = (int *)malloc(20*sizeof(int)*(ne+1));     /* list of faces */
  l2 = (int *)malloc(20*sizeof(int)*(ne+1));     /* list of faces */

  e = m->e;
  for (i = 0; i < ne; ++i) {
    l1[0] = e[0]; l1[1] = e[1]; l1[2] = e[2];
    l1[3] = i; l1[4] = 0;
    sort3(l1); l1 += 5;

    l1[0] = e[0]; l1[1] = e[1]; l1[2] = e[3];
    l1[3] = i; l1[4] = 1;
    sort3(l1); l1 += 5;

    l1[0] = e[0]; l1[1] = e[2]; l1[2] = e[3];
    l1[3] = i; l1[4] = 2;
    sort3(l1); l1 += 5;

    l1[0] = e[1]; l1[1] = e[2]; l1[2] = e[3];
    l1[3] = i; l1[4] = 3;
    sort3(l1); l1 += 5;

    e += 4;
  }

  l1[0] = -1;
  l1[1] = -1;
  l1[2] = -1;
  l1[3] = -1;
  l1[4] = -1;

  l1 -= 20*ne;

  /* Need to initialize final elements of l2 because these are used later. */
  l2[20*ne + 0] = -1;
  l2[20*ne + 1] = -1;
  l2[20*ne + 2] = -1;
  l2[20*ne + 3] = -1;
  l2[20*ne + 4] = -1;

  /* Now perform the radix sort */
  radix3(l2, l1, ls, 2, nv, 4*ne);
  radix3(l1, l2, ls, 1, nv, 4*ne);
  radix3(l2, l1, ls, 0, nv, 4*ne);

  free(ls);
  free(l1);			/* Results in l2   */
  l1 = l2;			/* Instead of copy */

  /* Check compatability between the faces */
  for (i = 0; i < 4*ne; ++i) {
    if ((l1[0] == l1[5]) && (l1[1] == l1[6]) && (l1[2] == l1[7])) {
      /* Face must be in only two elements! */
      cnt = 1;

#ifdef CHECK_MESH
      /* Check compatibility conditions for the faces             */
      /* Make sure one is clockwise and other is counter clockwise*/
      e = m->e + 4*l1[3];

      switch(l1[4]) {
      case 0:
        l1[0] = e[0]; l1[1] = e[1]; l1[2] = e[2];
        break;

      case 1:
        l1[0] = e[0]; l1[1] = e[3]; l1[2] = e[1];
        break;

      case 2:
        l1[0] = e[0]; l1[1] = e[2]; l1[2] = e[3];
        break;

      case 3:
        l1[0] = e[1]; l1[1] = e[3]; l1[2] = e[2];
        break;
      }
      face3(l1);

      e = m->e + 4*l1[8];
      switch(l1[9]) {
      case 0:
        l1[5] = e[0]; l1[6] = e[2]; l1[7] = e[1];
        break;

      case 1:
        l1[5] = e[0]; l1[6] = e[1]; l1[7] = e[3];
        break;

      case 2:
        l1[5] = e[0]; l1[6] = e[3]; l1[7] = e[2];
        break;

      case 3:
        l1[5] = e[1]; l1[6] = e[2]; l1[7] = e[3];
        break;
      }
      face3(l1 + 5);

      if ((l1[0] != l1[5]) || (l1[1] != l1[6])  || (l1[2] != l1[7])) {
        printf("Face Compatiblity: %d.%d and %d.%d\n",
               l1[3], l1[4], l1[8], l1[9]);
        printf("Face Compatiblity: (%d %d %d) and (%d %d %d)\n",
               l1[0], l1[1], l1[2], l1[5], l1[6], l1[7]);
	++errs;
      }

      sort3(l1);
      sort3(l1 + 5);
#endif

      /* Remove all replicated vertices */
      while ((i < 4*ne) &&
             (l1[0] == l1[5]) && (l1[1] == l1[6]) && (l1[2] == l1[7])) {
        l1 += 5;
        ++i;
	++cnt;
      }

      if (cnt > 2) {
	printf("Face Count Problem: %d.%d = %d\n", l1[3], l1[4], cnt);
	++errs;
      }
    }
    else {
      if (0 == b[l1[0]]) {
        b[l1[0]] = 1;
      }
      if (0 == b[l1[1]]) {
        b[l1[1]] = 1;
      }
      if (0 == b[l1[2]]) {
        b[l1[2]] = 1;
      }
      ++nse;
    }
    l1 += 5;
  }
  l1 -= 20*ne;

  m->nse = nse;

#if defined(CHECK_MESH) || defined(SEGMENT_MESH)
  /* Create surface mesh from the exposed faces */
  m->se = (int *)malloc(3*sizeof(int)*nse);
  se = m->se;

  for (i = 0; i < 4*ne; ++i) {
    if ((l1[0] == l1[5]) && (l1[1] == l1[6]) && (l1[2] == l1[7])) {
      /* Face in two elements; remove replicated faces */
      while ((i < 4*ne) &&
             (l1[0] == l1[5]) && (l1[1] == l1[6]) && (l1[2] == l1[7])) {
        l1 += 5;
        ++i;
      }
    }
    else {
      e = m->e + 4*l1[3];
      switch(l1[4]) {
      case 0:
        se[0] = e[0]; se[1] = e[1]; se[2] = e[2];
        break;

      case 1:
        se[0] = e[0]; se[1] = e[3]; se[2] = e[1];
        break;

      case 2:
        se[0] = e[0]; se[1] = e[2]; se[2] = e[3];
        break;

      case 3:
        se[0] = e[1]; se[1] = e[3]; se[2] = e[2];
        break;
      }
      se += 3;
    }
    l1 += 5;
  }
  l1 -= 20*ne;
#endif

  free(l1);
  return errs;
}

#if defined(CHECK_MESH) || defined(SEGMENT_MESH)
static int checkSurfData(MeshData *m, int hasBndry)
{
  const int nv = m->nv;
  const int nse = m->nse;

  int *se = m->se;
  int *l1, *l2, *ls;

  int i, cnt = 0;
  int errs = 0;

#ifdef SEGMENT_MESH
  int *si;

  /* Allocate surface incidence list */
  if (NULL == m->si) {
    m->si = (int *)calloc(3*nse, sizeof(int));
  }
  si = m->si;
#endif

  /* Sort edges in the surface mesh */
  ls = (int *)malloc(   sizeof(int)*nv);          /* vertex starts */
  l1 = (int *)malloc(12*sizeof(int)*(nse+1));     /* list of faces */
  l2 = (int *)malloc(12*sizeof(int)*(nse+1));     /* list of faces */

  for (i = 0; i < nse; ++i) {
    l1[0] = se[0]; l1[1] = se[1]; 
    l1[2] = i; l1[3] = 0;
    sort2(l1); l1 += 4;

    l1[0] = se[0]; l1[1] = se[2]; 
    l1[2] = i; l1[3] = 1;
    sort2(l1); l1 += 4;

    l1[0] = se[1]; l1[1] = se[2]; 
    l1[2] = i; l1[3] = 2;
    sort2(l1); l1 += 4;

    se += 3;
  }

  l1[0] = -1;
  l1[1] = -1;
  l1[2] = -1;
  l1[3] = -1;

  l1 -= 12*nse;

  /* Now perform the radix sort */
  radix2(l2, l1, ls, 1, nv, 3*nse);
  radix2(l1, l2, ls, 0, nv, 3*nse);

  free(ls);
  free(l2);

  /* Check compatability between the faces */
  for (i = 0; i < 3*nse; ++i) {
    if ((l1[0] == l1[4]) && (l1[1] == l1[5])) {
      /* Edge must be in only two elements! */
      cnt = 1;

      /* Check compatibility conditions for the edges             */
      /* Make sure one is clockwise and other is counter clockwise*/
      se = m->se + 3*l1[2];

      switch(l1[3]) {
      case 0:
        l1[0] = se[0]; l1[1] = se[1];
        break;

      case 1:
        l1[0] = se[2]; l1[1] = se[0];
        break;

      case 2:
        l1[0] = se[1]; l1[1] = se[2];
        break;
      }

      se = m->se + 3*l1[6];
      switch(l1[7]) {
      case 0:
        l1[4] = se[1]; l1[5] = se[0];
        break;

      case 1:
        l1[4] = se[0]; l1[5] = se[2];
        break;

      case 2:
        l1[4] = se[2]; l1[5] = se[1];
        break;
      }

      if ((l1[0] != l1[4]) || (l1[1] != l1[5])) {
        printf("Edge Compatiblity: %d.%d and %d.%d\n",
               l1[2], l1[3], l1[6], l1[7]);
	/* printf("  Vertices %d, %d\n", l1[0], l1[1]); */
	/* printf("  Vertices %d, %d\n", l1[4], l1[5]); */
	++errs;
      }
      
      sort2(l1);
      sort2(l1 + 4);

#ifdef SEGMENT_MESH
      /* Update surface incidence list */
      si[3*l1[2] + l1[3]] = l1[6];
      si[3*l1[6] + l1[7]] = l1[2];
#endif

      /* Remove all replicated vertices */
      while ((i < 3*nse) && (l1[0] == l1[4]) && (l1[1] == l1[5])) {
        l1 += 4;
        ++i;
	++cnt;
      }

      if (cnt > 2) {
	printf("Edge Count Problem: %d.%d = %d\n", l1[2], l1[3], cnt);
	printf("  Vertices %d, %d\n", l1[0], l1[1]);
	// ++errs;
      }
    }
    else {
      if (!hasBndry) {
	printf("Closed Surface Problem: %d.%d\n", l1[2], l1[3]);
	++errs;
      }
    }
    l1 += 4;
  }

  l1 -= 12*nse;
  free(l1);
  return errs;
}
#endif

#ifdef SEGMENT_MESH
static int segmentSurfData(MeshData *m)
{
  const int *si = m->si;

  const int nv = m->nv;
  const int nse = m->nse;

  int *ss = m->ss;
  int *q;

  int nss = 0;

  int qs, qe;
  int i, idx, loc;

  if (nv <= 0) {
    return 0;
  }

  q = (int *)malloc(nse*sizeof(int));
  
  idx = 0;
  while (idx >= 0) {
    ++nss;
    ss[idx] = nss;

    qs = 0;
    qe = 0;

    loc = 3*idx;
    for (i = 0; i < 3; ++i) {
      if (!ss[si[loc+i]]) {
	/* Mark as being in the segment and add index to queue */
	ss[si[loc+i]] = nss;
	q[qe] = si[loc+i];
	qe = (qe+1) % nse;
      }
    }

    while (qs != qe) {
      loc = 3*q[qs];
      qs = (qs+1) % nse;

      for (i = 0; i < 3; ++i) {
	if (!ss[si[loc+i]]) {
	  /* Mark as being in the segment and add index to queue */
	  ss[si[loc+i]] = nss;
	  q[qe] = si[loc+i];
	  qe = (qe+1) % nse;
	}
      }
    }
    
    for (i = idx+1; i < nse; ++i) {
      if (!ss[i]) {
	idx = i;
	break;
      }
      idx = -1;
    }
  }

  m->nss = nss;
  free(q);
  return 0;
}

static void fillQ(long double Q[7][7], double coeff[3][3])
{
  Q[0][0] = 1;
  Q[0][1] = 0;
  Q[0][2] = 0;
  Q[0][3] = 0;
  Q[0][4] = -coeff[0][0];
  Q[0][5] = -coeff[1][0];
  Q[0][6] = -coeff[2][0];
  
  Q[1][0] = 0;
  Q[1][1] = 1;
  Q[1][2] = 0;
  Q[1][3] = 0;
  Q[1][4] = -coeff[0][1];
  Q[1][5] = -coeff[1][1];
  Q[1][6] = -coeff[2][1];
  
  Q[2][0] = 0;
  Q[2][1] = 0;
  Q[2][2] = 1;
  Q[2][3] = 0;
  Q[2][4] = -coeff[0][2];
  Q[2][5] = -coeff[1][2];
  Q[2][6] = -coeff[2][2];
  
  Q[3][0] = 0;
  Q[3][1] = 0;
  Q[3][2] = 0;
  Q[3][3] = 1;
  Q[3][4] = 1;
  Q[3][5] = 1;
  Q[3][6] = 1;
  
  Q[4][0] = coeff[0][0];
  Q[4][1] = coeff[0][1];
  Q[4][2] = coeff[0][2];
  Q[4][3] = -1;
  Q[4][4] = 0;
  Q[4][5] = 0;
  Q[4][6] = 0;
  
  Q[5][0] = coeff[1][0];
  Q[5][1] = coeff[1][1];
  Q[5][2] = coeff[1][2];
  Q[5][3] = -1;
  Q[5][4] = 0;
  Q[5][5] = 0;
  Q[5][6] = 0;
  
  Q[6][0] = coeff[2][0];
  Q[6][1] = coeff[2][1];
  Q[6][2] = coeff[2][2];
  Q[6][3] = -1;
  Q[6][4] = 0;
  Q[6][5] = 0;
  Q[6][6] = 0;
  return;
}

static void fillRhs(long double rhs[7])
{
  rhs[0] = 1;
  rhs[1] = 1;
  rhs[2] = 1;
  rhs[3] = 1;
  rhs[4] = 0;
  rhs[5] = 0;
  rhs[6] = 0;
  return;
}

static int factor(long double Q[7][7], long double R[7][7], 
		  long double v[7], long double w[7]) 
{
  /* Householder QR factorization */
  long double alpha, beta;
  int i, j, k;

  for (i = 0; i < 7; ++i) {
    for (j = 0; j < 7; ++j) {
      R[i][j] = Q[i][j];
    }
  }

  for (i = 0; i < 7; ++i) {
    /* Compute Householder */
    alpha = 0.0;
    for (j = i; j < 7; ++j) {
      alpha += R[j][i]*R[j][i];
    }
    alpha = sqrtl(alpha);

    if (R[i][i] >= 0) {
      alpha = R[i][i] + alpha;
    }
    else {
      alpha = R[i][i] - alpha;
    }

    v[i] = 1.0;
    beta = 1.0;
    for (j = i+1; j < 7; ++j) {
      v[j] = R[j][i] / alpha;
      beta += v[j]*v[j];
    }
    beta = -2.0 / beta;

    /* Apply rotation */
    for (j = i; j < 7; ++j) {
      w[j] = 0;
      for (k = i; k < 7; ++k) {
        w[j] += R[k][j]*v[k];
      }
      w[j] *= beta;
    }

    for (j = i; j < 7; ++j) {
      for (k = i; k < 7; ++k) {
        R[j][k] += v[j]*w[k];
      }
    }

    if (fabsl(R[i][i]) < FACTOR_TOL) {
      printf("Not invertible: %d\n", i);
      return -1;
    }

    /* Store vector */
    for (j = i+1; j < 7; ++j) {
      R[j][i] = v[j];
    }
  }

  for (i = 0; i < 7; ++i) {
    for (j = 0; j < 7; ++j) {
      Q[i][j] = 0.0;
    }
    Q[i][i] = 1.0;
  }

  for (i = 6; i >= 0; --i) {
    v[i] = 1.0;
    beta = 1.0;
    for (j = i+1; j < 7; ++j) {
      v[j] = R[j][i];
      beta += v[j]*v[j];
    }
    beta = -2.0 / beta;

    /* Apply rotation */
    for (j = i; j < 7; ++j) {
      w[j] = 0;
      for (k = i; k < 7; ++k) {
        w[j] += Q[k][j]*v[k];
      }
      w[j] *= beta;
    }

    for (j = i; j < 7; ++j) {
      for (k = i; k < 7; ++k) {
        Q[j][k] += v[j]*w[k];
      }
    }
  }
  return 0;
}

static void solve(long double Q[7][7], long double R[7][7], 
		  long double r[7], long double x[7]) 
{
  int i, j;

  for (i = 0; i < 7; ++i) {
    x[i] = 0;
    for (j = 0; j < 7; ++j) {
      x[i] += Q[j][i]*r[j];
    }
  }

  for (i = 6; i >= 0; --i) {
    x[i] /= R[i][i];
    for (j = i-1; j >= 0; --j) {
      x[j] -= R[j][i]*x[i];
    }
  }
  return;
}

static int planarSurfData(MeshData *m)
{
  const int *si = m->si;

  const int nv = m->nv;
  const int nse = m->nse;
  const int nss = m->nss;

  double *v;
  int *se;
  int *ss = m->ss;

  int *pe;			/* Elements defining plane           */
  int *ps;			/* Start of plane i in pe            */
  int *pb;			/* Boundary vertices of the plane    */
  int *pbc;			/* Number of times vertex appears in */
                                /* boundary of the planes            */

  int *q;

  double coeff[3][3];
  long double Q[7][7];
  long double R[7][7];
  long double r[7];
  long double w1[7];
  long double w2[7];

  int nps = 0;

  int cs;

  int qs, qe, st, en;
  int i, j, idx, loc;
  int nnse, nnps, npe;

  if (nv <= 0) {
    return 0;
  }

  /* Allocate pe and ps */
  pe = (int *)malloc(nse*sizeof(int));
  ps = (int *)malloc((nse+1)*sizeof(int));
  ps[0] = 0;

  q = (int *)malloc(nse*sizeof(int));

  /* Compute planar segments */
  for (cs = 1; cs <= nss; ++cs) {
    idx = -1;
    while (1) {
      /* Get an element in current segment that is not evaluated */
      for (i = idx+1; i < nse; ++i) {
	if (cs == ss[i]) {
	  idx = i;
	  break;
	}
      }

      /* No more elements in current segment to visit */
      if (nse == i) {
	break;
      }

      /* Compute a plane passing through the element */
      se = m->se + 3*idx;
      for (i = 0; i < 3; ++i) {
	v = m->v + 3*se[i];
	coeff[i][0] = v[0];
	coeff[i][1] = v[1];
	coeff[i][2] = v[2];
      }

      fillQ(Q, coeff);
      if (factor(Q, R, r, w1)) {
	printf("Matrix not invertible.\n");
	exit(-1);
      }

      fillRhs(r);
      solve(Q, R, r, w1);
      
      /* Add the current element to the plane */
      ++nps;
      ss[idx] = nss + nps;

      ps[nps] = ps[nps-1];
      pe[ps[nps]++] = idx;

      qs = 0;
      qe = 0;

      loc = 3*idx;
      for (i = 0; i < 3; ++i) {
	if (cs == ss[si[loc+i]]) {
	  /* This incident element has not been added */
	  /* Check to see if it is in the current plane */
	  se = m->se + 3*si[loc+i];
	  for (j = 0; j < 3; ++j) {
	    v = m->v + 3*se[j];
	    coeff[j][0] = v[0];
	    coeff[j][1] = v[1];
	    coeff[j][2] = v[2];
	  }

	  fillQ(Q, coeff);
	  if (factor(Q, R, r, w2)) {
	    printf("Matrix not invertible.\n");
	    exit(-1);
	  }
	  
	  fillRhs(r);
	  solve(Q, R, r, w2);
      
	  if (sqrtl((w2[0]-w1[0])*(w2[0]-w1[0]) +
		    (w2[1]-w1[1])*(w2[1]-w1[1]) +
		    (w2[2]-w1[2])*(w2[2]-w1[2]) +
		    (w2[3]-w1[3])*(w2[3]-w1[3])) < PLANE_TOL) {
	    /* Element is in the plane */
	    ss[si[loc+i]] = nss + nps;
	    q[qe] = si[loc+i];
	    qe = (qe+1) % nse;

	    pe[ps[nps]++] = si[loc+i];
	  }
	}
      }

      while (qs != qe) {
	loc = 3*q[qs];
	qs = (qs+1) % nse;

	for (i = 0; i < 3; ++i) {
	  if (cs == ss[si[loc+i]]) {
	    /* This incident element has not been added */
	    /* Check to see if it is in the current plane */
	    se = m->se + 3*si[loc+i];
	    for (j = 0; j < 3; ++j) {
	      v = m->v + 3*se[j];
	      coeff[j][0] = v[0];
	      coeff[j][1] = v[1];
	      coeff[j][2] = v[2];
	    }
	    
	    fillQ(Q, coeff);
	    if (factor(Q, R, r, w2)) {
	      printf("Matrix not invertible.\n");
	      exit(-1);
	    }
	    
	    fillRhs(r);
	    solve(Q, R, r, w2);
	    
	    if (sqrtl((w2[0]-w1[0])*(w2[0]-w1[0]) +
		      (w2[1]-w1[1])*(w2[1]-w1[1]) +
		      (w2[2]-w1[2])*(w2[2]-w1[2]) +
		      (w2[3]-w1[3])*(w2[3]-w1[3])) < PLANE_TOL) {
	      /* Element is in the plane */
	      ss[si[loc+i]] = nss + nps;
	      q[qe] = si[loc+i];
	      qe = (qe+1) % nse;

	      pe[ps[nps]++] = si[loc+i];
	    }
	  }
	}
      }
    }
  }

  free(q);

  pb  = (int *)calloc(nv, sizeof(int));
  pbc = (int *)calloc(nv, sizeof(int));

  /* Count number of times each vertex appears in a plane */
  for (i = 0; i < nps; ++i) {
    for (j = ps[i]; j < ps[i+1]; ++j) {
      se = m->se + 3*pe[j];
      
      ++pb[se[0]];
      ++pb[se[1]];
      ++pb[se[2]];
    }

    for (j = ps[i]; j < ps[i+1]; ++j) {
      se = m->se + 3*pe[j];
      
      if (pb[se[0]]) {
	++pbc[se[0]];
	pb[se[0]] = 0;
      }

      if (pb[se[1]]) {
	++pbc[se[1]];
	pb[se[1]] = 0;
      }

      if (pb[se[2]]) {
	++pbc[se[2]];
	pb[se[2]] = 0;
      }
    }
  }
  free(pb);

  /* Appears in no planes    -> interior vertex          */
  /* Appears in one plane    -> just mapped to the plane */
  /* Appears in two planes   -> one an edge of the plane */
  /* Appears in three planes -> a fixed point            */
  /* Now go through the planes and reduce the segments   */
  /* If all vertices are fixed remove the plane...       */
  nnps = 0;
  npe = 0;
  en = ps[0];
  for (i = 0; i < nps; ++i) {
    st = en;
    en = ps[i+1];

    for (j = st; j < en; ++j) {
      se = m->se + 3*pe[j];
      if ((pbc[se[0]] <= 2) || (pbc[se[1]] <= 2) || (pbc[se[2]] <= 2)) {
	break;
      }
    }

    if (j == ps[i+1]) {
      /* All are fixed; remove the plane description */
      continue;
    }

    ++nnps;			/* Increment number of planes */
    ps[nnps] = ps[nnps-1];      /* Patch start vector    */
    
    for (; j < en; ++j) {	/* Continue where left off */
      se = m->se + 3*pe[j];
      if ((pbc[se[0]] <= 2) || (pbc[se[1]] <= 2) || (pbc[se[2]] <= 2)) {
	pe[npe++] = pe[j];	/* Remap element description */
	++ps[nnps];
      }
    }
  }
  nps = nnps;			/* Fix number of planes */

  m->npe = npe;
  m->nps = nps;
  if (nps) {
    /* Allocate the data for the planar elements */

    m->ps = ps;
    m->pd = (double *)malloc(4*sizeof(double)*npe);
    m->pe = (int *)malloc(3*sizeof(int)*npe);

    for (i = 0; i < npe; ++i) {
      se = m->se + 3*pe[i];
      m->pe[3*i+0] = se[0];
      m->pe[3*i+1] = se[1];
      m->pe[3*i+2] = se[2];
    }
    free(pe);

    for (i = 0; i < nps; ++i) {
      /* Fill in plane data */
      pe = m->pe + 3*m->ps[i];
      for (j = 0; j < 3; ++j) {
	v = m->v + 3*pe[j];
	coeff[j][0] = v[0];
	coeff[j][1] = v[1];
	coeff[j][2] = v[2];
      }

      fillQ(Q, coeff);
      factor(Q, R, r, w2);
      fillRhs(r);
      solve(Q, R, r, w2);

      m->pd[4*i+0] = (double) w2[0];
      m->pd[4*i+1] = (double) w2[1];
      m->pd[4*i+2] = (double) w2[2];
      m->pd[4*i+3] = (double) w2[3];
    }

    /* Collapse nse, se, etc */
    nnse = 0;
    se = m->se;
    for (i = 0; i < nse; ++i) {
      if ((pbc[se[0]] > 2) && (pbc[se[1]] > 2) && (pbc[se[2]] > 2)) {
	m->se[nnse+0] = se[0];
	m->se[nnse+1] = se[1];
	m->se[nnse+2] = se[2];
	nnse += 3;
      }
      se += 3;
    }
    m->nse = nnse / 3;
    m->nss = 0;

    memset(m->ss, 0, m->nse);
    memset(m->si, 0, 3*m->nse);
  }
  else {
    free(ps);
    free(pe);
  }

  printf("Planar  Segments: %d elements: %d\n", m->nps, m->npe);
  free(pbc);
  return 0;
}
#endif

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

#ifdef CHECK_MESH
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
  
#ifdef USE_WEIGHT
  double weight[6] = {1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
  double target[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
#endif

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

#ifndef USE_WEIGHT
    if (o_fcn(&f, x)) {
      if (f < -epsilon) {
	printf("Invalid element: %d: inverted.\n", i);
      } 
      else {
	printf("Invalid element: %d: zero area: %5.4e.\n", i, f);
      }
      ++errs;
    }
#else
    if (o_fcn(&f, x, weight, target)) {
      if (f < -epsilon) {
	printf("Invalid element: %d: inverted.\n", i);
      } 
      else {
	printf("Invalid element: %d: zero area: %5.4e.\n", i, f);
      }
      ++errs;
    }
#endif

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

#if defined(CHECK_MESH) || defined(SEGMENT_MESH)
  if (checkSurfData(md, 0)) {
    printf("Surface check failure.\n");
    freeMeshData(&md);
    return -2;
  }
#endif

#ifdef SEGMENT_MESH
  md->ss = (int *)calloc(md->nse, sizeof(int));
  if (segmentSurfData(md)) {
    printf("Surface segmentation failure.\n");
    freeMeshData(&md);
    return -2;
  }
  printf("Surface Segments: %d elements: %d\n", md->nss, md->nse);

  if (planarSurfData(md)) {
    printf("Planar segmentation failure.\n");
    freeMeshData(&md);
    return -2;
  }

  checkSurfData(md, 1);		/* Hack to get incidences */
  if (segmentSurfData(md)) {
    printf("Surface segmentation failure.\n");
    freeMeshData(&md);
    return -2;
  }
  printf("Surface Segments: %d elements: %d\n", md->nss, md->nse);
#endif

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
    fprintf(fp, "% 15.14e % 15.14e % 15.14e\n", v[0], v[1], v[2]);
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
    fprintf(fp, "%6d % 15.14e % 15.14e % 15.14e\n", i+1, v[0], v[1], v[2]);
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
    fprintf(fp, "%6d % 15.14e % 15.14e % 15.14e\n", i+1, v[0], v[1], v[2]);
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

  md->v = (double *)malloc(3*sizeof(double)*nv);
  memcpy(md->v, verts, 3*sizeof(double)*nv);

  md->ne = ne;
  md->e = (int *)malloc(4*sizeof(int)*ne);
  memcpy(md->e, elems, 4*sizeof(int)*ne);

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

#if defined(LIBRARY)
#include <sys/time.h>
#include <sys/resource.h>

#include "opt.h"

#ifndef MICROSEC
#  define MICROSEC 1000000
#endif

void meshOpt(int max_iter, double conv_tol,
	     int nv, double *v, int ne, int *e, int nf, int *f)
{
  Mesh *m = NULL;
  MeshData *md = NULL;
  double *hv, *tv;

#ifdef REORDER
  int *tp;
#endif

  double m0, m1, m2;
  double s0, s1, s2;
  double t1, t2;

  struct rusage r0, r1, r2;

  int i;

#ifdef REORDER
  int loc;
#endif

  if (nv <= 0) {
    fprintf(stderr, "Invalid number of vertices: %d\n", nv);
    exit(-1);
  }

  if (ne <= 0) {
    fprintf(stderr, "Invalid number of elements: %d\n", ne);
    exit(-1);
  }

  if (nf <  0) {
    fprintf(stderr, "Invalid number of fixed vertices: %d\n", nf);
    exit(-1);
  }

  if (NULL == v) {
    fprintf(stderr, "Vertices undefined: %d\n", nv);
    exit(-1);
  }

  if (NULL == e) {
    fprintf(stderr, "Elements undefined: %d\n", nv);
    exit(-1);
  }

  if ((nf > 0) && (NULL == f)) {
    fprintf(stderr, "Fixed vertices undefined: %d\n", nv);
    exit(-1);
  }
  
  getrusage(RUSAGE_SELF, &r0);
  if (allocMeshData(&md)) {
    printf("Allocation failure.\n");
    exit(-1);
  }

  md->nv = nv;
  md->np = 1;
  md->no = 1;

  md->v = v;

  md->ne = ne;
  md->e = e;

  md->nf = nf;
  if (nf > 0) {
    md->f = f;
  }

  if (checkMeshData(md)) {
    printf("Index check failure.\n");
    freeMeshData(&md);
    exit(-2);
  }

#ifdef CHECK_MESH
  if (checkVertData(md)) {
    printf("Vertex check failure.\n");
    freeMeshData(&md);
    exit(-2);
  }
#endif

  if (boundaryMeshData(md)) {
    printf("Boundary calculation failure.\n");
    freeMeshData(&md);
    exit(-2);
  }

#if defined(CHECK_MESH) || defined(SEGMENT_MESH)
  if (checkSurfData(md, 0)) {
    printf("Surface check failure.\n");
    freeMeshData(&md);
    exit(-2);
  }
#endif

#ifdef SEGMENT_MESH
  if (segmentSurfData(md)) {
    printf("Segmentation segmentation failure.\n");
    freeMeshData(&md);
    exit(-2);
  }
#endif

  if (allocMesh(&m)) {
    printf("Mesh allocation failure.\n");
    freeMeshData(&md);
    exit(-1);
  }

  m->d = md;

  if (computeMesh(m)) {
    printf("Mesh construction failure.\n");
    freeMesh(&m);
    exit(-2);
  }

#ifdef CHECK_MESH
  if (checkMesh(m)) {
    printf("Inverted element failure.\n");
    freeMesh(&m);
    exit(-2);
  }
#endif

#ifdef REORDER
  if (reorderMesh(m)) {
    printf("Reordering failure.\n");
    freeMesh(&m);
    exit(-2);
  }
#endif

  if (finishMesh(m)) {
    printf("Mesh initialization failure.\n");
    freeMesh(&m);
    exit(-2);
  }

  getrusage(RUSAGE_SELF, &r1);

  m1 = (double) r1.ru_utime.tv_usec;
  m2 = (double) r0.ru_utime.tv_usec;
  m0 = m1 - m2;
    
  s1 = (double) r1.ru_utime.tv_sec;
  s2 = (double) r0.ru_utime.tv_sec;
  s0 = s1 - s2;

  t1 = s0 + m0 / MICROSEC;

  m1 = (double) r1.ru_stime.tv_usec;
  m2 = (double) r0.ru_stime.tv_usec;
  m0 = m1 - m2;

  s1 = (double) r1.ru_stime.tv_sec;
  s2 = (double) r0.ru_stime.tv_sec;
  s0 = s1 - s2;

  t2 = s0 + m0 / MICROSEC;

  printf("Read: %5.4e System: %5.4e Total: %5.4e\n", t1, t2, t1+t2);

#ifndef REFINE
  getrusage(RUSAGE_SELF, &r1);
  optMesh(m, max_iter, conv_tol, 2); 
  getrusage(RUSAGE_SELF, &r2);

  m1 = (double) r2.ru_utime.tv_usec;
  m2 = (double) r1.ru_utime.tv_usec;
  m0 = m1 - m2;
    
  s1 = (double) r2.ru_utime.tv_sec;
  s2 = (double) r1.ru_utime.tv_sec;
  s0 = s1 - s2;

  t1 = s0 + m0 / MICROSEC;

  m1 = (double) r2.ru_stime.tv_usec;
  m2 = (double) r1.ru_stime.tv_usec;
  m0 = m1 - m2;

  s1 = (double) r2.ru_stime.tv_sec;
  s2 = (double) r1.ru_stime.tv_sec;
  s0 = s1 - s2;

  t2 = s0 + m0 / MICROSEC;

  printf("Optimize: User: %5.4e System: %5.4e Total: %5.4e\n", t1, t2, t1+t2);
#endif

  getrusage(RUSAGE_SELF, &r1);

  hv = md->v;
  tv = m->v;

#ifdef REORDER
  tp = m->per;

  for (i = 0; i < nv; ++i) {
    loc = 3*(*tp++);

    hv[0] = tv[loc  ];
    hv[1] = tv[loc+1];
    hv[2] = tv[loc+2];

    hv += 3;
  }
#else
  memcpy(hv, tv, 3*sizeof(double)*nv);
#endif

  md->v = NULL;
  md->e = NULL;
  md->f = NULL;
  freeMesh(&m);
  
  getrusage(RUSAGE_SELF, &r2);
  m1 = (double) r2.ru_utime.tv_usec;
  m2 = (double) r0.ru_utime.tv_usec;
  m0 = m1 - m2;
    
  s1 = (double) r2.ru_utime.tv_sec;
  s2 = (double) r0.ru_utime.tv_sec;
  s0 = s1 - s2;

  t1 = s0 + m0 / MICROSEC;

  m1 = (double) r2.ru_stime.tv_usec;
  m2 = (double) r0.ru_stime.tv_usec;
  m0 = m1 - m2;

  s1 = (double) r2.ru_stime.tv_sec;
  s2 = (double) r0.ru_stime.tv_sec;
  s0 = s1 - s2;

  t2 = s0 + m0 / MICROSEC;

  printf("Total: User: %5.4e System: %5.4e Total: %5.4e\n", t1, t2, t1+t2);
  return;
}
#endif
