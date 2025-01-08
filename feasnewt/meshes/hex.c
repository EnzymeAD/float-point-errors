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
  m->e = (int *)malloc(8*sizeof(int)*ne);

  /***************************************************************************/
  /* Read eight indices for each vertex.                                     */
  /***************************************************************************/

#ifdef SWEEP
  e = m->e;
  for (i = 0; i < ne; ++i) {
    buf[0] = '\0';
    fgets(buf, 1024, fp); ++line;
    if (8 != sscanf(buf, "%d%d%d%d%d%d%d%d",
                    e, e+1, e+3, e+2, e+4, e+5, e+7, e+6)) {
      printf("Read error on line %d: %s\n", line,
	     "expecting eight integer vertex indices");
      errs = 1;
    }
    e += 8;
  }
#else
  /* Hypercube ordering */
  e = m->e;
  for (i = 0; i < ne; ++i) {
    buf[0] = '\0';
    fgets(buf, 1024, fp); ++line;
    if (8 != sscanf(buf, "%d%d%d%d%d%d%d%d",
                    e, e+1, e+2, e+3, e+4, e+5, e+6, e+7)) {
      printf("Read error on line %d: %s\n", line,
	     "expecting eight integer vertex indices");
      errs = 1;
    }
    e += 8;
  }
#endif

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

  // Not that the method implemented can cause inverted elements because
  // hexahedral elements may not be convex.

  if (!errs) {
    double *newv;
    int *newe;

    double cent[3];

    newv = (double *)malloc(3*sizeof(double)*(nv + 8*ne));
    newe = (int *)malloc(56*sizeof(int)*ne);

    memcpy(newv, m->v, 3*sizeof(double)*nv);
    for (i = 0; i < ne; ++i) {
      cent[0] = 0.125*(m->v[3*m->e[8*i + 0] + 0] +
                       m->v[3*m->e[8*i + 1] + 0] +
                       m->v[3*m->e[8*i + 2] + 0] +
                       m->v[3*m->e[8*i + 3] + 0] +
                       m->v[3*m->e[8*i + 4] + 0] +
                       m->v[3*m->e[8*i + 5] + 0] +
                       m->v[3*m->e[8*i + 6] + 0] +
                       m->v[3*m->e[8*i + 7] + 0]);

      cent[1] = 0.125*(m->v[3*m->e[8*i + 0] + 1] +
                       m->v[3*m->e[8*i + 1] + 1] +
                       m->v[3*m->e[8*i + 2] + 1] +
                       m->v[3*m->e[8*i + 3] + 1] +
                       m->v[3*m->e[8*i + 4] + 1] +
                       m->v[3*m->e[8*i + 5] + 1] +
                       m->v[3*m->e[8*i + 6] + 1] +
                       m->v[3*m->e[8*i + 7] + 1]);

      cent[2] = 0.125*(m->v[3*m->e[8*i + 0] + 2] +
                       m->v[3*m->e[8*i + 1] + 2] +
                       m->v[3*m->e[8*i + 2] + 2] +
                       m->v[3*m->e[8*i + 3] + 2] +
                       m->v[3*m->e[8*i + 4] + 2] +
                       m->v[3*m->e[8*i + 5] + 2] +
                       m->v[3*m->e[8*i + 6] + 2] +
                       m->v[3*m->e[8*i + 7] + 2]);

      // New vertex a:  (v0 + center) / 2
      newv[3*(nv+8*i+0)+0] = (double) 0.5*(m->v[3*m->e[8*i + 0] + 0] + cent[0]);
      newv[3*(nv+8*i+0)+1] = (double) 0.5*(m->v[3*m->e[8*i + 0] + 1] + cent[1]);
      newv[3*(nv+8*i+0)+2] = (double) 0.5*(m->v[3*m->e[8*i + 0] + 2] + cent[2]);

      // New vertex b:  (v1 + center) / 2
      newv[3*(nv+8*i+1)+0] = (double) 0.5*(m->v[3*m->e[8*i + 1] + 0] + cent[0]);
      newv[3*(nv+8*i+1)+1] = (double) 0.5*(m->v[3*m->e[8*i + 1] + 1] + cent[1]);
      newv[3*(nv+8*i+1)+2] = (double) 0.5*(m->v[3*m->e[8*i + 1] + 2] + cent[2]);

      // New vertex c:  (v2 + center) / 2
      newv[3*(nv+8*i+2)+0] = (double) 0.5*(m->v[3*m->e[8*i + 2] + 0] + cent[0]);
      newv[3*(nv+8*i+2)+1] = (double) 0.5*(m->v[3*m->e[8*i + 2] + 1] + cent[1]);
      newv[3*(nv+8*i+2)+2] = (double) 0.5*(m->v[3*m->e[8*i + 2] + 2] + cent[2]);

      // New vertex d:  (v3 + center) / 2
      newv[3*(nv+8*i+3)+0] = (double) 0.5*(m->v[3*m->e[8*i + 3] + 0] + cent[0]);
      newv[3*(nv+8*i+3)+1] = (double) 0.5*(m->v[3*m->e[8*i + 3] + 1] + cent[1]);
      newv[3*(nv+8*i+3)+2] = (double) 0.5*(m->v[3*m->e[8*i + 3] + 2] + cent[2]);

      // New vertex e:  (v4 + center) / 2
      newv[3*(nv+8*i+4)+0] = (double) 0.5*(m->v[3*m->e[8*i + 4] + 0] + cent[0]);
      newv[3*(nv+8*i+4)+1] = (double) 0.5*(m->v[3*m->e[8*i + 4] + 1] + cent[1]);
      newv[3*(nv+8*i+4)+2] = (double) 0.5*(m->v[3*m->e[8*i + 4] + 2] + cent[2]);

      // New vertex f:  (v5 + center) / 2
      newv[3*(nv+8*i+5)+0] = (double) 0.5*(m->v[3*m->e[8*i + 5] + 0] + cent[0]);
      newv[3*(nv+8*i+5)+1] = (double) 0.5*(m->v[3*m->e[8*i + 5] + 1] + cent[1]);
      newv[3*(nv+8*i+5)+2] = (double) 0.5*(m->v[3*m->e[8*i + 5] + 2] + cent[2]);

      // New vertex g:  (v6 + center) / 2
      newv[3*(nv+8*i+6)+0] = (double) 0.5*(m->v[3*m->e[8*i + 6] + 0] + cent[0]);
      newv[3*(nv+8*i+6)+1] = (double) 0.5*(m->v[3*m->e[8*i + 6] + 1] + cent[1]);
      newv[3*(nv+8*i+6)+2] = (double) 0.5*(m->v[3*m->e[8*i + 6] + 2] + cent[2]);

      // New vertex h:  (v7 + center) / 2
      newv[3*(nv+8*i+7)+0] = (double) 0.5*(m->v[3*m->e[8*i + 7] + 0] + cent[0]);
      newv[3*(nv+8*i+7)+1] = (double) 0.5*(m->v[3*m->e[8*i + 7] + 1] + cent[1]);
      newv[3*(nv+8*i+7)+2] = (double) 0.5*(m->v[3*m->e[8*i + 7] + 2] + cent[2]);

      // New element 1: (0, 1, 2, 3, a, b, c, d)
      newe[56*i +  0] = m->e[8*i + 0];
      newe[56*i +  1] = m->e[8*i + 1];
      newe[56*i +  2] = m->e[8*i + 2];
      newe[56*i +  3] = m->e[8*i + 3];
      newe[56*i +  4] = nv + 8*i + 0;
      newe[56*i +  5] = nv + 8*i + 1;
      newe[56*i +  6] = nv + 8*i + 2;
      newe[56*i +  7] = nv + 8*i + 3;

      // New element 2: (4, 5, 0, 1, e, f, a, b)
      newe[56*i +  8] = m->e[8*i + 4];
      newe[56*i +  9] = m->e[8*i + 5];
      newe[56*i + 10] = m->e[8*i + 0];
      newe[56*i + 11] = m->e[8*i + 1];
      newe[56*i + 12] = nv + 8*i + 4;
      newe[56*i + 13] = nv + 8*i + 5;
      newe[56*i + 14] = nv + 8*i + 0;
      newe[56*i + 15] = nv + 8*i + 1;

      // New element 3: (5, 7, 1, 3, f, h, b, d)
      newe[56*i + 16] = m->e[8*i + 5];
      newe[56*i + 17] = m->e[8*i + 7];
      newe[56*i + 18] = m->e[8*i + 1];
      newe[56*i + 19] = m->e[8*i + 3];
      newe[56*i + 20] = nv + 8*i + 5;
      newe[56*i + 21] = nv + 8*i + 7;
      newe[56*i + 22] = nv + 8*i + 1;
      newe[56*i + 23] = nv + 8*i + 3;

      // New element 4: (6, 4, 2, 0, g, e, c, a)
      newe[56*i + 24] = m->e[8*i + 6];
      newe[56*i + 25] = m->e[8*i + 4];
      newe[56*i + 26] = m->e[8*i + 2];
      newe[56*i + 27] = m->e[8*i + 0];
      newe[56*i + 28] = nv + 8*i + 6;
      newe[56*i + 29] = nv + 8*i + 4;
      newe[56*i + 30] = nv + 8*i + 2;
      newe[56*i + 31] = nv + 8*i + 0;

      // New element 5: (7, 6, 3, 2, h, g, d, c)
      newe[56*i + 32] = m->e[8*i + 7];
      newe[56*i + 33] = m->e[8*i + 6];
      newe[56*i + 34] = m->e[8*i + 3];
      newe[56*i + 35] = m->e[8*i + 2];
      newe[56*i + 36] = nv + 8*i + 7;
      newe[56*i + 37] = nv + 8*i + 6;
      newe[56*i + 38] = nv + 8*i + 3;
      newe[56*i + 39] = nv + 8*i + 2;

      // New element 6: (6, 7, 4, 5, g, h, e, f)
      newe[56*i + 40] = m->e[8*i + 6];
      newe[56*i + 41] = m->e[8*i + 7];
      newe[56*i + 42] = m->e[8*i + 4];
      newe[56*i + 43] = m->e[8*i + 5];
      newe[56*i + 44] = nv + 8*i + 6;
      newe[56*i + 45] = nv + 8*i + 7;
      newe[56*i + 46] = nv + 8*i + 4;
      newe[56*i + 47] = nv + 8*i + 5;

      // New element 7: (a, b, c, d, e, f, g, h)
      newe[56*i + 48] = nv + 8*i + 0;
      newe[56*i + 49] = nv + 8*i + 1;
      newe[56*i + 50] = nv + 8*i + 2;
      newe[56*i + 51] = nv + 8*i + 3;
      newe[56*i + 52] = nv + 8*i + 4;
      newe[56*i + 53] = nv + 8*i + 5;
      newe[56*i + 54] = nv + 8*i + 6;
      newe[56*i + 55] = nv + 8*i + 7;
    }

    free(m->v);
    m->v = newv;
    m->nv = nv+8*ne;

    free(m->e);
    m->e = newe;
    m->ne = 7*ne;
  }

#else
  // This version is a "uniform" refinement of the mesh.  Nineteen
  // vertices are added per element, one vertex that bisects each edge (12),
  // one in the center of each face (6), and in in the center of the element.
  // The hexahedron is then split into eight pieces.  If we start with
  // regular hexahedron, then the refined one will have eight regular
  // ones.

  if (!errs) {
    double *newv;
    int *perm;
    int *newe;

    int edges[12][2] = {
      {0, 1}, {0, 2}, {0, 4}, {1, 3}, {1, 5}, {2, 3}, 
      {2, 6}, {3, 7}, {4, 5}, {4, 6}, {5, 7}, {6, 7}
    };

    int faces[6][4] = {
      {0, 1, 2, 3},
      {0, 1, 4, 5},
      {0, 2, 4, 6},
      {1, 3, 5, 7},
      {2, 3, 6, 7},
      {4, 5, 6, 7}
    };

    int vmap[8][8] = {
      {0,  0,  1, 12,  2, 13, 14, 18},
      {1,  3,  0, 12,  4, 15, 13, 18},
      {2,  1,  5, 12,  6, 14, 16, 18},
      {3,  5,  3, 12,  7, 16, 15, 18},
      {4,  8,  2, 13,  9, 17, 14, 18},
      {5, 10,  4, 15,  8, 17, 13, 18},
      {6,  9,  6, 14, 11, 17, 16, 18},
      {7, 11,  7, 16, 10, 17, 15, 18}
    };

    int elem[8];
    int nnv, vert, j;

    // We start off by adding the vertices to the mesh.  One vertex is
    // added per edge, face, and element.  Most edges appear in many 
    // elements.  The upper bound on the number of new vertices 19*ne, 
    // but the final count will be much lower once the duplicates are 
    // removed.  For the vertices, we store four values, three doubles 
    // for the coordinates and an integer for the initial position.  
    // We reallocate the data once we know the correct sizes.

    newv = (double *)malloc(4*sizeof(double)*(nv + 19*ne));
    perm = (int *)malloc(sizeof(int)*(nv + 19*ne));

    // There will be eight elements in the new discretization per element in
    // the old discretization; each element is a set of eight indices.
    // We do not need to reallocate this data.

    newe = (int *)malloc(64*sizeof(int)*ne);

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
      // Start with bisecting each edge.  The sort is done so that the
      // results of the computation are the same for each element where
      // the edge appears.

      for (j = 0; j < 12; ++j) {
        elem[0] = m->e[8*i+edges[j][0]];
        elem[1] = m->e[8*i+edges[j][1]];
        sort2(elem);

	newv[4*(nv+19*i+j+0)+0] = 0.5*(m->v[3*elem[0]+0]+m->v[3*elem[1]+0]);
	newv[4*(nv+19*i+j+0)+1] = 0.5*(m->v[3*elem[0]+1]+m->v[3*elem[1]+1]);
	newv[4*(nv+19*i+j+0)+2] = 0.5*(m->v[3*elem[0]+2]+m->v[3*elem[1]+2]);
        newv[4*(nv+19*i+j+0)+3] = nv + 19*i + j + 0;
        perm[nv + 19*i + j + 0] = -1;
      }

      // Now add a vertex to the center of each face.  The sort is
      // done because each face can appear in two elements.

      for (j = 0; j < 6; ++j) {
        elem[0] = m->e[8*i+faces[j][0]];
        elem[1] = m->e[8*i+faces[j][1]];
        elem[2] = m->e[8*i+faces[j][2]];
        elem[3] = m->e[8*i+faces[j][3]];
        sort4(elem);

	newv[4*(nv+19*i+j+12)+0] = 0.25*(m->v[3*elem[0]+0]+m->v[3*elem[1]+0]+
					 m->v[3*elem[2]+0]+m->v[3*elem[3]+0]);
	newv[4*(nv+19*i+j+12)+1] = 0.25*(m->v[3*elem[0]+1]+m->v[3*elem[1]+1]+
					 m->v[3*elem[2]+1]+m->v[3*elem[3]+1]);
	newv[4*(nv+19*i+j+12)+2] = 0.25*(m->v[3*elem[0]+2]+m->v[3*elem[1]+2]+
					 m->v[3*elem[2]+2]+m->v[3*elem[3]+2]);
        newv[4*(nv+19*i+j+12)+3] = nv + 19*i + j + 12;
        perm[nv + 19*i + j + 12] = -1;
      }

      // Now add a vertex to the center of the element.  The sort is
      // not needed in this case, but done for consistency.

      elem[0] = m->e[8*i+0];
      elem[1] = m->e[8*i+1];
      elem[2] = m->e[8*i+2];
      elem[3] = m->e[8*i+3];
      elem[4] = m->e[8*i+4];
      elem[5] = m->e[8*i+5];
      elem[6] = m->e[8*i+6];
      elem[7] = m->e[8*i+7];
      sort8(elem);

      newv[4*(nv+19*i+18)+0] = 0.125*(m->v[3*elem[0]+0]+m->v[3*elem[1]+0]+
                                      m->v[3*elem[2]+0]+m->v[3*elem[3]+0]+
                                      m->v[3*elem[4]+0]+m->v[3*elem[5]+0]+
                                      m->v[3*elem[6]+0]+m->v[3*elem[7]+0]);
      newv[4*(nv+19*i+18)+1] = 0.125*(m->v[3*elem[0]+1]+m->v[3*elem[1]+1]+
                                      m->v[3*elem[2]+1]+m->v[3*elem[3]+1]+
                                      m->v[3*elem[4]+1]+m->v[3*elem[5]+1]+
                                      m->v[3*elem[6]+1]+m->v[3*elem[7]+1]);
      newv[4*(nv+19*i+18)+2] = 0.125*(m->v[3*elem[0]+2]+m->v[3*elem[1]+2]+
                                      m->v[3*elem[2]+2]+m->v[3*elem[3]+2]+
                                      m->v[3*elem[4]+2]+m->v[3*elem[5]+2]+
                                      m->v[3*elem[6]+2]+m->v[3*elem[7]+2]);
      newv[4*(nv+19*i+18)+3] = nv + 19*i + 18;
      perm[nv + 19*i + 18] = -1;

      // Add the elements now; the first vertex in each element description
      // correspondes to a corner of the original mesh, which the remaining
      // ones are new vertices.

      for (j = 0; j < 8; ++j) {
        newe[8*(8*i+j)+0] = m->e[8*i + vmap[j][0]];
        newe[8*(8*i+j)+1] = nv+19*i+vmap[j][1];
        newe[8*(8*i+j)+2] = nv+19*i+vmap[j][2];
        newe[8*(8*i+j)+3] = nv+19*i+vmap[j][3];
        newe[8*(8*i+j)+4] = nv+19*i+vmap[j][4];
        newe[8*(8*i+j)+5] = nv+19*i+vmap[j][5];
        newe[8*(8*i+j)+6] = nv+19*i+vmap[j][6];
        newe[8*(8*i+j)+7] = nv+19*i+vmap[j][7];
      }
    }

    // Now sort the vertices by their coordinates and index number
    qsort(newv, nv + 19*ne, 4*sizeof(double), vert4);

    // Compute the permutation
    nnv = nv;
    for (i = 0; i < nv + 19*ne; ) {
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
        for (j = i; j < nv + 19*ne; ++j) {
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
    for (; i < nv + 19*ne; ++i) {
      assert(perm[i] > 0 && perm[i] < nnv);
    }

    // Now permute the vertices
    free(m->v);
    m->v = (double *)malloc(3*sizeof(double)*nnv);
    m->nv = nnv;

    for (i = 0; i < nv + 19*ne; ++i) {
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
      m->e[8*i + 0] = perm[m->e[8*i + 0]];
      m->e[8*i + 1] = perm[m->e[8*i + 1]];
      m->e[8*i + 2] = perm[m->e[8*i + 2]];
      m->e[8*i + 3] = perm[m->e[8*i + 3]];
      m->e[8*i + 4] = perm[m->e[8*i + 4]];
      m->e[8*i + 5] = perm[m->e[8*i + 5]];
      m->e[8*i + 6] = perm[m->e[8*i + 6]];
      m->e[8*i + 7] = perm[m->e[8*i + 7]];
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
    for (j = 0; j < 8; ++j) {
      if ((e[j] < 0) || (e[j] >= nv)) {
        printf("Index error for element %d\n", i);
	errs = 1;
	break;
      }
    }
    e += 8;
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
  qsort(vd, nv, 4*sizeof(double), vert3);

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
      for (k = 0; k < 8; ++k) {
	if (ed[k] != perm[ed[k]]) {
	  ed[k] = perm[ed[k]];
	  printf("Fixing %d.%d\n", j, k);
	}
      }
      ed += 8;
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
    for (j = 0; j < 8; ++j) {
      ++b[e[j]];
    }
    e += 8;
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
      for (j = 0; j < 8; ++j) {
	e[j] = b[e[j]];
      }
      e += 8;
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

  f = m->f;
  b = m->b;
  for (i = 0; i < nf; ++i) {
    b[f[0]] = 2;
    f += 1;
  }

  /* Finish with boundary vertices; compute and sort the faces */
  ls = (int *)malloc(   sizeof(int)*nv);         /* vertex starts */
  l1 = (int *)malloc(36*sizeof(int)*(ne+1));     /* list of faces */
  l2 = (int *)malloc(36*sizeof(int)*(ne+1));     /* list of faces */

  e = m->e;
  for (i = 0; i < ne; ++i) {
    /* Hypercube ordering for the vertices */

    l1[0] = e[0]; l1[1] = e[1]; l1[2] = e[2]; l1[3] = e[3];
    l1[4] = i; l1[5] = 0;
    sort4(l1); l1 += 6;

    l1[0] = e[0]; l1[1] = e[1]; l1[2] = e[4]; l1[3] = e[5];
    l1[4] = i; l1[5] = 1;
    sort4(l1); l1 += 6;

    l1[0] = e[0]; l1[1] = e[2]; l1[2] = e[4]; l1[3] = e[6];
    l1[4] = i; l1[5] = 2;
    sort4(l1); l1 += 6;

    l1[0] = e[1]; l1[1] = e[3]; l1[2] = e[5]; l1[3] = e[7];
    l1[4] = i; l1[5] = 3;
    sort4(l1); l1 += 6;

    l1[0] = e[2]; l1[1] = e[3]; l1[2] = e[6]; l1[3] = e[7];
    l1[4] = i; l1[5] = 4;
    sort4(l1); l1 += 6;

    l1[0] = e[4]; l1[1] = e[5]; l1[2] = e[6]; l1[3] = e[7];
    l1[4] = i; l1[5] = 5;
    sort4(l1); l1 += 6;

    e += 8;
  }

  l1[0] = -1;
  l1[1] = -1;
  l1[2] = -1;
  l1[3] = -1;
  l1[4] = -1;
  l1[5] = -1;

  l1 -= 36*ne;

  /* Now perform the radix sort */
  radix4(l2, l1, ls, 3, nv, 6*ne);
  radix4(l1, l2, ls, 2, nv, 6*ne);
  radix4(l2, l1, ls, 1, nv, 6*ne);
  radix4(l1, l2, ls, 0, nv, 6*ne);

  free(ls);
  free(l2);

  for (i = 0; i < 6*ne; ++i) {
    if ((l1[0] == l1[6]) && (l1[1] == l1[7]) &&
        (l1[2] == l1[8]) && (l1[3] == l1[9])) {
      /* Face must be in only two elements! */
      cnt = 1;

#ifdef CHECK_MESH
      /* Check compatibility conditions for the faces             */
      /* Make sure one is clockwise and other is counter clockwise*/

      e = m->e + 8*l1[4];
      switch(l1[5]) {
      case 0:
        l1[0] = e[0]; l1[1] = e[2]; l1[2] = e[3]; l1[3] = e[1];
        break;

      case 1:
        l1[0] = e[0]; l1[1] = e[1]; l1[2] = e[5]; l1[3] = e[4];
        break;

      case 2:
        l1[0] = e[0]; l1[1] = e[4]; l1[2] = e[6]; l1[3] = e[2];
        break;

      case 3:
        l1[0] = e[1]; l1[1] = e[3]; l1[2] = e[7]; l1[3] = e[5];
        break;

      case 4:
        l1[0] = e[2]; l1[1] = e[6]; l1[2] = e[7]; l1[3] = e[3];
        break;

      case 5:
        l1[0] = e[4]; l1[1] = e[5]; l1[2] = e[7]; l1[3] = e[6];
        break;
      }
      face4(l1);

      e = m->e + 8*l1[10];
      switch(l1[11]) {
      case 0:
        l1[6] = e[0]; l1[7] = e[1]; l1[8] = e[3]; l1[9] = e[2];
        break;

      case 1:
        l1[6] = e[0]; l1[7] = e[4]; l1[8] = e[5]; l1[9] = e[1];
        break;

      case 2:
        l1[6] = e[0]; l1[7] = e[2]; l1[8] = e[6]; l1[9] = e[4];
        break;

      case 3:
        l1[6] = e[1]; l1[7] = e[5]; l1[8] = e[7]; l1[9] = e[3];
        break;

      case 4:
        l1[6] = e[2]; l1[7] = e[3]; l1[8] = e[7]; l1[9] = e[6];
        break;

      case 5:
        l1[6] = e[4]; l1[7] = e[6]; l1[8] = e[7]; l1[9] = e[5];
        break;
      }
      face4(l1 + 6);

      if ((l1[0] != l1[6]) || (l1[1] != l1[7])  ||
          (l1[2] != l1[8]) || (l1[3] != l1[9])) {
        printf("Face Compatiblity: %d.%d and %d.%d\n",
               l1[4], l1[5], l1[10], l1[11]);
	++errs;
      }

      sort4(l1);
      sort4(l1 + 6);
#endif

      /* Remove all replicated vertices */
      while ((i < 6*ne) &&
             (l1[0] == l1[6]) && (l1[1] == l1[7]) &&
             (l1[2] == l1[8]) && (l1[3] == l1[9])) {
        l1 += 6;
        ++i;
	++cnt;
      }

      if (cnt > 2) {
	printf("Face Count Problem: %d.%d = %d\n", l1[4], l1[5], cnt);
	++errs;
      }
    }
    else {
      /* This face is on the boundary */
      if (0 == b[l1[0]]) {
        b[l1[0]] = 1;
      }
      if (0 == b[l1[1]]) {
        b[l1[1]] = 1;
      }
      if (0 == b[l1[2]]) {
        b[l1[2]] = 1;
      }
      if (0 == b[l1[3]]) {
        b[l1[3]] = 1;
      }
      ++nse;
    }
    l1 += 6;
  }
  l1 -= 36*ne;

  m->nse = nse;

#if defined(CHECK_MESH) || defined(SEGMENT_MESH)
  /* Create surface mesh from the exposed faces */
  m->se = (int *)malloc(4*sizeof(int)*nse);
  se = m->se;

  for (i = 0; i < 6*ne; ++i) {
    if ((l1[0] == l1[6]) && (l1[1] == l1[7]) &&
        (l1[2] == l1[8]) && (l1[3] == l1[9])) {
      /* Face in two elements; remove replicated faces */
      while ((i < 6*ne) &&
             (l1[0] == l1[6]) && (l1[1] == l1[7]) &&
             (l1[2] == l1[8]) && (l1[3] == l1[9])) {
        l1 += 6;
        ++i;
      }
    }
    else {
      /* This face is on the boundary */
      e = m->e + 8*l1[4];
      switch(l1[5]) {
      case 0:
        se[0] = e[0]; se[1] = e[2]; se[2] = e[1]; se[3] = e[3];
        break;

      case 1:
        se[0] = e[0]; se[1] = e[1]; se[2] = e[4]; se[3] = e[5];
        break;

      case 2:
        se[0] = e[0]; se[1] = e[4]; se[2] = e[2]; se[3] = e[6];
        break;

      case 3:
        se[0] = e[1]; se[1] = e[3]; se[2] = e[5]; se[3] = e[7];
        break;

      case 4:
        se[0] = e[2]; se[1] = e[6]; se[2] = e[3]; se[3] = e[7];
        break;

      case 5:
        se[0] = e[4]; se[1] = e[5]; se[2] = e[6]; se[3] = e[7];
        break;
      }
      se += 4;
    }
    l1 += 6;
  }
  l1 -= 36*ne;
#endif

  free(l1);
  return errs;
}

#if defined(CHECK_MESH) || defined(SEGMENT_MESH)
static int checkSurfData(MeshData *m)
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
  si = (int *)calloc(4*nse, sizeof(int));
  m->si = si;
#endif

  /* Sort edges in the surface mesh */
  ls = (int *)malloc(   sizeof(int)*nv);          /* vertex starts */
  l1 = (int *)malloc(16*sizeof(int)*(nse+1));     /* list of faces */
  l2 = (int *)malloc(16*sizeof(int)*(nse+1));     /* list of faces */

  for (i = 0; i < nse; ++i) {
    /* Hypercube ordering for the vertices */

    l1[0] = se[0]; l1[1] = se[1];
    l1[2] = i; l1[3] = 0;
    sort2(l1); l1 += 4;

    l1[0] = se[0]; l1[1] = se[2];
    l1[2] = i; l1[3] = 1;
    sort2(l1); l1 += 4;

    l1[0] = se[1]; l1[1] = se[3];
    l1[2] = i; l1[3] = 2;
    sort2(l1); l1 += 4;

    l1[0] = se[2]; l1[1] = se[3];
    l1[2] = i; l1[3] = 3;
    sort2(l1); l1 += 4;

    se += 4;
  }

  l1[0] = -1;
  l1[1] = -1;
  l1[2] = -1;
  l1[3] = -1;

  l1 -= 16*nse;

  /* Now perform the radix sort */
  radix2(l2, l1, ls, 1, nv, 4*nse);
  radix2(l1, l2, ls, 0, nv, 4*nse);

  free(ls);
  free(l2);

  for (i = 0; i < 4*nse; ++i) {
    if ((l1[0] == l1[4]) && (l1[1] == l1[5])) {
      /* Edge must be in only two elements! */
      cnt = 1;

      /* Check compatibility conditions for the edges             */
      /* Make sure one is clockwise and other is counter clockwise*/
      se = m->se + 4*l1[2];

      switch(l1[3]) {
      case 0:
        l1[0] = se[0]; l1[1] = se[1];
        break;

      case 1:
        l1[0] = se[2]; l1[1] = se[0];
        break;

      case 2:
        l1[0] = se[1]; l1[1] = se[3];
        break;

      case 3:
        l1[0] = se[3]; l1[1] = se[2];
        break;
      }

      se = m->se + 4*l1[6];
      switch(l1[7]) {
      case 0:
        l1[4] = se[1]; l1[5] = se[0];
        break;

      case 1:
        l1[4] = se[0]; l1[5] = se[2];
        break;

      case 2:
        l1[4] = se[3]; l1[5] = se[1];
        break;

      case 3:
        l1[4] = se[2]; l1[5] = se[3];
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
      si[4*l1[2] + l1[3]] = l1[6];
      si[4*l1[6] + l1[7]] = l1[2];
#endif

      /* Remove all replicated vertices */
      while ((i < 4*nse) && (l1[0] == l1[4]) && (l1[1] == l1[5])) {
        l1 += 4;
        ++i;
	++cnt;
      }

      if (cnt > 2) {
	printf("Edge Count Problem: %d.%d = %d\n", l1[2], l1[3], cnt);
	/* printf("  Vertices %d, %d\n", l1[0], l1[1]); */
	++errs;
      }
    }
    else {
      printf("Closed Surface Problem: %d.%d\n", l1[2], l1[3]);
      ++errs;
    }
    l1 += 4;
  }

  l1 -= 16*nse;
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

  int *ss;
  int *q;

  int nss = 0;

  int qs, qe;
  int i, idx, loc;

  if (nv <= 0) {
    return 0;
  }

  ss = (int *)calloc(nse, sizeof(int));
  m->ss = ss;

  q = (int *)malloc(nse*sizeof(int));
  
  idx = 0;
  while (idx >= 0) {
    ++nss;
    ss[idx] = nss;

    qs = 0;
    qe = 0;

    loc = 4*idx;
    for (i = 0; i < 4; ++i) {
      if (!ss[si[loc+i]]) {
	/* Mark as being in the segment and add index to queue */
	ss[si[loc+i]] = nss;
	q[qe] = si[loc+i];
	qe = (qe+1) % nse;
      }
    }

    while (qs != qe) {
      loc = 4*q[qs];
      qs = (qs+1) % nse;

      for (i = 0; i < 4; ++i) {
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
  printf("Surface Segments: %d\n", nss);
  free(q);

#if 1
  printf("after segments b[5064] = %d\n", m->b[5064]);

  for (i = 0; i < m->nv; ++i) {
    m->b[i] = -m->b[i];
  }

  printf("after negation b[5064] = %d\n", m->b[5064]);

  for (i = 0; i < nse; ++i) {
    m->b[m->se[4*i+0]] = ss[i];
    m->b[m->se[4*i+1]] = ss[i];
    m->b[m->se[4*i+2]] = ss[i];
    m->b[m->se[4*i+3]] = ss[i];
  }

  printf("after segment b[5064] = %d\n", m->b[5064]);

  for (i = 0; i < m->nv; ++i) {
    printf("%4d\n", m->b[i]);
  }
#endif

#if 0
  {
    FILE *fp;
    char fname[256];
    int j;

    for (i = 1; i <= nss; ++i) {
      sprintf(fname, "s%d.dat", i);
      fp = fopen(fname, "w");

      for (j = 0; j < nse; ++j) {
        if (i == ss[j]) {
          fprintf(fp, "%10.9e %10.9e %10.9e %10.9e %10.9e\n",
		  m->v[3*m->se[4*j+0]], m->v[3*m->se[4*j+1]],
		  m->v[3*m->se[4*j+3]], m->v[3*m->se[4*j+2]],
		  m->v[3*m->se[4*j+0]]);
          fprintf(fp, "%10.9e %10.9e %10.9e %10.9e %10.9e\n",
		  m->v[3*m->se[4*j+0]+1], m->v[3*m->se[4*j+1]+1],
		  m->v[3*m->se[4*j+3]+1], m->v[3*m->se[4*j+2]+1],
		  m->v[3*m->se[4*j+0]+1]);
          fprintf(fp, "%10.9e %10.9e %10.9e %10.9e %10.9e\n",
		  m->v[3*m->se[4*j+0]+2], m->v[3*m->se[4*j+1]+2],
		  m->v[3*m->se[4*j+3]+2], m->v[3*m->se[4*j+2]+2],
		  m->v[3*m->se[4*j+0]+2]);
        }
      }
      fclose(fp);
    }
  }
#endif

  return 0;
}
#endif

static int computeMesh(Mesh *m) 
{
  const MeshData *md = m->d;

  const int hn = md->nv;
  const int he = md->ne;

  int *h = md->e;
  int *t = (int *)malloc(32*sizeof(int)*he);

  int i;

  m->nv = hn;
  m->ne = 8*he;

  m->e = t;
  for (i = 0; i < he; ++i) {
    t[0] = h[0]; t[1] = h[1]; t[2] = h[2]; t[3] = h[4]; t += 4;
    t[0] = h[1]; t[1] = h[3]; t[2] = h[0]; t[3] = h[5]; t += 4;
    t[0] = h[2]; t[1] = h[0]; t[2] = h[3]; t[3] = h[6]; t += 4;
    t[0] = h[3]; t[1] = h[2]; t[2] = h[1]; t[3] = h[7]; t += 4;
    t[0] = h[4]; t[1] = h[6]; t[2] = h[5]; t[3] = h[0]; t += 4;
    t[0] = h[5]; t[1] = h[4]; t[2] = h[7]; t[3] = h[1]; t += 4;
    t[0] = h[6]; t[1] = h[7]; t[2] = h[4]; t[3] = h[2]; t += 4;
    t[0] = h[7]; t[1] = h[5]; t[2] = h[6]; t[3] = h[3]; t += 4;
    h += 8;
  }

  m->v = (double *)malloc(3*sizeof(double)*hn);
  m->p = (int    *)malloc(1*sizeof(int   )*hn);

  memcpy(m->v, md->v, 3*sizeof(double)*hn);
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
  int     elem[8];

  int     v1, v2, v3, v4;
  int     i, j, pos, neg, deg;
  int     errs = 0;
  
#ifdef CHECK_PLANAR
  int     plane[24] = {0, 1, 3, 4,
		       0, 1, 5, 4,
                       1, 3, 7, 5,
                       3, 2, 6, 7,
                       2, 0, 4, 6,
		       4, 5, 7, 6};
#endif

#ifdef USE_WEIGHT
  double weight[6] = {1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
  double target[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
#endif

  for (i = 0; i < ne; ++i) {
    pos = 0;
    neg = 0;
    deg = 0;

    for (j = 0; j < 8; ++j) {
      elem[j] = e[0];

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
	  printf("Neg: %5.4e\n", f);
	  ++neg;
	}
	else {
	  ++deg;
	}
      } 
      else {
	++pos;
      }
#else
      if (o_fcn(&f, x, weight, target)) {
	if (f < -epsilon) {
	  printf("Neg: %5.4e\n", f);
	  ++neg;
	}
	else {
	  ++deg;
	}
      } 
      else {
	++pos;
      }
#endif
    }

    if (8 != pos) {
      if (8 == deg) {
	printf("Invalid element: %d: zero area.\n", i);
      }
      else if (8 == neg) {
	printf("Invalid element: %d: inverted.\n", i);
      }
      else {
	printf("Invalid element: %d: twisted: "
               "%d positive, %d negative, %d degenerate.\n", 
	       i, pos, neg, deg);
      }
      ++errs;
    }

#ifdef CHECK_PLANAR
    pos = 0;
    neg = 0;
    deg = 0;

    for (j = 0; j < 6; ++j) {
      v1 = elem[plane[4*j + 0]];
      v2 = elem[plane[4*j + 1]];
      v3 = elem[plane[4*j + 2]];
      v4 = elem[plane[4*j + 3]];

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
	  printf("Neg face: %5.4e\n", f);
	  ++neg;
	}
	else {
	  ++deg;
	}
      } 
      else {
	printf("Pos face: %5.4e\n", f);
	++pos;
      }
#else
      if (o_fcn(&f, x, weight, target)) {
	if (f < -epsilon) {
	  printf("Neg face: %5.4e\n", f);
	  ++neg;
	}
	else {
	  ++deg;
	}
      } 
      else {
	printf("Pos face: %5.4e\n", f);
	++pos;
      }
#endif
    }

    if (6 != deg) {
      printf("Invalid element: %d: faces not planar: %d\n", i, neg + pos);
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
  if (checkSurfData(md)) {
    printf("Surface check failure.\n");
    freeMeshData(&md);
    return -2;
  }
#endif

#ifdef SEGMENT_MESH
  if (segmentSurfData(md)) {
    printf("Segmentation segmentation failure.\n");
    freeMeshData(&md);
    return -2;
  }
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

#ifdef SWEEP
  e = m->e;
  for (i = 0; i < ne; ++i) {
    fprintf(fp, "%d %d %d %d %d %d %d %d\n",
            e[0], e[1], e[3], e[2], e[4], e[5], e[7], e[6]);
    e += 8;
  }
#else
  e = m->e;
  for (i = 0; i < ne; ++i) {
    fprintf(fp, "%d %d %d %d %d %d %d %d\n",
            e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7]);
    e += 8;
  }
#endif

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
  fprintf(fp, "param E := %d;\n", 8*ne);

  fprintf(fp, "\nvar x : 1 2 3 =\n");
  v = m->v;
  for (i = 0; i < nv; ++i) {
    fprintf(fp, "%6d % 15.14e % 15.14e % 15.14e\n", i+1, v[0], v[1], v[2]);
    v += 3;
  }

  fprintf(fp, ";\n\nparam TETS : 1 2 3 4=\n");
  e = m->e;
  for (i = 0; i < ne; ++i) {
    fprintf(fp, "%6d %6d %6d %6d %6d\n", 8*i+1,e[0]+1,e[1]+1,e[2]+1,e[4]+1);
    fprintf(fp, "%6d %6d %6d %6d %6d\n", 8*i+2,e[1]+1,e[3]+1,e[0]+1,e[5]+1);
    fprintf(fp, "%6d %6d %6d %6d %6d\n", 8*i+3,e[2]+1,e[0]+1,e[3]+1,e[6]+1);
    fprintf(fp, "%6d %6d %6d %6d %6d\n", 8*i+4,e[3]+1,e[2]+1,e[1]+1,e[7]+1);
    fprintf(fp, "%6d %6d %6d %6d %6d\n", 8*i+5,e[4]+1,e[6]+1,e[5]+1,e[0]+1);
    fprintf(fp, "%6d %6d %6d %6d %6d\n", 8*i+6,e[5]+1,e[4]+1,e[7]+1,e[1]+1);
    fprintf(fp, "%6d %6d %6d %6d %6d\n", 8*i+7,e[6]+1,e[7]+1,e[4]+1,e[2]+1);
    fprintf(fp, "%6d %6d %6d %6d %6d\n", 8*i+8,e[7]+1,e[5]+1,e[6]+1,e[3]+1);
    e += 8;
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
  md->e = (int *)malloc(8*sizeof(int)*ne);
  memcpy(md->e, elems, 8*sizeof(int)*ne);

#ifdef SWEEP
 {
   int *e = md->e;
   int i, tmpv;

   for (i = 0; i < ne; ++i) {
     tmpv = e[2];
     e[2] = e[3];
     e[3] = tmpv;

     tmpv = e[6];
     e[6] = e[7];
     e[7] = tmpv;
     e += 8;
   }
 }
#endif

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


#include <sys/time.h>
#include <sys/resource.h>

#include "opt.h"

#ifndef MICROSEC
#  define MICROSEC 1000000
#endif

#if defined(LIBRARY)
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

#ifdef SWEEP
  int tmp;
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

#ifdef SWEEP
  for (i = 0; i < ne; ++i) {
    tmp = e[2];
    e[2] = e[3];
    e[3] = tmp;

    tmp = e[6];
    e[6] = e[7];
    e[7] = tmp;

    e += 8;
  }
  e = md->e;
#endif

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
  if (checkSurfData(md)) {
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

#ifdef SWEEP
  for (i = 0; i < ne; ++i) {
    tmp = e[2];
    e[2] = e[3];
    e[3] = tmp;

    tmp = e[6];
    e[6] = e[7];
    e[7] = tmp;

    e += 8;
  }
  e = md->e;
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
