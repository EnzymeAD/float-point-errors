#ifndef MESH_H
#define MESH_H

/*****************************************************************************/
/* Hexahedral mesh in hypercube format                                       */
/*****************************************************************************/

typedef struct _MeshData 
{
  double  *v;   /* Vertices of the mesh                                      */
  int     *e;   /* Elements of the mesh (quad or hex)                        */
  int     *f;   /* Fixed vertices in the mesh                                */
  int     *b;   /* Boundary vetrices computed for the mesh                   */

  int    *se;   /* Surface elements in the mesh				     */
  int    *si;	/* Surface incidence list				     */
  int    *ss;	/* Surface segment indicator 				     */

  int    *pe;   /* Planar elements in mesh                                   */
  int    *ps;   /* Start of plane segment i in pe                            */
  double *pd;   /* Data for plane segment i                                  */

  int    *he;   /* Spherical elements in mesh                                */
  int    *hs;   /* Start of spherical segment i in he                        */
  double *hd;   /* Data for spherical segment i                              */

  int *part;    /* Partition number (not used)				     */

  int nv;	/* Number of nodes in mesh                                   */
  int ne;	/* Number of elements in mesh                                */
  int nf;	/* Number of fixed nodes in mesh                             */
  int np;	/* Number of partitions	(not used)			     */
  int no;	/* Amount of overlap (not used) 			     */

  int nse;	/* Number of surface elements in mesh                        */
  int nss;	/* Number of surface segments in mesh                        */

  int npe;      /* Number of planar elements in mesh                         */
  int nps;      /* Number of planar segments in mesh                         */

  int nhe;      /* Number of spherical elements in mesh                      */
  int nhs;      /* Number of spherical segments in mesh                      */
} MeshData;

/*****************************************************************************/
/* Mesh structure.  Defines a tri or tet mesh depending upon the code        */
/* linked in.                                                            */
/*****************************************************************************/

typedef struct _Mesh 
{
  MeshData *d;  /* Data for the mesh (quad and hex)                          */

  double *v;	/* Vertices of the mesh                                      */
  int    *e;    /* Elements of the mesh (tri or tet)                         */

  double *w;	/* Weight matrices for each element			     */
		/* 16 values per element				     */

  int    *p;	/* Permutation vector -- offset in v                         */
                /* -- nonnegative value gives location in compacted vector   */
		/* -- negative value means the coordinate is fixed           */
  int    *i;    /* Inverse permutation                                       */

  double *g;	/* Gradient vector                                           */

  int    *len;	/* Hessian -- length of the row (in blocks)                  */
  int    *col;	/* Hessian -- column number of the row (in blocks)           */
  double *dat;	/* Hessian -- data of the row                                */
  int    *inst;	/* Hessian -- accumulation instructions                      */

  int    *per;	/* Reordering -- vertex permutation                          */

  int nv;	/* Number of nodes in mesh                                   */
  int ne;	/* Number of elements in mesh                                */

  int nf;	/* Number of fixed nodes in mesh                             */
  int nn;	/* Number of non-fixed nodes in mesh                         */

  int nb;	/* Number of blocks in hessian calculation                   */
  int ndb;	/* Number of diagonal blocks in hessian calculation          */
  int nob;	/* Number of off-diagonal blocks in hessian calculation      */

  int nz;	/* Total number of nonzeros in Hessian                       */
} Mesh;

int allocMesh(Mesh **m);
int allocMeshData(MeshData **m);

int freeMeshData(MeshData **m);

int readMesh(const char *fname, Mesh **m);
int writeMesh(const char *fname, Mesh *m);
int writeMesh_AMPL(const char *fname, const char *fname_opt, Mesh *m);

/* Sorts used for checking mesh. */
int vert2(const void *p1, const void *p2);
int vert3(const void *p1, const void *p2);
int vert4(const void *p1, const void *p2);

void face3(int *l);
void face4(int *l);

void sort2(int *l);
void sort3(int *l);
void sort4(int *l);
void sort8(int *l);

void radix2(int *dest, int *src, int *len,
	    const int idx, const int nv, const int ne);
void radix3(int *dest, int *src, int *len,
	    const int idx, const int nv, const int ne);
void radix4(int *dest, int *src, int *len,
	    const int idx, const int nv, const int ne);

/* Used for reading and writing. */
int reorderMesh(Mesh *m);
int finishMesh(Mesh *m);

/* API calls */
int createMesh(double *verts, const int nv, 
	       int    *elems, const int ne, 
	       int    *fixed, const int nf, 
	       Mesh **m);
int freeMesh(Mesh **m);

int setVertices(double *verts, int nv, Mesh *m);
int getVertices(double *verts, int nv, Mesh *m);

#endif
