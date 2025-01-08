#ifndef MESH_H
#define MESH_H

/*****************************************************************************/
/* Hexahedral mesh in hypercube format                                       */
/*****************************************************************************/

typedef struct _MeshData 
{
  double *v;    /* Vertices of the mesh                                      */
  int    *e;    /* Elements of the mesh (quad or hex)                        */
  int    *f;    /* Fixed vertices in the mesh                                */
  int    *b;    /* Boundary vetrices computed for the mesh                   */
  int    *part;	/* Partition number for the vertex                           */

  double *edat; /* Amount to squash for ideal element                        */
  int    *ecat; /* Ideal element category:                                   */
                /*   1 -- squash only in z direction                         */
	        /*   2 -- squash only in y, z directions                     */
	        /*   3 -- squash in all directions                           */

  int nv;	/* Number of nodes in mesh                                   */
  int ne;	/* Number of elements in mesh                                */
  int nf;	/* Number of fixed nodes in mesh                             */
  int np;	/* Number of partitions					     */
  int no;	/* Amount of overlap                                         */
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

  double *edat; /* Amount to squash for ideal element                        */
  int    *ecat; /* Ideal element category:                                   */
                /*   1 -- squash only in z direction                         */
	        /*   2 -- squash only in y, z directions                     */
	        /*   3 -- squash in all directions                           */

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
  int    *pel;	/* Reordering -- element permutation                         */
  int    *iper;	/* Reordering -- inverse vertex permutation                  */

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

int freeMesh(Mesh **m);
int freeMeshData(MeshData **m);

int readMesh(const char *fname, Mesh **m);
int writeMesh(const char *fname, Mesh *m);
int writeMesh_AMPL(const char *fname, const char *fname_opt, Mesh *m);

/* In aux -- used for reading and writing. */
int reorderMesh(Mesh *m);
int finishMesh(Mesh *m);

#endif
