#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mesh.h"

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

  (*m)->edat = NULL;
  (*m)->ecat = NULL;
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

  (*m)->edat = NULL;
  (*m)->ecat = NULL;

  (*m)->p = NULL;
  (*m)->i = NULL;
  (*m)->g = NULL;
  
  (*m)->len = NULL;
  (*m)->col = NULL;
  (*m)->dat = NULL;
  (*m)->inst = NULL;
  
  (*m)->per = NULL;
  (*m)->pel = NULL;
  (*m)->iper = NULL;
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

  if (NULL != (*m)->edat) { free((*m)->edat); }
  if (NULL != (*m)->ecat) { free((*m)->ecat); }

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
  if (NULL != (*m)->edat) { free((*m)->edat); }
  if (NULL != (*m)->ecat) { free((*m)->ecat); }

  if (NULL != (*m)->p) { free((*m)->p); }
  if (NULL != (*m)->i) { free((*m)->i); }
  if (NULL != (*m)->g) { free((*m)->g); }
  if (NULL != (*m)->len) { free((*m)->len); }
  if (NULL != (*m)->col) { free((*m)->col); }
  if (NULL != (*m)->dat) { free((*m)->dat); }
  if (NULL != (*m)->inst) { free((*m)->inst); }
  if (NULL != (*m)->per) { free((*m)->per); }
  if (NULL != (*m)->pel) { free((*m)->pel); }
  if (NULL != (*m)->iper) { free((*m)->iper); }

  free(*m);
  (*m) = NULL;
  return 0;
}

