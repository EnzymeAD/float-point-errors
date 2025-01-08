#ifndef PRE_H
#define PRE_H

#include "mesh.h"

typedef struct _Precond {
  void (*destroy)(struct _Precond *);
  void (*calc   )(struct _Precond *, const Mesh *);
  void (*apply  )(double *, double *, const struct _Precond *, const Mesh *);
  
  double *data;

  int     max_nn;
  int     max_nz;
} Precond;

Precond *preCreate(const int p, const int max_nn, const int max_nz);

#endif
