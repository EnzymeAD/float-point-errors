#ifndef FEASNEWT_H
#define FEASNEWT_H

#include "mesh.h"

int optMesh(Mesh *m, int max_iter, double conv_tol, int precond);

#endif

