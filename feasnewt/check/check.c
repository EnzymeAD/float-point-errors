#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "check.h"
#include "fcn.h"

int check_mesh(Mesh *mesh)
{
  int    nn = mesh->nn;

  if (nn <= 0) {
    /* No nodes!  Just return.                                               */
    return 0;
  }

  if (cFcn(mesh)) {
    fprintf(stderr, "Invalid starting point.\n");
    exit(-1);
  }
  return 0;
}

