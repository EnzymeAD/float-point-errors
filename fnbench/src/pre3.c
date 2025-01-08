#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mesh.h"
#include "pre.h"

/*****************************************************************************/
/* Preconditioners:                                                          */
/*   1 -- diagonal                                                           */
/*   2 -- block cholesky                                                     */
/*****************************************************************************/

static void preDestroy0(Precond *p)
{
  free(p);
  return;
}

static void preCalc0(Precond *p, const Mesh *mesh)
{
  return;
}

static void preApply0(double *z, double *v, const Precond *p, const Mesh *mesh)
{
  memcpy(z, v, mesh->nn*sizeof(double));
  return;
}

static Precond *preCreate0(const int max_nn, const int max_nz)
{
  Precond *p;

  p = (Precond *)malloc(sizeof(Precond));
  p->destroy = preDestroy0;
  p->calc    = preCalc0;
  p->apply   = preApply0;
  p->data    = NULL;

  p->max_nn = max_nn;
  p->max_nz = max_nz;
  return p;
}

static void preDestroy1(Precond *p)
{
  free(p->data);
  free(p);
  return;
}

static void preCalc1(Precond *p, const Mesh *mesh)
{
  int    *l = mesh->len;
  double *d = mesh->dat;
  double *pd = p->data;

  const int nn = mesh->nn;
  int i;

  for (i = 0; i < nn; ++i) {
    pd[0] = 1.0 / (d[0] + d[3] + d[5]);

    pd += 1;
    d  += 6 + 9*(*l++ - 1);
  }
  return; 
}

static void preApply1(double *z, double *v, const Precond *p, const Mesh *mesh)
{
  double *pd = p->data;

  const int nn = mesh->nn;
  int i;

  for (i = 0; i < nn; ++i) {
    z[0] = v[0]*pd[0];
    z[1] = v[1]*pd[0];
    z[2] = v[2]*pd[0];

    pd += 1;
    v  += 3;
    z  += 3;
  }
  return;
}

static Precond *preCreate1(const int max_nn, const int max_nz)
{
  Precond *p;

  p = (Precond *)malloc(sizeof(Precond));
  p->destroy = preDestroy1;
  p->calc    = preCalc1;
  p->apply   = preApply1;
  p->data    = (double *)malloc(max_nn*sizeof(double));

  p->max_nn = max_nn;
  p->max_nz = max_nz;
  return p;
}

static void preDestroy2(Precond *p)
{
  free(p->data);
  free(p);
  return;
}

static void preCalc2(Precond *p, const Mesh *mesh)
{
  int    *l = mesh->len;
  double *d = mesh->dat;
  double *pd = p->data;

  const int nn = mesh->nn;
  int i;

  for (i = 0; i < nn; ++i) {
    pd[0] = 1.0  / d[0];
    pd[1] = d[1] * pd[0];
    pd[2] = d[2] * pd[0];

    pd[3] = 1.0  / (d[3]-d[1]*pd[1]);

    pd[5] = d[4] - d[2]*pd[1];
    pd[4] = pd[3] * pd[5];

    pd[5] = 1.0  / (d[5]-d[2]*pd[2]-pd[4]*pd[5]);

    pd += 6;
    d  += 6 + 9*(*l++ - 1);
  }
  return; 
}

static void preApply2(double *z, double *v, const Precond *p, const Mesh *mesh)
{
  double *pd = p->data;

  const int nn = mesh->nn;
  int i;

  for (i = 0; i < nn; ++i) {
    z[0] = v[0]; 
    z[1] = v[1] - pd[1]*z[0];
    z[2] = v[2] - pd[2]*z[0] - pd[4]*z[1];
 
    z[0] *= pd[0]; 
    z[1] *= pd[3];
    z[2] *= pd[5];
 
    z[1] -= pd[4]*z[2];
    z[0] -= pd[1]*z[1] + pd[2]*z[2];

    pd += 6;
    v  += 3;
    z  += 3;
  }
  return;
}

static Precond *preCreate2(const int max_nn, const int max_nz)
{
  Precond *p;

  p = (Precond *)malloc(sizeof(Precond));
  p->destroy = preDestroy2;
  p->calc    = preCalc2;
  p->apply   = preApply2;
  p->data    = (double *)malloc(6*max_nn*sizeof(double));

  p->max_nn = max_nn;
  p->max_nz = max_nz;
  return p;
}

Precond *preCreate(const int p, const int max_nn, const int max_nz)
{
  switch(p) {
  case 0:
    return preCreate0(max_nn, max_nz);

  case 1:
    return preCreate1(max_nn, max_nz);

  case 2:
    return preCreate2(max_nn, max_nz);
  }
  return NULL;
}

