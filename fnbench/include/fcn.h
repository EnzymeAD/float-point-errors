#ifndef FCN_H
#define FCN_H

#include "mesh.h"

/*****************************************************************************/
/* Tolerance used in the code for small/negative area.                       */
/*****************************************************************************/

#define epsilon 1.0e-14

/*****************************************************************************/
/* Functions to initialize the mesh so that the gradient and Hessian can     */
/* be correctly accumulated.                                                 */
/*****************************************************************************/

void gMesh(Mesh *m);
void hMesh(Mesh *m);

void hMesh_d(Mesh *m);

double gNorm(const Mesh *m);

/*****************************************************************************/
/* Functions to obtain the global objective, gradient, and Hessian.          */
/*****************************************************************************/

int cFcn(const Mesh *m);

int oFcn(double *obj, const Mesh *m);
int gFcn(double *obj, const Mesh *m);
int hFcn(double *obj, const Mesh *m);

int hFcn_d(double *obj, const Mesh *m);

int oMax(double *max, int *idx, const Mesh *m);

/*****************************************************************************/
/* Functions to obtain the global gradient or Hessian.  These assume that    */
/* the evaluation point is valid.  Only the gradient or Hessian is           */
/* calculated in an optimized procedure.                                     */
/*****************************************************************************/

void gOnly(const Mesh *m);
void hOnly(const Mesh *m);

/*****************************************************************************/
/* Functions to obtain the element objective, gradient, and Hessian.         */
/*****************************************************************************/

int o_fcn(double *obj, const double *x);
int o_fcn_gm(double *obj, const double *x);
int o_fcn_hm(double *obj, const double *x);

int g_fcn(double *obj, double *g_obj, const double *x);
int h_fcn(double *obj, double *g_obj, double *h_obj, const double *x);

#ifdef CHECK
int o_cond(double *obj, const double *x);
#endif

/*****************************************************************************/
/* Functions to obtain the element Hessian.  Assumes that the point is       */
/* valid.  No such routine for the gradient, since this does not save any    */
/* operations.                                                               */
/*****************************************************************************/

void h_only(double *h_obj, const double *x);

/*****************************************************************************/
/* Functions to obtain scaled element objective, gradient, and Hessian.      */
/*****************************************************************************/

int o_fcn_d(double *obj, const double *x, const double *d);
int g_fcn_d(double *obj, double *g_obj, const double *x, const double *d);
int h_fcn_d(double *obj, double *g_obj, double *h_obj, 
            const double *x, const double *d);

void h_only_d(double *h_obj, const double *x, const double *d);

/*****************************************************************************/
/* Functions to obtain the local objective, gradient, and Hessian.           */
/*   vert  : the vertex at the center of the local patch                     */
/*   elems : the indices of the elements in which the vertex is contained    */
/*   nelems: the number of such indices                                      */
/*****************************************************************************/

int oFcnl(double *obj, const int  vert,
          const int *elems, const int nelems, const Mesh *m);
int gFcnl(double *obj, double *g_obj, const int  vert,
          const int *elems, const int nelems, const Mesh *m);
int hFcnl(double *obj, double *g_obj, double *h_obj, const int  vert,
          const int *elems, const int nelems, const Mesh *m);

/*****************************************************************************/
/* Functions to obtain the element objective, gradient, and Hessian          */
/* specialized for the local objective.  Only the derivatives with respect   */
/* to the first coordinate are calculated.                                   */
/*****************************************************************************/

int g_fcnl_0(double *obj, double *g_obj, const double *x);
int g_fcnl_1(double *obj, double *g_obj, const double *x);
int g_fcnl_2(double *obj, double *g_obj, const double *x);
int g_fcnl_3(double *obj, double *g_obj, const double *x);

int h_fcnl_0(double *obj, double *g_obj, double *h_obj, const double *x);
int h_fcnl_1(double *obj, double *g_obj, double *h_obj, const double *x);
int h_fcnl_2(double *obj, double *g_obj, double *h_obj, const double *x);
int h_fcnl_3(double *obj, double *g_obj, double *h_obj, const double *x);

#endif

