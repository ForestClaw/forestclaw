/*
  Copyright (c) 2012 Carsten Burstedde, Donna Calhoun
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef FC3D_CLAWPACK5_H
#define FC3D_CLAWPACK5_H

#include <fclaw_clawpatch3.h>
#include <fclaw2d_vtable.h>
#include <fclaw_package.h>
#include <fclaw2d_global.h>

#include "fc3d_clawpack5_options.h"

#include "clawpack5_user_fort.h"


#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


typedef void (*fc3d_clawpack5_setprob_t)(void);

typedef void (*fc3d_clawpack5_bc2_t)(const int* meqn, const int* mbc,
                                     const int* mx, const int* my,
                                     const double* xlower, const double* ylower,
                                     const double* dx, const double* dy,
                                     const double q[], const int* maux,
                                     const double aux[], const double* t,
                                     const double* dt, const int mthbc[]);

typedef  void (*fc3d_clawpack5_qinit_t)(const int* meqn,const int* mbc,
                                        const int* mx, const int* my, const int* mz,
                                        const double* xlower, const double* ylower, const double* zlower,
                                        const double* dx, const double* dy, const double* dz,
                                        double q[], const int* maux, double aux[]);

typedef void (*fc3d_clawpack5_setaux_t)(const int* mbc,
                                        const int* mx, const int* my, const int* mz,
                                        const double* xlower, const double* ylower,
                                        const double* zlower,
                                        const double* dx, const double* dy,
                                        const double* dz,
                                        const int* maux, double aux[]);

typedef void (*fc3d_clawpack5_b4step2_t)(const int* mbc,
                                         const int* mx, const int* my, const int* meqn,
                                         double q[], const double* xlower,
                                         const double* ylower,
                                         const double* dx, const double* dy,
                                         const double* t, const double* dt,
                                         const int* maux, double aux[]);

typedef void (*fc3d_clawpack5_src2_t)(const int* meqn,
                                      const int* mbc, const int* mx,const int* my,
                                      const double* xlower, const double* ylower,
                                      const double* dx, const double* dy, double q[],
                                      const int* maux, double aux[], const double* t,
                                      const double* dt);


typedef void (*fc3d_clawpack5_rpn2_t)(const int* ixy,const int* maxm, const int* meqn,
                                      const int* mwaves, const int* maux,
                                      const int* mbc,const int* mx,
                                      double ql[], double qr[], double auxl[], double auxr[],
                                      double wave[], double s[],double amdq[], double apdq[]);

typedef void (*fc3d_clawpack5_rpt2_t)(const int* ixy, const int* imp, const int* maxm, const int* meqn,
                                       const int* mwaves, const int* maux, const int* mbc,const int* mx,
                                       double ql[], double qr[], double aux1[], double aux2[],
                                       double aux3[],  double asdq[],
                                       double bmasdq[], double bpasdq[]);


typedef void (*fc3d_clawpack5_flux2_t)(const int* ixy,const int* maxm, const int* meqn,
                                        const int* maux,const int* mbc,const int* mx,
                                        double q1d[], double dtdx1d[],
                                        double aux1[], double aux2[], double aux3[],
                                        double faddm[],double faddp[], double gaddm[],
                                        double gaddp[],double cfl1d[], double wave[],
                                        double s[], double amdq[],double apdq[],double cqxx[],
                                        double bmasdq[], double bpasdq[],
                                        fc3d_clawpack5_rpn2_t rpn2,
                                        fc3d_clawpack5_rpt2_t rpt2);

typedef void (*fc3d_clawpack5_fluxfun_t)(const int* meqn, double q[], double aux[],
                                          double fq[]);


typedef struct fc3d_clawpack5_vtable
{
    fc3d_clawpack5_setprob_t setprob;
    fc3d_clawpack5_bc2_t bc2;
    fc3d_clawpack5_qinit_t qinit;
    fc3d_clawpack5_setaux_t setaux;
    fc3d_clawpack5_b4step2_t b4step2;
    fc3d_clawpack5_src2_t src2;
    fc3d_clawpack5_rpn2_t rpn2;
    fc3d_clawpack5_rpt2_t rpt2;
    fc3d_clawpack5_fluxfun_t fluxfun;
} fc3d_clawpack5_vtable_t;

#if 0
void fc3d_clawpack5_set_vtable(const fc3d_clawpack5_vtable_t vt);

void fc3d_clawpack5_set_vtable_defaults(fclaw2d_vtable_t *fclaw_vt,
                                        fc3d_clawpack5_vtable_t* vt);
#endif

void fc3d_clawpack5_set_vtable_defaults();

fc3d_clawpack5_vtable_t* fc3d_clawpack5_vt();

#define CLAWPACK5_BC2_DEFAULT FCLAW_F77_FUNC(clawpack5_bc2_default, \
                                             CLAWPACK5_BC2_DEFAULT)
void CLAWPACK5_BC2_DEFAULT(const int* meqn, const int* mbc,
                           const int* mx, const int* my,
                           const double* xlower, const double* ylower,
                           const double* dx, const double* dy,
                           const double q[], const int* maux,
                           const double aux[], const double* t,
                           const double* dt, const int mthbc[]);

/* --------------------------------------------------------------------
   Time stepping
   -------------------------------------------------------------------- */

#define CLAWPACK5_STEP2_WRAP FCLAW_F77_FUNC(clawpack5_step2_wrap,CLAWPACK5_STEP2_WRAP)
void CLAWPACK5_STEP2_WRAP(const int* maxm, const int* meqn, const int* maux,
                            const int* mbc, const int method[], const int mthlim[],
                            const int* mcapa, const int* mwaves, const int* mx,
                            const int* my, double qold[], double auxold[],
                            const double* dx, const double* dy, const double* dt,
                            const double* cfl, double work[], const int* mwork,
                            const double* xlower, const double* ylower, const int* level,
                            const double* t, double fp[], double fm[], double gp[],
                            double gm[],
                            fc3d_clawpack5_rpn2_t rpn2,
                            fc3d_clawpack5_rpt2_t rpt2,
                            fc3d_clawpack5_flux2_t flux2,
                            int block_corner_count[],int* ierror);

#define CLAWPACK5_STEP2 FCLAW_F77_FUNC(clawpack5_step2,CLAWPACK5_STEP2)
void CLAWPACK5_STEP2(const int* maxm, const int* meqn, const int* maux,
                            const int* mbc, const int* mx,
                            const int* my, double qold[], double aux[],
                            const double* dx, const double* dy, const double* dt,
                            const double* cflgrid, double fm[], double fp[], double gm[],
                            double gp[],
                            fc3d_clawpack5_rpn2_t rpn2,
                            fc3d_clawpack5_rpt2_t rpt2);

#define CLAWPACK5_FLUX2 FCLAW_F77_FUNC(clawpack5_flux2,CLAWPACK5_FLUX2)
void CLAWPACK5_FLUX2(const int* ixy,const int* maxm, const int* meqn,
                      const int* maux,const int* mbc,const int* mx,
                      double q1d[], double dtdx1d[],
                      double aux1[], double aux2[], double aux3[],
                      double faddm[],double faddp[], double gaddm[],
                      double gaddp[],double cfl1d[], double wave[],
                      double s[], double amdq[],double apdq[],double cqxx[],
                      double bmasdq[], double bpasdq[],
                      fc3d_clawpack5_rpn2_t rpn2,fc3d_clawpack5_rpt2_t rpt2);
/*
#define CLAWPACK5_FLUX2FW FCLAW_F77_FUNC(clawpack5_flux2fw,CLAWPACK5_FLUX2FW)
void CLAWPACK5_FLUX2FW(const int* ixy,const int* maxm, const int* meqn, //
                        const int* maux,const int* mbc,const int* mx,
                        double q1d[], double dtdx1d[],
                        double aux1[], double aux2[], double aux3[],
                        double faddm[],double faddp[], double gaddm[],
                        double gaddp[],double cfl1d[], double fwave[],
                        double s[], double amdq[],double apdq[],double cqxx[],
                        double bmasdq[], double bpasdq[],
                        fc3d_clawpack5_rpn2_t rpn2,fc3d_clawpack5_rpt2_t rpt2,
                        const int* mwaves, const int* mcapa,
                        int method[], int mthlim[]);*/

#define CLAWPACK5_SET_CAPACITY FCLAW_F77_FUNC(clawpack5_set_capacity,CLAWPACK5_SET_CAPACITY)
void CLAWPACK5_SET_CAPACITY(const int* mx, const int *my, const int *mbc,
                             const double *dx, const double* dy, double area[],
                             const int *mcapa, const int* maux, double aux[]);

#define FC3D_CLAWPACK5_FORT_TAG4REFINEMENT FCLAW_F77_FUNC(fc3d_clawpack5_fort_tag4refinement, \
                                                          FC3D_CLAWPACK5_FORT_TAG4REFINEMENT)

void FC3D_CLAWPACK5_FORT_TAG4REFINEMENT(const int* mx,const int* my,const int* mz,
                                        const int* mbc,const int* meqn,
                                        const double* xlower, const double* ylower,
                                        const double* zlower,
                                        const double* dx, const double* dy, const double* dz,
                                        const int* blockno,
                                        double q[],
                                        const double* tag_threshold,
                                        const int* init_flag,
                                        int* tag_patch);



#define FC3D_CLAWPACK5_FORT_TAG4COARSENING FCLAW_F77_FUNC(fc3d_clawpack5_fort_tag4coarsening, \
                                                          FC3D_CLAWPACK5_FORT_TAG4COARSENING)

void FC3D_CLAWPACK5_FORT_TAG4COARSENING(const int* mx, const int* my, const int* mz,
                                        const int* mbc, const int* meqn,
                                        const double* xlower, const double* ylower,
                                        const double* zlower, 
                                        const double* dx, const double* dy, const double* dz,
                                        const int* blockno,
                                        double q0[],double q1[],
                                        double q2[],double q3[],
                                        const double* tag_threshold,
                                        int* tag_patch);

#define FC3D_CLAWPACK5_FORT_INTERPOLATE2FINE FCLAW_F77_FUNC(fc3d_clawpack5_fort_interpolate2fine, \
                                                FC3D_CLAWPACK5_FORT_INTERPOLATE2FINE)
void FC3D_CLAWPACK5_FORT_INTERPOLATE2FINE(const int* mx,const int* my, const int* mz,
                                          const int* mbc, const int* meqn,
                                          double qcoarse[], double qfine[],
                                          double areacoarse[], double areafine[],
                                          const int* igrid, const int* manifold);

#define FC3D_CLAWPACK5_FORT_AVERAGE2COARSE FCLAW_F77_FUNC(fc3d_clawpack5_fort_average2coarse, \
                                                FC3D_CLAWPACK5_FORT_AVERAGE2COARSE)
void FC3D_CLAWPACK5_FORT_AVERAGE2COARSE(const int* mx, const int* my, const int* mz,
                                        const int* mbc, const int* meqn,
                                        double qcoarse[],double qfine[],
                                        double areacoarse[],double areafine[],
                                        const int* igrid, const int* manifold);

#define FC3D_CLAWPACK5_FORT_COPY_FACE FCLAW_F77_FUNC(fc3d_clawpack5_fort_copy_face, \
                                                     FC3D_CLAWPACK5_FORT_COPY_FACE)

void FC3D_CLAWPACK5_FORT_COPY_FACE(const int* mx, const int* my, const int* mz,
                                   const int* mbc, const int* meqn,
                                   double qthis[],double qneighbor[], const int* a_idir,
                                   fclaw2d_transform_data_t** transform_cptr);


#define FC3D_CLAWPACK5_FORT_AVERAGE_FACE FCLAW_F77_FUNC(fc3d_clawpack5_fort_average_face, \
                                                        FC3D_CLAWPACK5_FORT_AVERAGE_FACE)
void FC3D_CLAWPACK5_FORT_AVERAGE_FACE(const int* mx, const int* my, const int*mz,
                                      const int* mbc, const int* meqn,
                                      double qcoarse[],double qfine[],
                                      double areacoarse[], double areafine[],
                                      const int* idir, const int* iside,
                                      const int* num_neighbors,
                                      const int* refratio, const int* igrid,
                                      const int* manifold, fclaw2d_transform_data_t** transform_cptr);

#define FC3D_CLAWPACK5_FORT_INTERPOLATE_FACE FCLAW_F77_FUNC(fc3d_clawpack5_fort_interpolate_face, \
                                                            FC3D_CLAWPACK5_FORT_INTERPOLATE_FACE)
void FC3D_CLAWPACK5_FORT_INTERPOLATE_FACE(const int* mx, const int* my, const int* mz, const int* mbc,
                                          const int* meqn,
                                          double qcoarse[],double qfine[],
                                          const int* idir, const int* iside,
                                          const int* num_neighbors,
                                          const int* refratio, const int* igrid,
                                          fclaw2d_transform_data_t** transform_cptr);

#define FC3D_CLAWPACK5_FORT_COPY_CORNER FCLAW_F77_FUNC(fc3d_clawpack5_fort_copy_corner, \
                                                       FC3D_CLAWPACK5_FORT_COPY_CORNER)
void FC3D_CLAWPACK5_FORT_COPY_CORNER(const int* mx, const int* my, const int* mz, const int* mbc,
                                     const int* meqn, double this_q[],double neighbor_q[],
                                     const int* a_corner,fclaw2d_transform_data_t** transform_cptr);

#define FC3D_CLAWPACK5_FORT_AVERAGE_CORNER FCLAW_F77_FUNC(fc3d_clawpack5_fort_average_corner, \
                                                          FC3D_CLAWPACK5_FORT_AVERAGE_CORNER)
void FC3D_CLAWPACK5_FORT_AVERAGE_CORNER(const int* mx, const int* my, const int* mz, const int* mbc,
                                        const int* meqn, const int* a_refratio,
                                        double qcoarse[], double qfine[],
                                        double areacoarse[], double areafine[],
                                        const int* manifold,
                                        const int* a_corner, fclaw2d_transform_data_t** transform_cptr);

#define FC3D_CLAWPACK5_FORT_INTERPOLATE_CORNER FCLAW_F77_FUNC(fc3d_clawpack5_fort_interpolate_corner, \
                                                             FC3D_CLAWPACK5_FORT_INTERPOLATE_CORNER)
void FC3D_CLAWPACK5_FORT_INTERPOLATE_CORNER(const int* mx, const int* my, const int* mz, const int* mbc,
                                            const int* meqn, const int* a_refratio, double this_q[],
                                            double neighbor_q[], const int* a_corner,
                                            fclaw2d_transform_data_t** transform_cptr);

#define  FC3D_CLAWPACK5_FORT_WRITE_FILE FCLAW_F77_FUNC(fc3d_clawpack5_fort_write_file, \
                                                       FC3D_CLAWPACK5_FORT_WRITE_FILE)
void  FC3D_CLAWPACK5_FORT_WRITE_FILE(char* matname1,
                                     int* mx,        int* my,        int* mz,
                                     int* meqn,      int* mbc,
                                     double* xlower, double* ylower, double* zlower,
                                     double* dx,     double* dy,     double* dz,
                                     double q[],
                                     int* patch_num, int* level,
                                     int* blockno,   int* mpirank);

#define FC3D_CLAWPACK5_FORT_WRITE_HEADER FCLAW_F77_FUNC(fc3d_clawpack5_fort_write_header, \
                                                        FC3D_CLAWPACK5_FORT_WRITE_HEADER)
void FC3D_CLAWPACK5_FORT_WRITE_HEADER(int* iframe, double* time, int* meqn, int* maux,int* ngrids);


#define FC3D_CLAWPACK5_FORT_CONSERVATION_CHECK FCLAW_F77_FUNC(fc3d_clawpack5_fort_conservation_check, \
                                                              FC3D_CLAWPACK5_FORT_CONSERVATION_CHECK)
void FC3D_CLAWPACK5_FORT_CONSERVATION_CHECK(int *mx, int *my, int* mz, int* mbc, int* meqn,
                                            double *dx, double *dy, double *dz,
                                            double* area, double *q, double* sum);

#define FC3D_CLAWPACK5_FORT_COMPUTE_PATCH_AREA FCLAW_F77_FUNC(fc3d_clawpack5_fort_compute_patch_area, \
                                                              FC3D_CLAWPACK5_FORT_COMPUTE_PATCH_AREA)

double FC3D_CLAWPACK5_FORT_COMPUTE_PATCH_AREA(int *mx, int* my, int* mz, int*mbc, double* dx,
                                              double* dy, double* dz, double area[]);


#define FC3D_CLAWPACK5_FORT_COMPUTE_ERROR_NORM FCLAW_F77_FUNC(fc3d_clawpack5_fort_compute_error_norm, \
                                                              FC3D_CLAWPACK5_FORT_COMPUTE_ERROR_NORM)

void FC3D_CLAWPACK5_FORT_COMPUTE_ERROR_NORM(int *mx, int *my, int *mz, int *mbc, int *meqn,
                                            double *dx, double *dy, double *dz, double area[],
                                            double error[], double error_norm[]);

#if 0
#define FC3D_CLAWPACK5_FORT_COMPUTE_ERROR FCLAW_F77_FUNC(fc3d_clawpack5_fort_compute_error, \
                                                         FC3D_CLAWPACK5_FORT_COMPUTE_ERROR)

void FC3D_CLAWPACK5_FORT_COMPUTE_ERROR(int *mx, int *my, int *mbc, int *meqn,
                                       double *dx, double *dy, double area[],
                                       double error[], double error_norm[]);

#endif

#define FC3D_CLAWPACK5_FORT_GHOSTPACK_QAREA FCLAW_F77_FUNC(fc3d_clawpack5_fort_ghostpack_qarea, \
                                                           FC3D_CLAWPACK5_FORT_GHOSTPACK_QAREA)
void  FC3D_CLAWPACK5_FORT_GHOSTPACK_QAREA(int *mx, int *my, int* mz, int *mbc,
                                          int *meqn, int *mint,
                                          double qdata[], double area[],
                                          double qpack[], int *psize,
                                          int *packmode, int *ierror);

#define FC3D_CLAWPACK5_FORT_TIMEINTERP FCLAW_F77_FUNC (fc3d_clawpack5_fort_timeinterp, \
                                                       FC3D_CLAWPACK5_FORT_TIMEINTERP)
void FC3D_CLAWPACK5_FORT_TIMEINTERP(const int *mx, const int* my, const int* mz, const int* mbc,
                                    const int *meqn, const int* psize,
                                    double qcurr[], double qlast[],
                                    double qinterp[],const double* alpha,
                                    const int* ierror);

#define CLAWPACK5_SET_BLOCK FCLAW_F77_FUNC(clawpack5_set_block,CLAWPACK5_SET_BLOCK)
void CLAWPACK5_SET_BLOCK(int* blockno);

#define FC3D_CLAWPACK5_GET_BLOCK FCLAW_F77_FUNC(fc3d_clawpack5_get_block, \
                                                FC3D_CLAWPACK5_GET_BLOCK)
int FC3D_CLAWPACK5_GET_BLOCK();


#define CLAWPACK5_UNSET_BLOCK FCLAW_F77_FUNC(clawpack5_unset_block, \
                                              CLAWPACK5_UNSET_BLOCK)
void CLAWPACK5_UNSET_BLOCK();

/***************************** MINIMAL API ******************************/

void fc3d_clawpack5_register_vtable (fclaw_package_container_t *
                                      pkg_container,
                                      fc3d_clawpack5_options_t *
                                      clawopt);

/* -------------------------------------------------------------------------
   New routines
   ------------------------------------------------------------------------- */
void fc3d_clawpack5_aux_data(fclaw2d_global_t *glob,
                              fclaw_patch_t *this_patch,
                              double **aux, int* maux);

// void fc3d_clawpack5_register (fclaw_app_t* app, const char *configfile, fclaw2d_global_t* glob);
void fc3d_clawpack5_set_options (fclaw2d_global_t* glob, fc3d_clawpack5_options_t* clawopt);

void fc3d_clawpack5_package_register(fclaw_app_t* app,
                                      fc3d_clawpack5_options_t* clawopt);

void fc3d_clawpack5_output_header_ascii(fclaw2d_global_t* glob,
                                        int iframe);

fc3d_clawpack5_options_t* fc3d_clawpack5_get_options(fclaw2d_global_t *glob);

/* -------------------------------------------------------------------------
   Routines that won't change
   ------------------------------------------------------------------------- */
void
    fc3d_clawpack5_setprob(fclaw2d_global_t* glob);

void
    fc3d_clawpack5_setaux(fclaw2d_global_t *glob,
                           fclaw_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx);

void
    fc3d_clawpack5_set_capacity(fclaw2d_global_t *glob,
                                 fclaw_patch_t *this_patch,
                                 int this_block_idx,
                                 int this_patch_idx);

void
    fc3d_clawpack5_qinit(fclaw2d_global_t *glob,
                          fclaw_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx);

void
    fc3d_clawpack5_b4step2(fclaw2d_global_t *glob,
                            fclaw_patch_t *this_patch,
                            int this_block_idx,
                            int this_patch_idx,
                            double t,
                            double dt);

void
    fc3d_clawpack5_bc2(fclaw2d_global_t *glob,
                        fclaw_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx,
                        double t,
                        double dt,
                        fclaw_bool intersects_bc[],
                        fclaw_bool time_interp);

void
    fc3d_clawpack5_src2(fclaw2d_global_t *glob,
                         fclaw_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx,
                         double t,
                         double dt);


/* A single step method that advances the solution a single step on a single grid
   using a time step dt determined by the subcycle manager */
double fc3d_clawpack5_step2(fclaw2d_global_t *glob,
                            fclaw_patch_t *this_patch,
                            int this_block_idx,
                            int this_patch_idx,
                            double t,
                            double dt);

/* Use this ro return only the right hand side of the clawpack algorithm */
double
    fc3d_clawpack5_step2_rhs(fclaw2d_global_t *glob,
                              fclaw_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx,
                              double t,
                              double *rhs);

double fc3d_clawpack5_update(fclaw2d_global_t *glob,
                             fclaw_patch_t *this_patch,
                             int this_block_idx,
                             int this_patch_idx,
                             double t,
                             double dt);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif /* !FC3D_CLAWPACH5_H */
