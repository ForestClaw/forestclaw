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

#ifndef FC2D_CLAWPACK5_H
#define FC2D_CLAWPACK5_H

#include <fclaw2d_forestclaw.h>
#include <fclaw2d_clawpatch.h>
#include <fclaw2d_vtable.h>
#include <fclaw_package.h>

#include "fc2d_clawpack5_options.h"

#include "clawpack5_user_fort.h"


#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


typedef void (*fc2d_clawpack5_setprob_t)();

typedef void (*fc2d_clawpack5_bc2_t)(const int* meqn, const int* mbc,
                                     const int* mx, const int* my,
                                     const double* xlower, const double* ylower,
                                     const double* dx, const double* dy,
                                     const double q[], const int* maux,
                                     const double aux[], const double* t,
                                     const double* dt, const int mthbc[]);

typedef  void (*fc2d_clawpack5_qinit_t)(const int* meqn,const int* mbc,
                                        const int* mx, const int* my,
                                        const double* xlower, const double* ylower,
                                        const double* dx, const double* dy,
                                        double q[], const int* maux, double aux[]);

typedef void (*fc2d_clawpack5_setaux_t)(const int* mbc,
                                        const int* mx, const int* my,
                                        const double* xlower, const double* ylower,
                                        const double* dx, const double* dy,
                                        const int* maux, double aux[]);

typedef void (*fc2d_clawpack5_b4step2_t)(const int* mbc,
                                         const int* mx, const int* my, const int* meqn,
                                         double q[], const double* xlower,
                                         const double* ylower,
                                         const double* dx, const double* dy,
                                         const double* t, const double* dt,
                                         const int* maux, double aux[]);

typedef void (*fc2d_clawpack5_src2_t)(const int* meqn,
                                      const int* mbc, const int* mx,const int* my,
                                      const double* xlower, const double* ylower,
                                      const double* dx, const double* dy, double q[],
                                      const int* maux, double aux[], const double* t,
                                      const double* dt);


typedef void (*fc2d_clawpack5_rpn2_t)(const int* ixy,const int* maxm, const int* meqn,
                                      const int* mwaves, const int* maux,
                                      const int* mbc,const int* mx,
                                      double ql[], double qr[], double auxl[], double auxr[],
                                      double wave[], double s[],double amdq[], double apdq[]);

typedef void (*fc2d_clawpack5_rpt2_t)(const int* ixy, const int* imp, const int* maxm, const int* meqn,
                                       const int* mwaves, const int* maux, const int* mbc,const int* mx,
                                       double ql[], double qr[], double aux1[], double aux2[],
                                       double aux3[],  double asdq[],
                                       double bmasdq[], double bpasdq[]);


typedef void (*fc2d_clawpack5_flux2_t)(const int* ixy,const int* maxm, const int* meqn,
                                        const int* maux,const int* mbc,const int* mx,
                                        double q1d[], double dtdx1d[],
                                        double aux1[], double aux2[], double aux3[],
                                        double faddm[],double faddp[], double gaddm[],
                                        double gaddp[],double cfl1d[], double wave[],
                                        double s[], double amdq[],double apdq[],double cqxx[],
                                        double bmasdq[], double bpasdq[],
                                        fc2d_clawpack5_rpn2_t rpn2,
                                        fc2d_clawpack5_rpt2_t rpt2);

typedef void (*fc2d_clawpack5_fluxfun_t)(const int* meqn, double q[], double aux[],
                                          double fq[]);


typedef struct fc2d_clawpack5_vtable
{
    // fclaw2d_clawpatch_vtable_t  clawpatch_vt;

    fc2d_clawpack5_setprob_t setprob;
    fc2d_clawpack5_bc2_t bc2;
    fc2d_clawpack5_qinit_t qinit;
    fc2d_clawpack5_setaux_t setaux;
    fc2d_clawpack5_b4step2_t b4step2;
    fc2d_clawpack5_src2_t src2;
    fc2d_clawpack5_rpn2_t rpn2;
    fc2d_clawpack5_rpt2_t rpt2;
    fc2d_clawpack5_fluxfun_t fluxfun;
} fc2d_clawpack5_vtable_t;

#if 0
void fc2d_clawpack5_set_vtable(const fc2d_clawpack5_vtable_t vt);

void fc2d_clawpack5_set_vtable_defaults(fclaw2d_vtable_t *fclaw_vt,
                                        fc2d_clawpack5_vtable_t* vt);
#endif

void fc2d_clawpack5_set_vtable_defaults();

fc2d_clawpack5_vtable_t* fc2d_clawpack5_vt();

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
                            fc2d_clawpack5_rpn2_t rpn2,
                            fc2d_clawpack5_rpt2_t rpt2,
                            fc2d_clawpack5_flux2_t flux2,
                            int block_corner_count[],int* ierror);

#define CLAWPACK5_STEP2 FCLAW_F77_FUNC(clawpack5_step2,CLAWPACK5_STEP2)
void CLAWPACK5_STEP2(const int* maxm, const int* meqn, const int* maux,
                            const int* mbc, const int* mx,
                            const int* my, double qold[], double aux[],
                            const double* dx, const double* dy, const double* dt,
                            const double* cflgrid, double fm[], double fp[], double gm[],
                            double gp[],
                            fc2d_clawpack5_rpn2_t rpn2,
                            fc2d_clawpack5_rpt2_t rpt2);

#define CLAWPACK5_FLUX2 FCLAW_F77_FUNC(clawpack5_flux2,CLAWPACK5_FLUX2)
void CLAWPACK5_FLUX2(const int* ixy,const int* maxm, const int* meqn,
                      const int* maux,const int* mbc,const int* mx,
                      double q1d[], double dtdx1d[],
                      double aux1[], double aux2[], double aux3[],
                      double faddm[],double faddp[], double gaddm[],
                      double gaddp[],double cfl1d[], double wave[],
                      double s[], double amdq[],double apdq[],double cqxx[],
                      double bmasdq[], double bpasdq[],
                      fc2d_clawpack5_rpn2_t rpn2,fc2d_clawpack5_rpt2_t rpt2);
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
                        fc2d_clawpack5_rpn2_t rpn2,fc2d_clawpack5_rpt2_t rpt2,
                        const int* mwaves, const int* mcapa,
                        int method[], int mthlim[]);*/

#define CLAWPACK5_SET_CAPACITY FCLAW_F77_FUNC(clawpack5_set_capacity,CLAWPACK5_SET_CAPACITY)
void CLAWPACK5_SET_CAPACITY(const int* mx, const int *my, const int *mbc,
                             const double *dx, const double* dy, double area[],
                             const int *mcapa, const int* maux, double aux[]);

#define FC2D_CLAWPACK5_FORT_TAG4REFINEMENT FCLAW_F77_FUNC(fc2d_clawpack5_fort_tag4refinement, \
                                                          FC2D_CLAWPACK5_FORT_TAG4REFINEMENT)

void FC2D_CLAWPACK5_FORT_TAG4REFINEMENT(const int* mx,const int* my,
                                        const int* mbc,const int* meqn,
                                        const double* xlower, const double* ylower,
                                        const double* dx, const double* dy,
                                        const int* blockno,
                                        double q[],
                                        const double* tag_threshold,
                                        const int* init_flag,
                                        int* tag_patch);



#define FC2D_CLAWPACK5_FORT_TAG4COARSENING FCLAW_F77_FUNC(fc2d_clawpack5_fort_tag4coarsening, \
                                                          FC2D_CLAWPACK5_FORT_TAG4COARSENING)

void FC2D_CLAWPACK5_FORT_TAG4COARSENING(const int* mx, const int* my,
                                 const int* mbc, const int* meqn,
                                 const double* xlower, const double* ylower,
                                 const double* dx, const double* dy,
                                 const int* blockno,
                                 double q0[],double q1[],
                                 double q2[],double q3[],
                                 const double* tag_threshold,
                                 int* tag_patch);

#define FC2D_CLAWPACK5_FORT_INTERPOLATE2FINE FCLAW_F77_FUNC(fc2d_clawpack5_fort_interpolate2fine, \
                                                FC2D_CLAWPACK5_FORT_INTERPOLATE2FINE)
void FC2D_CLAWPACK5_FORT_INTERPOLATE2FINE(const int* mx,const int* my,
                                   const int* mbc, const int* meqn,
                                   double qcoarse[], double qfine[],
                                   double areacoarse[], double areafine[],
                                   const int* igrid, const int* manifold);

#define FC2D_CLAWPACK5_FORT_AVERAGE2COARSE FCLAW_F77_FUNC(fc2d_clawpack5_fort_average2coarse, \
                                                FC2D_CLAWPACK5_FORT_AVERAGE2COARSE)
void FC2D_CLAWPACK5_FORT_AVERAGE2COARSE(const int* mx, const int* my,
                                 const int* mbc, const int* meqn,
                                 double qcoarse[],double qfine[],
                                 double areacoarse[],double areafine[],
                                 const int* igrid, const int* manifold);

#define FC2D_CLAWPACK5_FORT_COPY_FACE FCLAW_F77_FUNC(fc2d_clawpack5_fort_copy_face, \
                                                     FC2D_CLAWPACK5_FORT_COPY_FACE)

void FC2D_CLAWPACK5_FORT_COPY_FACE(const int* mx, const int* my, const int* mbc, const int* meqn,
                                   double qthis[],double qneighbor[], const int* a_idir,
                                   fclaw2d_transform_data_t** transform_cptr);


#define FC2D_CLAWPACK5_FORT_AVERAGE_FACE FCLAW_F77_FUNC(fc2d_clawpack5_fort_average_face, \
                                                        FC2D_CLAWPACK5_FORT_AVERAGE_FACE)
void FC2D_CLAWPACK5_FORT_AVERAGE_FACE(const int* mx, const int* my, const int* mbc,
                                      const int* meqn,
                                      double qcoarse[],double qfine[],
                                      double areacoarse[], double areafine[],
                                      const int* idir, const int* iside,
                                      const int* num_neighbors,
                                      const int* refratio, const int* igrid,
                                      const int* manifold, fclaw2d_transform_data_t** transform_cptr);

#define FC2D_CLAWPACK5_FORT_INTERPOLATE_FACE FCLAW_F77_FUNC(fc2d_clawpack5_fort_interpolate_face, \
                                                            FC2D_CLAWPACK5_FORT_INTERPOLATE_FACE)
void FC2D_CLAWPACK5_FORT_INTERPOLATE_FACE(const int* mx, const int* my, const int* mbc,
                                          const int* meqn,
                                          double qcoarse[],double qfine[],
                                          const int* idir, const int* iside,
                                          const int* num_neighbors,
                                          const int* refratio, const int* igrid,
                                          fclaw2d_transform_data_t** transform_cptr);

#define FC2D_CLAWPACK5_FORT_COPY_CORNER FCLAW_F77_FUNC(fc2d_clawpack5_fort_copy_corner, \
                                                       FC2D_CLAWPACK5_FORT_COPY_CORNER)
void FC2D_CLAWPACK5_FORT_COPY_CORNER(const int* mx, const int* my, const int* mbc,
                                     const int* meqn, double this_q[],double neighbor_q[],
                                     const int* a_corner,fclaw2d_transform_data_t** transform_cptr);

#define FC2D_CLAWPACK5_FORT_AVERAGE_CORNER FCLAW_F77_FUNC(fc2d_clawpack5_fort_average_corner, \
                                                          FC2D_CLAWPACK5_FORT_AVERAGE_CORNER)
void FC2D_CLAWPACK5_FORT_AVERAGE_CORNER(const int* mx, const int* my, const int* mbc,
                                        const int* meqn, const int* a_refratio,
                                        double qcoarse[], double qfine[],
                                        double areacoarse[], double areafine[],
                                        const int* manifold,
                                        const int* a_corner, fclaw2d_transform_data_t** transform_cptr);

#define FC2D_CLAWPACK5_FORT_INTERPOLATE_CORNER FCLAW_F77_FUNC(fc2d_clawpack5_fort_interpolate_corner, \
                                                             FC2D_CLAWPACK5_FORT_INTERPOLATE_CORNER)
void FC2D_CLAWPACK5_FORT_INTERPOLATE_CORNER(const int* mx, const int* my, const int* mbc,
                                            const int* meqn, const int* a_refratio, double this_q[],
                                            double neighbor_q[], const int* a_corner,
                                            fclaw2d_transform_data_t** transform_cptr);

#define  FC2D_CLAWPACK5_FORT_WRITE_FILE FCLAW_F77_FUNC(fc2d_clawpack5_fort_write_file, \
                                                       FC2D_CLAWPACK5_FORT_WRITE_FILE)
void  FC2D_CLAWPACK5_FORT_WRITE_FILE(char* matname1,
                                     int* mx,        int* my,
                                     int* meqn,      int* mbc,
                                     double* xlower, double* ylower,
                                     double* dx,     double* dy,
                                     double q[],
                                     int* patch_num, int* level,
                                     int* blockno,   int* mpirank);

#define FC2D_CLAWPACK5_FORT_WRITE_HEADER FCLAW_F77_FUNC(fc2d_clawpack5_fort_write_header, \
                                                        FC2D_CLAWPACK5_FORT_WRITE_HEADER)
void FC2D_CLAWPACK5_FORT_WRITE_HEADER(int* iframe, double* time, int* meqn, int* maux,int* ngrids);


#define FC2D_CLAWPACK5_FORT_CONSERVATION_CHECK FCLAW_F77_FUNC(fc2d_clawpack5_fort_conservation_check, \
                                                              FC2D_CLAWPACK5_FORT_CONSERVATION_CHECK)
void FC2D_CLAWPACK5_FORT_CONSERVATION_CHECK(int *mx, int *my, int* mbc, int* meqn,
                                            double *dx, double *dy,
                                            double* area, double *q, double* sum);

#define FC2D_CLAWPACK5_FORT_COMPUTE_PATCH_AREA FCLAW_F77_FUNC(fc2d_clawpack5_fort_compute_patch_area, \
                                                              FC2D_CLAWPACK5_FORT_COMPUTE_PATCH_AREA)

double FC2D_CLAWPACK5_FORT_COMPUTE_PATCH_AREA(int *mx, int* my, int*mbc, double* dx,
                                              double* dy, double area[]);


#define FC2D_CLAWPACK5_FORT_COMPUTE_ERROR_NORM FCLAW_F77_FUNC(fc2d_clawpack5_fort_compute_error_norm, \
                                                              FC2D_CLAWPACK5_FORT_COMPUTE_ERROR_NORM)

void FC2D_CLAWPACK5_FORT_COMPUTE_ERROR_NORM(int *mx, int *my, int *mbc, int *meqn,
                                            double *dx, double *dy, double area[],
                                            double error[], double error_norm[]);

#define FC2D_CLAWPACK5_FORT_GHOSTPACK_QAREA FCLAW_F77_FUNC(fc2d_clawpack5_fort_ghostpack_qarea, \
                                                           FC2D_CLAWPACK5_FORT_GHOSTPACK_QAREA)
void  FC2D_CLAWPACK5_FORT_GHOSTPACK_QAREA(int *mx, int *my, int *mbc,
                                          int *meqn, int *mint,
                                          double qdata[], double area[],
                                          double qpack[], int *psize,
                                          int *packmode, int *ierror);

#define FC2D_CLAWPACK5_FORT_TIMEINTERP FCLAW_F77_FUNC (fc2d_clawpack5_fort_timeinterp, \
                                                       FC2D_CLAWPACK5_FORT_TIMEINTERP)
void FC2D_CLAWPACK5_FORT_TIMEINTERP(const int *mx, const int* my, const int* mbc,
                                    const int *meqn, const int* psize,
                                    double qcurr[], double qlast[],
                                    double qinterp[],const double* alpha,
                                    const int* ierror);

#define CLAWPACK5_SET_BLOCK FCLAW_F77_FUNC(clawpack5_set_block,CLAWPACK5_SET_BLOCK)
void CLAWPACK5_SET_BLOCK(int* blockno);

#define FC2D_CLAWPACK5_GET_BLOCK FCLAW_F77_FUNC(fc2d_clawpack5_get_block, \
                                                 FC2D_CLAWPACK5_GET_BLOCK)
int FC2D_CLAWPACK5_GET_BLOCK();


#define CLAWPACK5_UNSET_BLOCK FCLAW_F77_FUNC(clawpack5_unset_block, \
                                              CLAWPACK5_UNSET_BLOCK)
void CLAWPACK5_UNSET_BLOCK();

/***************************** MINIMAL API ******************************/

void fc2d_clawpack5_register_vtable (fclaw_package_container_t *
                                      pkg_container,
                                      fc2d_clawpack5_options_t *
                                      clawopt);

/* -------------------------------------------------------------------------
   New routines
   ------------------------------------------------------------------------- */
void fc2d_clawpack5_define_auxarray(fclaw2d_domain_t* domain,
                                    fclaw2d_patch_t* this_patch);

void fc2d_clawpack5_aux_data(fclaw2d_domain_t* domain,
                              fclaw2d_patch_t *this_patch,
                              double **aux, int* maux);

int fc2d_clawpack5_get_maux(fclaw2d_domain_t* domain);
void fc2d_clawpack5_maux(fclaw2d_domain_t* domain, int* maux);

void fc2d_clawpack5_register (fclaw_app_t* app, const char *configfile);

void fc2d_clawpack5_package_register(fclaw_app_t* app,
                                      fc2d_clawpack5_options_t* clawopt);

int fc2d_clawpack5_get_package_id (void);

void fc2d_clawpack5_output_header_ascii(fclaw2d_domain_t* domain,
                                        int iframe);

fc2d_clawpack5_options_t* fc2d_clawpack5_get_options(fclaw2d_domain_t *domain);

/* -------------------------------------------------------------------------
   Routines that won't change
   ------------------------------------------------------------------------- */
void
    fc2d_clawpack5_setprob(fclaw2d_domain_t* domain);

void
    fc2d_clawpack5_setaux(fclaw2d_domain_t *domain,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx);

void
    fc2d_clawpack5_set_capacity(fclaw2d_domain_t *domain,
                                 fclaw2d_patch_t *this_patch,
                                 int this_block_idx,
                                 int this_patch_idx);

void
    fc2d_clawpack5_qinit(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx);

void
    fc2d_clawpack5_b4step2(fclaw2d_domain_t *domain,
                            fclaw2d_patch_t *this_patch,
                            int this_block_idx,
                            int this_patch_idx,
                            double t,
                            double dt);

void
    fc2d_clawpack5_bc2(fclaw2d_domain_t *domain,
                        fclaw2d_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx,
                        double t,
                        double dt,
                        fclaw_bool intersects_bc[],
                        fclaw_bool time_interp);

void
    fc2d_clawpack5_src2(fclaw2d_domain_t *domain,
                         fclaw2d_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx,
                         double t,
                         double dt);


/* A single step method that advances the solution a single step on a single grid
   using a time step dt determined by the subcycle manager */
double
    fc2d_clawpack5_step2(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt);

/* Use this ro return only the right hand side of the clawpack algorithm */
double
    fc2d_clawpack5_step2_rhs(fclaw2d_domain_t *domain,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx,
                              double t,
                              double *rhs);

double
fc2d_clawpack5_update(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
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


#endif /* !FC2D_CLAWPACH5_H */
