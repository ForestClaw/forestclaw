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

#ifndef FCLAW_CLAWPATCH3_H
#define FCLAW_CLAWPATCH3_H

#include <fclaw2d_patch.h>
#include <fclaw2d_global.h>
#include <fclaw2d_domain.h>

#include <fclaw_clawpatch3_options.h>
#include <fclaw_clawpatch3_regrid.h>
#include <fclaw_clawpatch3_output_ascii.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

typedef struct fclaw_clawpatch3_vtable fclaw_clawpatch3_vtable_t;

void* fclaw_clawpatch3_new_patch();   /* Called in fclaw2d_patch */
void fclaw_clawpatch3_delete_patch(void *cp);


void fclaw_clawpatch3_grid_data(fclaw2d_global_t* glob,
                                 fclaw2d_patch_t* this_patch,
                                 int* mx, int* my, int* mz, int* mbc,
                                 double* xlower, double* ylower, double* zlower,
                                 double* dx, double* dy, double* dz);


void fclaw_clawpatch3_metric_data(fclaw2d_global_t* glob,
                                   fclaw2d_patch_t* this_patch,
                                   double **xp, double **yp, double **zp,
                                   double **xd, double **yd, double **zd,
                                   double **area);

void fclaw_clawpatch3_metric_data2(fclaw2d_global_t* glob,
                                    fclaw2d_patch_t* this_patch,
                                    double **xnormals, double **ynormals,
                                    double **xtangents, double **ytangents,
                                    double **surfnormals, double ** edgelengths,
                                    double **curvature);

double* fclaw_clawpatch3_get_area(fclaw2d_global_t* glob,
                                   fclaw2d_patch_t* this_patch);

void fclaw_clawpatch3_soln_data(fclaw2d_global_t* glob,
                                 fclaw2d_patch_t* this_patch,
                                 double **q, int* meqn);

void fclaw_clawpatch3_aux_data(fclaw2d_global_t *glob,
                                fclaw2d_patch_t *this_patch,
                                double **aux, int* maux);

double* fclaw_clawpatch3_get_q(fclaw2d_global_t* glob,
                                fclaw2d_patch_t* this_patch);


double* fclaw_clawpatch3_get_error(fclaw2d_global_t* glob,
                                    fclaw2d_patch_t* this_patch);

size_t fclaw_clawpatch3_size(fclaw2d_global_t *glob);

void fclaw_clawpatch3_timesync_data(fclaw2d_global_t* glob,
                                     fclaw2d_patch_t* this_patch,
                                     fclaw_bool time_interp,
                                     double **q, int* meqn);

double* fclaw_clawpatch3_get_q_timesync(fclaw2d_global_t* glob,
                                         fclaw2d_patch_t* this_patch,
                                         int time_interp);

void fclaw_clawpatch3_save_current_step(fclaw2d_global_t* glob,
                                         fclaw2d_patch_t* this_patch);

void fclaw_clawpatch3_restore_step(fclaw2d_global_t* glob,
                                    fclaw2d_patch_t* this_patch);

void fclaw_clawpatch3_save_step(fclaw2d_global_t* glob,
                                 fclaw2d_patch_t* this_patch);

#if 0
int* fclaw_clawpatch3_block_corner_count(fclaw2d_domain_t* domain,
                                          fclaw2d_patch_t* this_patch);

void fclaw_clawpatch3_set_block_corner_count(fclaw2d_domain_t* domain,
                                              fclaw2d_patch_t* this_patch,
                                              int icorner, int block_corner_count);
#endif

void fclaw_clawpatch3_setup_timeinterp(fclaw2d_global_t *glob,
                                        fclaw2d_patch_t *this_patch,
                                        double alpha);
#if 0
void fclaw_clawpatch3_finegrid_neighbors(fclaw2d_domain_t* domain);

int fclaw_clawpatch3_has_finegrid_neighbors(fclaw2d_domain_t* domain,
                                             fclaw2d_patch_t* this_patch);
#endif

/* -----------------------------------------------------
   Define/build clawpatches
   ---------------------------------------------------- */

void fclaw_clawpatch3_define(fclaw2d_global_t *glob,
                              fclaw2d_patch_t *this_patch,
                              int blockno, int patchno,
                              fclaw2d_build_mode_t build_mode);

/* A callback for building the domain and repartitioning */
void fclaw_clawpatch3_build(fclaw2d_global_t *glob,
                             fclaw2d_patch_t *this_patch,
                             int this_block_idx,
                             int this_patch_idx,
                             void *user);

void fclaw_clawpatch3_build_ghost(fclaw2d_global_t *glob,
                                   fclaw2d_patch_t *this_patch,
                                   int blockno,
                                   int patchno,
                                   void *user);


void fclaw_clawpatch3_build_from_fine(fclaw2d_global_t *glob,
                                       fclaw2d_patch_t *fine_patches,
                                       fclaw2d_patch_t *coarse_patch,
                                       int blockno,
                                       int coarse_patchno,
                                       int fine0_patchno,
                                       fclaw2d_build_mode_t build_mode);

void fclaw_clawpatch3_set_boundary_to_nan(fclaw2d_domain_t* domain,
                                           int minlevel,
                                           int maxlevel,
                                           int time_interp);

void fclaw_clawpatch3_set_boundary_to_value(fclaw2d_domain_t* domain,
                                             int minlevel,
                                             int maxlevel,
                                             int time_interp,
                                             double value);

void fclaw_clawpatch3_set_corners_to_value(fclaw2d_domain_t* domain,
                                            int minlevel,
                                            int maxlevel,
                                            int time_interp,
                                            double value);

void fclaw_clawpatch3_set_corners_to_nan(fclaw2d_domain_t* domain,
                                          int minlevel,
                                          int maxlevel,
                                          int time_interp);


/* -----------------------------------------------------
   Build/pack/size for partitioning
   ---------------------------------------------------- */
size_t fclaw_clawpatch3_partition_packsize(fclaw2d_global_t* glob);


void fclaw_clawpatch3_partition_pack(fclaw2d_global_t *glob,
                                      fclaw2d_patch_t *this_patch,
                                      int this_block_idx,
                                      int this_patch_idx,
                                      void *user);

void fclaw_clawpatch3_partition_unpack(fclaw2d_global_t *glob,
                                        fclaw2d_domain_t *new_domain,
                                        fclaw2d_patch_t *this_patch,
                                        int this_block_idx,
                                        int this_patch_idx,
                                        void *user);
#if 0
void fclaw_clawpatch3_initialize_after_partition(fclaw2d_domain_t* domain,
                                                  fclaw2d_patch_t* this_patch,
                                                  int this_block_idx,
                                                  int this_patch_idx);
void fclaw_clawpatch3_initialize_after_regrid(fclaw2d_domain_t* domain,
                                               fclaw2d_patch_t* this_patch,
                                               int this_block_idx,
                                               int this_patch_idx);
#endif


/* -----------------------------------------------------
   Build/pack/size for ghost exchange
   ---------------------------------------------------- */
size_t fclaw_clawpatch3_ghost_packsize(fclaw2d_global_t* glob);

void fclaw_clawpatch3_ghost_unpack(fclaw2d_global_t* glob,
                                    fclaw2d_patch_t* this_patch,
                                    int this_block_idx, int this_patch_idx,
                                    double *qdata, fclaw_bool time_interp);


void fclaw_clawpatch3_ghost_pack(fclaw2d_global_t *glob,
                                  fclaw2d_patch_t *this_patch,
                                  double *patch_data,
                                  int time_interp);

void fclaw_clawpatch3_local_ghost_alloc(fclaw2d_global_t* glob,
                                         fclaw2d_patch_t* this_patch,
                                         void **q);

void fclaw_clawpatch3_local_ghost_free(fclaw2d_global_t* glob,
                                        void **q);


#if 0
void fclaw_clawpatch3_ghost_comm(fclaw2d_domain_t* domain,
                                  fclaw2d_patch_t* this_patch,
                                  double *qpack, int time_interp,
                                  int packmode);
#endif

void fclaw_clawpatch3_init_vtable_defaults();

void fclaw_clawpatch3_copy_face(fclaw2d_global_t *glob,
                                 fclaw2d_patch_t *this_patch,
                                 fclaw2d_patch_t *neighbor_patch,
                                 int iface,
                                 int time_interp,
                                 fclaw2d_transform_data_t *transform_data);

void fclaw_clawpatch3_average_face(fclaw2d_global_t *glob,
                                    fclaw2d_patch_t *coarse_patch,
                                    fclaw2d_patch_t *fine_patch,
                                    int idir,
                                    int iface_coarse,
                                    int p4est_refineFactor,
                                    int refratio,
                                    fclaw_bool time_interp,
                                    int igrid,
                                    fclaw2d_transform_data_t* transform_data);

void fclaw_clawpatch3_interpolate_face(fclaw2d_global_t *glob,
                                        fclaw2d_patch_t *coarse_patch,
                                        fclaw2d_patch_t *fine_patch,
                                        int idir,
                                        int iside,
                                        int p4est_refineFactor,
                                        int refratio,
                                        fclaw_bool a_time_interp,
                                        int igrid,
                                        fclaw2d_transform_data_t* transform_data);



void fclaw_clawpatch3_copy_corner(fclaw2d_global_t *glob,
                                   fclaw2d_patch_t *this_patch,
                                   fclaw2d_patch_t *corner_patch,
                                   int icorner,
                                   int time_interp,
                                   fclaw2d_transform_data_t *transform_data);

void fclaw_clawpatch3_average_corner(fclaw2d_global_t *glob,
                                      fclaw2d_patch_t *coarse_patch,
                                      fclaw2d_patch_t *fine_patch,
                                      int coarse_corner,
                                      int refratio,
                                      fclaw_bool time_interp,
                                      fclaw2d_transform_data_t* transform_data);

void fclaw_clawpatch3_interpolate_corner(fclaw2d_global_t* glob,
                                          fclaw2d_patch_t* coarse_patch,
                                          fclaw2d_patch_t* fine_patch,
                                          int coarse_corner,
                                          int refratio,
                                          fclaw_bool a_time_interp,
                                          fclaw2d_transform_data_t* transform_data);

typedef void (*fclaw_fort_copy_face_t)(const int* mx, const int* my, const int* mz, 
                                         const int* mbc, const int* meqn,
                                         double qthis[],double qneighbor[], const int* a_idir,
                                         fclaw2d_transform_data_t** transform_cptr);

typedef void (*fclaw_fort_average_face_t)(const int* mx, const int* my, const int* mz, const int* mbc,
                                            const int* meqn,
                                            double qcoarse[],double qfine[],
                                            double areacoarse[], double areafine[],
                                            const int* idir, const int* iside,
                                            const int* num_neighbors,
                                            const int* refratio, const int* igrid,
                                            const int* manifold, fclaw2d_transform_data_t** transform_cptr);

typedef void (*fclaw_fort_interpolate_face_t)(const int* mx, const int* my, const int* mz, const int* mbc,
                                                const int* meqn,
                                                double qcoarse[],double qfine[],
                                                const int* idir, const int* iside,
                                                const int* num_neighbors,
                                                const int* refratio, const int* igrid,
                                                fclaw2d_transform_data_t** transform_cptr);

typedef void (*fclaw_fort_copy_corner_t)(const int* mx, const int* my, const int* mz, const int* mbc,
                                           const int* meqn, double this_q[],double neighbor_q[],
                                           const int* a_corner,fclaw2d_transform_data_t** transform_cptr);

typedef void (*fclaw_fort_average_corner_t)(const int* mx, const int* my, const int* mz, const int* mbc,
                                              const int* meqn, const int* a_refratio,
                                              double qcoarse[], double qfine[],
                                              double areacoarse[], double areafine[],
                                              const int* manifold,
                                              const int* a_corner, fclaw2d_transform_data_t** transform_cptr);

typedef void (*fclaw_fort_interpolate_corner_t)(const int* mx, const int* my, const int* mz, const int* mbc,
                                                  const int* meqn, const int* a_refratio, double this_q[],
                                                  double neighbor_q[], const int* a_corner,
                                                  fclaw2d_transform_data_t** transform_cptr);

typedef void (*fclaw_fort_timeinterp_t)(const int *mx, const int* my, const int* mz, const int* mbc,
                                          const int *meqn, const int* psize,
                                          double qcurr[], double qlast[],
                                          double qinterp[],const double* alpha,
                                          const int* ierror);

typedef void (*fclaw_fort_ghostpack_qarea_t)(int *mx, int *my, int *mz, int *mbc,
                                               int *meqn, int *mint,
                                               double qdata[], double area[],
                                               double qpack[], int *psize,
                                               int *packmode, int *ierror);

typedef void (*fclaw_ghostpack_extra_t)(fclaw2d_global_t *glob,
                                          fclaw2d_patch_t *this_patch,
                                          int mint,
                                          double qpack[], int extrasize,
                                          int packmode, int* ierror);

/* Diagnostic functions */

typedef void (*fclaw_fort_error_t)(int* blockno, int *mx, int *my, int* mz, int *mbc,int *meqn,
                                     double *dx, double *dy, double *dz, double *xlower,
                                     double *ylower, double *zlower, double *t, double q[],
                                     double error[]);

typedef void (*fclaw_fort_conscheck_t)(int *mx, int *my, int* mz, int* mbc, int* meqn,
                                         double *dx, double *dy, double *dz,
                                         double area[], double q[], double sum[]);

typedef double (*fclaw_fort_area_t)(int *mx, int* my, int* mz, int*mbc, double* dx,
                                      double* dy, double *dz, double area[]);

typedef void (*fclaw_fort_norm_t)(int *mx, int *my, int* mz, int *mbc,int *meqn,
                                    double *dx, double *dy, double *dz, double area[],
                                    double error[], double error_norm[]);




void fclaw_clawpatch3_diagnostics_initialize(fclaw2d_global_t *glob,
                                              void** patch_acc);

void fclaw_clawpatch3_diagnostics_compute(fclaw2d_global_t *glob,
                                           void* patch_acc);

void fclaw_clawpatch3_diagnostics_gather(fclaw2d_global_t *glob,
                                          void* patch_acc, int init_flag);

void fclaw_clawpatch3_diagnostics_reset(fclaw2d_global_t *glob,
                                         void* patch_acc);

void fclaw_clawpatch3_diagnostics_finalize(fclaw2d_global_t *glob,
                                            void** patch_acc);



struct fclaw_clawpatch3_vtable
{
    /* regridding functions */
    fclaw_fort_tag4refinement_t      fort_tag4refinement;
    fclaw_fort_tag4coarsening_t      fort_tag4coarsening;
    fclaw_fort_average2coarse_t      fort_average2coarse;
    fclaw_fort_interpolate2fine_t    fort_interpolate2fine;

    /* ghost filling functions */
    fclaw_fort_copy_face_t           fort_copy_face;
    fclaw_fort_average_face_t        fort_average_face;
    fclaw_fort_interpolate_face_t    fort_interpolate_face;
    fclaw_fort_copy_corner_t         fort_copy_corner;
    fclaw_fort_average_corner_t      fort_average_corner;
    fclaw_fort_interpolate_corner_t  fort_interpolate_corner;

    /* output functions */
    fclaw_fort_write_header_t        fort_write_header;
    fclaw_fort_write_file_t          fort_write_file;

    /* Time interpolation functions */
    fclaw_fort_timeinterp_t          fort_timeinterp;

    /* ghost patch functions */
    fclaw_fort_ghostpack_qarea_t     fort_ghostpack_qarea;
    fclaw_ghostpack_extra_t          ghostpack_extra;

    /* diagnostic functions */
    fclaw_fort_error_t               fort_compute_patch_error;
    fclaw_fort_conscheck_t           fort_conservation_check;
    fclaw_fort_norm_t                fort_compute_error_norm;
    fclaw_fort_area_t                fort_compute_patch_area;

    int defaults_set;
};

fclaw_clawpatch3_vtable_t* fclaw_clawpatch3_vt();

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW2D_CLAWPATCH_H */
