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

#ifndef FCLAW2D_CLAWPATCH_H
#define FCLAW2D_CLAWPATCH_H

/* Needed for function pointer typedefs in fclaw2d_clawpatch_t */
#include <fclaw2d_patch.h>  /* Needed to get enum for build modes */
#include <fclaw2d_clawpatch_regrid.h>
#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_clawpatch_output.h> 

struct fclaw2d_transform_data;

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

typedef struct fclaw2d_clawpatch_vtable fclaw2d_clawpatch_vtable_t;

void fclaw2d_clawpatch_vtable_initialize();

fclaw2d_clawpatch_vtable_t* fclaw2d_clawpatch_vt();


void* fclaw2d_clawpatch_new();   /* Called in fclaw2d_patch */

void fclaw2d_clawpatch_delete(void *cp);


void fclaw2d_clawpatch_grid_data(struct fclaw2d_global* glob,
                                 struct fclaw2d_patch* this_patch,
                                 int* mx, int* my, int* mbc,
                                 double* xlower, double* ylower,
                                 double* dx, double* dy);


void fclaw2d_clawpatch_metric_data(struct fclaw2d_global* glob,
                                   struct fclaw2d_patch* this_patch,
                                   double **xp, double **yp, double **zp,
                                   double **xd, double **yd, double **zd,
                                   double **area);

void fclaw2d_clawpatch_metric_data2(struct fclaw2d_global* glob,
                                    struct fclaw2d_patch* this_patch,
                                    double **xnormals, double **ynormals,
                                    double **xtangents, double **ytangents,
                                    double **surfnormals, double ** edgelengths,
                                    double **curvature);

double* fclaw2d_clawpatch_get_area(struct fclaw2d_global* glob,
                                   struct fclaw2d_patch* this_patch);

void fclaw2d_clawpatch_soln_data(struct fclaw2d_global* glob,
                                 struct fclaw2d_patch* this_patch,
                                 double **q, int* meqn);

void fclaw2d_clawpatch_aux_data(struct fclaw2d_global *glob,
                                struct fclaw2d_patch *this_patch,
                                double **aux, int* maux);

double* fclaw2d_clawpatch_get_q(struct fclaw2d_global* glob,
                                struct fclaw2d_patch* this_patch);


double* fclaw2d_clawpatch_get_error(struct fclaw2d_global* glob,
                                    struct fclaw2d_patch* this_patch);

size_t fclaw2d_clawpatch_size(struct fclaw2d_global *glob);

void fclaw2d_clawpatch_timesync_data(struct fclaw2d_global* glob,
                                     struct fclaw2d_patch* this_patch,
                                     int time_interp,
                                     double **q, int* meqn);

double* fclaw2d_clawpatch_get_q_timesync(struct fclaw2d_global* glob,
                                         struct fclaw2d_patch* this_patch,
                                         int time_interp);

void fclaw2d_clawpatch_save_current_step(struct fclaw2d_global* glob,
                                         struct fclaw2d_patch* this_patch);

void fclaw2d_clawpatch_restore_step(struct fclaw2d_global* glob,
                                    struct fclaw2d_patch* this_patch);

void fclaw2d_clawpatch_save_step(struct fclaw2d_global* glob,
                                 struct fclaw2d_patch* this_patch);

#if 0
int* fclaw2d_clawpatch_block_corner_count(fclaw2d_domain_t* domain,
                                          struct fclaw2d_patch* this_patch);

void fclaw2d_clawpatch_set_block_corner_count(fclaw2d_domain_t* domain,
                                              struct fclaw2d_patch* this_patch,
                                              int icorner, int block_corner_count);
#endif

void fclaw2d_clawpatch_setup_timeinterp(struct fclaw2d_global *glob,
                                        struct fclaw2d_patch *this_patch,
                                        double alpha);
#if 0
void fclaw2d_clawpatch_finegrid_neighbors(fclaw2d_domain_t* domain);

int fclaw2d_clawpatch_has_finegrid_neighbors(fclaw2d_domain_t* domain,
                                             struct fclaw2d_patch* this_patch);
#endif

/* -----------------------------------------------------
   Define/build clawpatches
   ---------------------------------------------------- */

void fclaw2d_clawpatch_define(struct fclaw2d_global *glob,
                              struct fclaw2d_patch *this_patch,
                              int blockno, int patchno,
                              fclaw2d_build_mode_t build_mode);

/* A callback for building the domain and repartitioning */
void fclaw2d_clawpatch_build(struct fclaw2d_global *glob,
                             struct fclaw2d_patch *this_patch,
                             int this_block_idx,
                             int this_patch_idx,
                             void *user);

void fclaw2d_clawpatch_build_ghost(struct fclaw2d_global *glob,
                                   struct fclaw2d_patch *this_patch,
                                   int blockno,
                                   int patchno,
                                   void *user);


void fclaw2d_clawpatch_build_from_fine(struct fclaw2d_global *glob,
                                       struct fclaw2d_patch *fine_patches,
                                       struct fclaw2d_patch *coarse_patch,
                                       int blockno,
                                       int coarse_patchno,
                                       int fine0_patchno,
                                       fclaw2d_build_mode_t build_mode);

void fclaw2d_clawpatch_set_boundary_to_nan(fclaw2d_domain_t* domain,
                                           int minlevel,
                                           int maxlevel,
                                           int time_interp);

void fclaw2d_clawpatch_set_boundary_to_value(fclaw2d_domain_t* domain,
                                             int minlevel,
                                             int maxlevel,
                                             int time_interp,
                                             double value);

void fclaw2d_clawpatch_set_corners_to_value(fclaw2d_domain_t* domain,
                                            int minlevel,
                                            int maxlevel,
                                            int time_interp,
                                            double value);

void fclaw2d_clawpatch_set_corners_to_nan(fclaw2d_domain_t* domain,
                                          int minlevel,
                                          int maxlevel,
                                          int time_interp);


/* -----------------------------------------------------
   Build/pack/size for partitioning
   ---------------------------------------------------- */
size_t fclaw2d_clawpatch_partition_packsize(struct fclaw2d_global* glob);


void fclaw2d_clawpatch_partition_pack(struct fclaw2d_global *glob,
                                      struct fclaw2d_patch *this_patch,
                                      int this_block_idx,
                                      int this_patch_idx,
                                      void *user);

void fclaw2d_clawpatch_partition_unpack(struct fclaw2d_global *glob,
                                        fclaw2d_domain_t *new_domain,
                                        struct fclaw2d_patch *this_patch,
                                        int this_block_idx,
                                        int this_patch_idx,
                                        void *user);

/* -----------------------------------------------------
   Build/pack/size for ghost exchange
   ---------------------------------------------------- */
size_t fclaw2d_clawpatch_ghost_packsize(struct fclaw2d_global* glob);

void fclaw2d_clawpatch_ghost_unpack(struct fclaw2d_global* glob,
                                    struct fclaw2d_patch* this_patch,
                                    int this_block_idx, int this_patch_idx,
                                    double *qdata, int time_interp);


void fclaw2d_clawpatch_ghost_pack(struct fclaw2d_global *glob,
                                  struct fclaw2d_patch *this_patch,
                                  double *patch_data,
                                  int time_interp);

void fclaw2d_clawpatch_local_ghost_alloc(struct fclaw2d_global* glob,
                                         struct fclaw2d_patch* this_patch,
                                         void **q);

void fclaw2d_clawpatch_local_ghost_free(struct fclaw2d_global* glob,
                                        void **q);


void fclaw2d_clawpatch_copy_face(struct fclaw2d_global *glob,
                                 struct fclaw2d_patch *this_patch,
                                 struct fclaw2d_patch *neighbor_patch,
                                 int iface,
                                 int time_interp,
                                 struct fclaw2d_transform_data *transform_data);

void fclaw2d_clawpatch_average_face(struct fclaw2d_global *glob,
                                    struct fclaw2d_patch *coarse_patch,
                                    struct fclaw2d_patch *fine_patch,
                                    int idir,
                                    int iface_coarse,
                                    int p4est_refineFactor,
                                    int refratio,
                                    int time_interp,
                                    int igrid,
                                    struct fclaw2d_transform_data* transform_data);

void fclaw2d_clawpatch_interpolate_face(struct fclaw2d_global *glob,
                                        struct fclaw2d_patch *coarse_patch,
                                        struct fclaw2d_patch *fine_patch,
                                        int idir,
                                        int iside,
                                        int p4est_refineFactor,
                                        int refratio,
                                        int a_time_interp,
                                        int igrid,
                                        struct fclaw2d_transform_data* transform_data);



void fclaw2d_clawpatch_copy_corner(struct fclaw2d_global *glob,
                                   struct fclaw2d_patch *this_patch,
                                   struct fclaw2d_patch *corner_patch,
                                   int icorner,
                                   int time_interp,
                                   struct fclaw2d_transform_data *transform_data);

void fclaw2d_clawpatch_average_corner(struct fclaw2d_global *glob,
                                      struct fclaw2d_patch *coarse_patch,
                                      struct fclaw2d_patch *fine_patch,
                                      int coarse_corner,
                                      int refratio,
                                      int time_interp,
                                      struct fclaw2d_transform_data* transform_data);

void fclaw2d_clawpatch_interpolate_corner(struct fclaw2d_global* glob,
                                          struct fclaw2d_patch* coarse_patch,
                                          struct fclaw2d_patch* fine_patch,
                                          int coarse_corner,
                                          int refratio,
                                          int a_time_interp,
                                          struct fclaw2d_transform_data* transform_data);

typedef void (*fclaw2d_fort_copy_face_t)(const int* mx, const int* my, const int* mbc, 
                                         const int* meqn,
                                         double qthis[],double qneighbor[], 
                                         const int* a_idir,
                                         struct fclaw2d_transform_data** transform_cptr);

typedef void (*fclaw2d_fort_average_face_t)(const int* mx, const int* my, const int* mbc,
                                                const int* meqn,
                                                double qcoarse[],double qfine[],
                                                double areacoarse[], double areafine[],
                                                const int* idir, const int* iside,
                                                const int* num_neighbors,
                                                const int* refratio, const int* igrid,
                                                const int* manifold, 
                                                struct fclaw2d_transform_data** transform_cptr);

typedef void (*fclaw2d_fort_interpolate_face_t)(const int* mx, const int* my, const int* mbc,
                                                const int* meqn,
                                                double qcoarse[],double qfine[],
                                                const int* idir, const int* iside,
                                                const int* num_neighbors,
                                                const int* refratio, const int* igrid,
                                                struct fclaw2d_transform_data** transform_cptr);

typedef void (*fclaw2d_fort_copy_corner_t)(const int* mx, const int* my, const int* mbc,
                                     const int* meqn, double this_q[],double neighbor_q[],
                                     const int* a_corner,
                                     struct fclaw2d_transform_data** transform_cptr);

typedef void (*fclaw2d_fort_average_corner_t)(const int* mx, const int* my, const int* mbc,
                                        const int* meqn, const int* a_refratio,
                                        double qcoarse[], double qfine[],
                                        double areacoarse[], double areafine[],
                                        const int* manifold,
                                        const int* a_corner, 
                                        struct fclaw2d_transform_data** transform_cptr);

typedef void (*fclaw2d_fort_interpolate_corner_t)(const int* mx, const int* my, const int* mbc,
                                                  const int* meqn, const int* a_refratio, 
                                                  double this_q[],
                                                  double neighbor_q[], const int* a_corner,
                                                  struct fclaw2d_transform_data** transform_cptr);

typedef void (*fclaw2d_fort_timeinterp_t)(const int *mx, const int* my, const int* mbc,
                                          const int *meqn, const int* psize,
                                          double qcurr[], double qlast[],
                                          double qinterp[],const double* alpha,
                                          const int* ierror);

typedef void (*fclaw2d_fort_ghostpack_qarea_t)(int *mx, int *my, int *mbc,
                                               int *meqn, int *mint,
                                               double qdata[], double area[],
                                               double qpack[], int *psize,
                                               int *packmode, int *ierror);

typedef void (*fclaw2d_ghostpack_extra_t)(struct fclaw2d_global *glob,
                                          struct fclaw2d_patch *this_patch,
                                          int mint,
                                          double qpack[], int extrasize,
                                          int packmode, int* ierror);

/* Diagnostic functions */

typedef void (*fclaw2d_fort_error_t)(int* blockno, int *mx, int *my, int *mbc,int *meqn,
                                     double *dx, double *dy, double *xlower,
                                     double *ylower, double *t, double q[],
                                     double error[]);

typedef void (*fclaw2d_fort_conscheck_t)(int *mx, int *my, int* mbc, int* meqn,
                                         double *dx, double *dy,
                                         double area[], double q[], double sum[]);

typedef double (*fclaw2d_fort_area_t)(int *mx, int* my, int*mbc, double* dx,
                                      double* dy, double area[]);

typedef void (*fclaw2d_fort_norm_t)(int *mx, int *my, int *mbc,int *meqn,
                                    double *dx, double *dy, double area[],
                                    double error[], double error_norm[]);




void fclaw2d_clawpatch_diagnostics_initialize(struct fclaw2d_global *glob,
                                              void** patch_acc);

void fclaw2d_clawpatch_diagnostics_compute(struct fclaw2d_global *glob,
                                           void* patch_acc);

void fclaw2d_clawpatch_diagnostics_gather(struct fclaw2d_global *glob,
                                          void* patch_acc, int init_flag);

void fclaw2d_clawpatch_diagnostics_reset(struct fclaw2d_global *glob,
                                         void* patch_acc);

void fclaw2d_clawpatch_diagnostics_finalize(struct fclaw2d_global *glob,
                                            void** patch_acc);


/* -----------------------------------------------------
   Patch virtual table
   ---------------------------------------------------- */

struct fclaw2d_clawpatch_vtable
{

    /* These types should all have 'clawpatch' in the name */

    /* regridding functions */
    fclaw2d_fort_tag4refinement_t      fort_tag4refinement;
    fclaw2d_fort_tag4coarsening_t      fort_tag4coarsening;
    fclaw2d_fort_average2coarse_t      fort_average2coarse;
    fclaw2d_fort_interpolate2fine_t    fort_interpolate2fine;

    /* ghost filling functions */
    fclaw2d_fort_copy_face_t           fort_copy_face;
    fclaw2d_fort_average_face_t        fort_average_face;
    fclaw2d_fort_interpolate_face_t    fort_interpolate_face;
    fclaw2d_fort_copy_corner_t         fort_copy_corner;
    fclaw2d_fort_average_corner_t      fort_average_corner;
    fclaw2d_fort_interpolate_corner_t  fort_interpolate_corner;

    /* output functions (ascii) */
    fclaw2d_fort_header_ascii_t        fort_header_ascii;

    fclaw2d_patch_callback_t           cb_output_ascii;    
    fclaw2d_fort_output_ascii_t        fort_output_ascii;

    /* Time interpolation functions */
    fclaw2d_fort_timeinterp_t          fort_timeinterp;

    /* ghost patch functions */
    fclaw2d_fort_ghostpack_qarea_t     fort_ghostpack_qarea;
    fclaw2d_ghostpack_extra_t          ghostpack_extra;

    /* diagnostic functions */
    fclaw2d_fort_error_t       fort_compute_patch_error;
    fclaw2d_fort_conscheck_t   fort_conservation_check;
    fclaw2d_fort_norm_t        fort_compute_error_norm;
    fclaw2d_fort_area_t        fort_compute_patch_area;

    int is_set;
};

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW2D_CLAWPATCH_H */
