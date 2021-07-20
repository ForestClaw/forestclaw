/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun
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

#ifndef FCLAW2D_CLAWPATCH46_FORT3_H
#define FCLAW2D_CLAWPATCH46_FORT3_H

#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix syntax highlighting */
#endif



/* ------------------------------ Time stepping functions ----------------------------- */
#define FCLAW2D_CLAWPATCH46_FORT3_TIMEINTERP \
            FCLAW_F77_FUNC (fclaw2d_clawpatch46_fort3_timeinterp, \
                            FCLAW2D_CLAWPATCH46_FORT3_TIMEINTERP)
void FCLAW2D_CLAWPATCH46_FORT3_TIMEINTERP(const int *mx, const int* my, const int *mz,
                                         const int* mbc,
                                         const int *meqn, const int* psize,
                                         double qcurr[], double qlast[],
                                         double qinterp[],const double* alpha,
                                         const int* ierror);

/* --------------------------------- Regridding functions ----------------------------- */

#define FCLAW2D_CLAWPATCH46_FORT3_TAG4REFINEMENT \
           FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort3_tag4refinement, \
                          FCLAW2D_CLAWPATCH46_FORT3_TAG4REFINEMENT)

void FCLAW2D_CLAWPATCH46_FORT3_TAG4REFINEMENT(const int* mx,const int* my,const int *mz,
                                             const int* mbc,const int* meqn,
                                             const double* xlower, const double* ylower,
                                             const double* zlower,
                                             const double* dx, const double* dy,
                                             const double* dz,
                                             const int* blockno,
                                             double q[],
                                             const double* tag_threshold,
                                             const int* init_flag,
                                             int* tag_patch);



#define FCLAW2D_CLAWPATCH46_FORT3_TAG4COARSENING \
           FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort3_tag4coarsening, \
                          FCLAW2D_CLAWPATCH46_FORT3_TAG4COARSENING)

void FCLAW2D_CLAWPATCH46_FORT3_TAG4COARSENING(const int* mx, const int* my,
                                              const int *mz,
                                              const int* mbc, const int* meqn,
                                              double xlower[], double ylower[],
                                              const double* zlower,
                                              const double* dx, const double* dy,
                                              const double* dz,
                                              const int* blockno,
                                              double q0[], double q1[],
                                              double q2[], double q3[],
                                              const double* tag_threshold,
                                              const int* initflag,
                                              int* tag_patch);

#define FCLAW2D_CLAWPATCH46_FORT3_INTERPOLATE2FINE \
           FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort3_interpolate2fine, \
                          FCLAW2D_CLAWPATCH46_FORT3_INTERPOLATE2FINE)
void FCLAW2D_CLAWPATCH46_FORT3_INTERPOLATE2FINE(const int* mx,const int* my,const int* mz,
                                               const int* mbc, const int* meqn,
                                               double qcoarse[], double qfine[],
                                               double areacoarse[], double areafine[],
                                               const int* igrid, const int* manifold);
  
#define FCLAW2D_CLAWPATCH46_FORT3_AVERAGE2COARSE \
      FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort3_average2coarse, \
                     FCLAW2D_CLAWPATCH46_FORT3_AVERAGE2COARSE)
void FCLAW2D_CLAWPATCH46_FORT3_AVERAGE2COARSE(const int* mx, const int* my,const int* mz,
                                             const int* mbc, const int* meqn,
                                             double qcoarse[],double qfine[],
                                             double areacoarse[],double areafine[],
                                             const int* igrid, const int* manifold);



/* ---------------------------------- Ghost filling  ---------------------------------- */

#define FCLAW2D_CLAWPATCH46_FORT3_COPY_FACE \
         FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort3_copy_face, \
                        FCLAW2D_CLAWPATCH46_FORT3_COPY_FACE)

void FCLAW2D_CLAWPATCH46_FORT3_COPY_FACE(const int* mx, const int* my, const int* mz,
                                        const int* mbc, const int* meqn,
                                        double qthis[],double qneighbor[], 
                                        const int* a_idir,
                                        struct fclaw2d_patch_transform_data** 
                                        transform_cptr);


#define FCLAW2D_CLAWPATCH46_FORT3_AVERAGE_FACE \
             FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort3_average_face, \
                            FCLAW2D_CLAWPATCH46_FORT3_AVERAGE_FACE)
void FCLAW2D_CLAWPATCH46_FORT3_AVERAGE_FACE(const int* mx, const int* my, const int* mz,
                                           const int* mbc, const int* meqn,
                                           double qcoarse[],double qfine[],
                                           double areacoarse[], double areafine[],
                                           const int* idir, const int* iside,
                                           const int* num_neighbors,
                                           const int* refratio, const int* igrid,
                                           const int* manifold, 
                                           struct fclaw2d_patch_transform_data** 
                                           transform_cptr);
  
#define FCLAW2D_CLAWPATCH46_FORT3_INTERPOLATE_FACE \
              FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort3_interpolate_face, \
                             FCLAW2D_CLAWPATCH46_FORT3_INTERPOLATE_FACE)
void FCLAW2D_CLAWPATCH46_FORT3_INTERPOLATE_FACE(const int* mx, const int* my, const int* mz,
                                               const int* mbc,const int* meqn,
                                               double qcoarse[],double qfine[],
                                               const int* idir, const int* iside,
                                               const int* num_neighbors,
                                               const int* refratio, const int* igrid,
                                               struct fclaw2d_patch_transform_data** 
                                               transform_cptr);

#define FCLAW2D_CLAWPATCH46_FORT3_COPY_CORNER \
             FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort3_copy_corner, \
                            FCLAW2D_CLAWPATCH46_FORT3_COPY_CORNER)
void FCLAW2D_CLAWPATCH46_FORT3_COPY_CORNER(const int* mx, const int* my, const int* mz,
                                          const int* mbc,
                                          const int* meqn, double this_q[],
                                          double neighbor_q[],const int* a_corner,
                                          struct fclaw2d_patch_transform_data** 
                                          transform_cptr);

#define FCLAW2D_CLAWPATCH46_FORT3_AVERAGE_CORNER \
      FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort3_average_corner, \
                     FCLAW2D_CLAWPATCH46_FORT3_AVERAGE_CORNER)
void FCLAW2D_CLAWPATCH46_FORT3_AVERAGE_CORNER(const int* mx, const int* my, const int* mz,
                                             const int* mbc,
                                             const int* meqn, const int* a_refratio,
                                             double qcoarse[], double qfine[],
                                             double areacoarse[], double areafine[],
                                             const int* manifold,const int* a_corner, 
                                             struct fclaw2d_patch_transform_data** 
                                             transform_cptr);
  
#define FCLAW2D_CLAWPATCH46_FORT3_INTERPOLATE_CORNER \
      FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort3_interpolate_corner, \
                     FCLAW2D_CLAWPATCH46_FORT3_INTERPOLATE_CORNER)
void FCLAW2D_CLAWPATCH46_FORT3_INTERPOLATE_CORNER(const int* mx, const int* my, const int* mz,
                                                 const int* mbc,const int* meqn, 
                                                 const int* a_refratio, double this_q[],
                                                 double neighbor_q[], const int* a_corner,
                                                 struct fclaw2d_patch_transform_data** 
                                                 transform_cptr);



/* ------------------------- Pillow grid block corner routines  ----------------------- */

#if 0
/* Doesn't yet work in 3d */
#define FCLAW2D_PILLOW46_COPY_BLOCK_CORNER \
               FCLAW_F77_FUNC(fclaw2d_pillow46_copy_block_corner, \
                              FCLAW2D_PILLOW46_COPY_BLOCK_CORNER)

void FCLAW2D_PILLOW46_COPY_BLOCK_CORNER(const int* mx, const int* my,const int *mz,
                                        const int* mbc, const int* meqn,
                                        double qthis[], 
                                        double qneighbor[], 
                                        const int* icorner,
                                        const int* iblock);

#define FCLAW2D_PILLOW46_AVERAGE_BLOCK_CORNER \
          FCLAW_F77_FUNC(fclaw2d_pillow46_average_block_corner,\
                         FCLAW2D_PILLOW46_AVERAGE_BLOCK_CORNER)

void  FCLAW2D_PILLOW46_AVERAGE_BLOCK_CORNER(const int* mx, const int* my, const int* mz, 
                                            const int* mbc,
                                            const int* meqn, 
                                            const int* refratio, 
                                            double qcoarse[],
                                            double qfine[], 
                                            double areacoarse[], 
                                            double areafine[],
                                            const int* a_coarse_corner,
                                            const int* blockno);

// Averaging at block boundaries between coarse and fine grids.
#define FCLAW2D_PILLOW46_INTERPOLATE_BLOCK_CORNER \
          FCLAW_F77_FUNC(fclaw2d_pillow46_interpolate_block_corner, \
                         FCLAW2D_PILLOW46_INTERPOLATE_BLOCK_CORNER)

void  FCLAW2D_PILLOW46_INTERPOLATE_BLOCK_CORNER(const int* mx, const int* my, const int* mz,
                                                int* mbc,
                                                int* meqn, int* refratio,
                                                double qcoarse[],
                                                double qfine[], 
                                                int* icoarse_corner,
                                                int* blockno);

#endif

/* ------------------------------------ Output functions ------------------------------ */

#define  FCLAW2D_CLAWPATCH46_FORT3_OUTPUT_ASCII \
           FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort3_output_ascii, \
                          FCLAW2D_CLAWPATCH46_FORT3_OUTPUT_ASCII)
void FCLAW2D_CLAWPATCH46_FORT3_OUTPUT_ASCII(const char *matname1,
                                            const int* mx,const int* my, const int* mz,
                                            const int* meqn, const int* mbc,
                                            const double* xlower, const double* ylower,
                                            const double* zlower,
                                            const double* dx,     const double* dy,
                                            const double* dz,
                                            double q[],
                                            const int* patch_num, const int* level,
                                            const int* blockno,   const int* mpirank);

#define FCLAW2D_CLAWPATCH46_FORT3_HEADER_ASCII \
         FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort3_header_ascii, \
                        FCLAW2D_CLAWPATCH46_FORT3_HEADER_ASCII)
void FCLAW2D_CLAWPATCH46_FORT3_HEADER_ASCII(const char* matname1, const char* matname2,
                                           const double* time, const int* meqn, 
                                           const int* maux, 
                                           const int* ngrids);
/* ----------------------------- Diagnostics functions -------------------------------- */

#define FCLAW2D_CLAWPATCH46_FORT3_CONSERVATION_CHECK \
          FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort3_conservation_check, \
                         FCLAW2D_CLAWPATCH46_FORT3_CONSERVATION_CHECK)
void FCLAW2D_CLAWPATCH46_FORT3_CONSERVATION_CHECK(const int *mx, const int *my, 
                                                 const int *mz, 
                                                 const int* mbc, const int* meqn,
                                                 const double *dx, const double *dy,
                                                 const double *dz,
                                                 double* area, double *q, double* sum,
                                                 double* c_kahan);

#define FCLAW2D_CLAWPATCH46_FORT3_COMPUTE_PATCH_AREA \
          FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort3_compute_patch_area, \
                         FCLAW2D_CLAWPATCH46_FORT3_COMPUTE_PATCH_AREA)

double FCLAW2D_CLAWPATCH46_FORT3_COMPUTE_PATCH_AREA(const int *mx, const int* my, 
                                                   const int* mz, const int *mbc, 
                                                   const double* dx, const double *dy,
                                                   double* dz, double area[]);

#define FCLAW2D_CLAWPATCH46_FORT3_COMPUTE_ERROR_NORM \
           FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort3_compute_error_norm, \
                          FCLAW2D_CLAWPATCH46_FORT3_COMPUTE_ERROR_NORM)

void FCLAW2D_CLAWPATCH46_FORT3_COMPUTE_ERROR_NORM (const int* blockno, 
                                                  const int* mx, 
                                                  const int* my,
                                                  const int* mz,
                                                  const int* mbc,
                                                  const int* meqn, 
                                                  const double* dx,
                                                  const double* dy,
                                                  const double* dz,
                                                  double area[], double error[],
                                                  double* error_norm);

/* ----------------------------- Parallel ghost patches  ------------------------------ */


#define FCLAW2D_CLAWPATCH46_FORT3_LOCAL_GHOST_PACK \
          FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort3_local_ghost_pack, \
                         FCLAW2D_CLAWPATCH46_FORT3_LOCAL_GHOST_PACK)
void FCLAW2D_CLAWPATCH46_FORT3_LOCAL_GHOST_PACK(const int *mx, 
                                                const int *my, 
                                                const int *mz, 
                                                const int *mbc,
                                                const int *meqn, 
                                                const int *mint,
                                                double qdata[], double area[],
                                                double qpack[], int *psize,
                                                const int *packmode, int *ierror);



/* ----------------------------- User convenience headers ----------------------------- */
#define CLAWPATCH46_TAG4REFINEMENT FCLAW_F77_FUNC(clawpatch46_tag4refinement, \
                                                 CLAWPATCH46_TAG4REFINEMENT)

void CLAWPATCH46_TAG4REFINEMENT(const int* mx,const int* my, const int* mz,
                               const int* mbc,const int* meqn,
                               const double* xlower, const double* ylower,
                               const double* zlower, 
                               const double* dx, const double* dy, const double *dz,
                               const int* blockno,
                               double q[],
                               const double* tag_threshold,
                               const int* init_flag,
                               int* tag_patch);



#define CLAWPATCH46_TAG4COARSENING FCLAW_F77_FUNC(clawpatch46_tag4coarsening, \
                                                CLAWPATCH46_TAG4COARSENING)

void CLAWPATCH46_TAG4COARSENING(const int* mx, const int* my, const int *mz,
                               const int* mbc, const int* meqn,
                               const double* xlower, const double* ylower,
                               const double* zlower,
                               const double* dx, const double* dy, const double *dz,
                               const int* blockno,
                               double *qfine[],
                               const double* tag_threshold,
                               const int* initflag,
                               int* tag_patch);



#ifdef __cplusplus
}
#endif

#endif
