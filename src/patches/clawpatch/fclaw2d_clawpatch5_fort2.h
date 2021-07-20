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

#ifndef FCLAW2D_CLAWPATCH5_FORT2_H
#define FCLAW2D_CLAWPATCH5_FORT2_H

#ifdef __cplusplus
extern "C"
{
#endif


#if 0
/* Fix syntax highlighting */
#endif

struct fclaw2d_patch_transform_data;  /* Should be replaced by long int?  */

/* ------------------------------ Time stepping functions ----------------------------- */
#define FCLAW2D_CLAWPATCH5_FORT_TIMEINTERP \
                     FCLAW_F77_FUNC (fclaw2d_clawpatch5_fort_timeinterp, \
                     FCLAW2D_CLAWPATCH5_FORT_TIMEINTERP)

void FCLAW2D_CLAWPATCH5_FORT_TIMEINTERP(const int *mx, const int* my, const int* mbc,
                                         const int *meqn, const int* psize,
                                         double qcurr[], double qlast[],
                                         double qinterp[],const double* alpha,
                                         const int* ierror);

/* --------------------------------- Regridding functions ----------------------------- */

#define FCLAW2D_CLAWPATCH5_FORT_TAG4REFINEMENT \
                   FCLAW_F77_FUNC(fclaw2d_clawpatch5_fort_tag4refinement, \
                   FCLAW2D_CLAWPATCH5_FORT_TAG4REFINEMENT)

void FCLAW2D_CLAWPATCH5_FORT_TAG4REFINEMENT(const int* mx,const int* my,
                                             const int* mbc,const int* meqn,
                                             const double* xlower, const double* ylower,
                                             const double* dx, const double* dy,
                                             const int* blockno,
                                             double q[],
                                             const double* tag_threshold,
                                             const int* initflag,
                                             int* tag_patch);



#define FCLAW2D_CLAWPATCH5_FORT_TAG4COARSENING \
                FCLAW_F77_FUNC(fclaw2d_clawpatch5_fort_tag4coarsening, \
                FCLAW2D_CLAWPATCH5_FORT_TAG4COARSENING)

void FCLAW2D_CLAWPATCH5_FORT_TAG4COARSENING(const int* mx, const int* my,
                                             const int* mbc, const int* meqn,
                                             double xlower[], double ylower[],
                                             const double* dx, const double* dy,
                                             const int* blockno,
                                             double q0[],double q1[],
                                             double q2[],double q3[],
                                             const double* tag_threshold,
                                             const int* initflag,
                                             int* tag_patch);

#define FCLAW2D_CLAWPATCH5_FORT_INTERPOLATE2FINE \
                   FCLAW_F77_FUNC(fclaw2d_clawpatch5_fort_interpolate2fine, \
                   FCLAW2D_CLAWPATCH5_FORT_INTERPOLATE2FINE)

void FCLAW2D_CLAWPATCH5_FORT_INTERPOLATE2FINE(const int* mx,const int* my,
                                               const int* mbc, const int* meqn,
                                               double qcoarse[], double qfine[],
                                               double areacoarse[], double areafine[],
                                               const int* igrid, const int* manifold);
  
#define FCLAW2D_CLAWPATCH5_FORT_AVERAGE2COARSE \
                   FCLAW_F77_FUNC(fclaw2d_clawpatch5_fort_average2coarse, \
                   FCLAW2D_CLAWPATCH5_FORT_AVERAGE2COARSE)

void FCLAW2D_CLAWPATCH5_FORT_AVERAGE2COARSE(const int* mx, const int* my,
                                             const int* mbc, const int* meqn,
                                             double qcoarse[],double qfine[],
                                             double areacoarse[],double areafine[],
                                             const int* igrid, const int* manifold);



/* ---------------------------------- Ghost filling  ---------------------------------- */

#define FCLAW2D_CLAWPATCH5_FORT_COPY_FACE \
                   FCLAW_F77_FUNC(fclaw2d_clawpatch5_fort_copy_face, \
                   FCLAW2D_CLAWPATCH5_FORT_COPY_FACE)

void FCLAW2D_CLAWPATCH5_FORT_COPY_FACE(const int* mx, const int* my, const int* mbc, 
                                       const int* meqn,
                                        double qthis[],double qneighbor[], const int* a_idir,
                                        struct fclaw2d_patch_transform_data** 
                                        transform_cptr);


#define FCLAW2D_CLAWPATCH5_FORT_AVERAGE_FACE \
                    FCLAW_F77_FUNC(fclaw2d_clawpatch5_fort_average_face, \
                    FCLAW2D_CLAWPATCH5_FORT_AVERAGE_FACE)

void FCLAW2D_CLAWPATCH5_FORT_AVERAGE_FACE(const int* mx, const int* my, const int* mbc,
                                           const int* meqn,
                                           double qcoarse[],double qfine[],
                                           double areacoarse[], double areafine[],
                                           const int* idir, const int* iside,
                                           const int* num_neighbors,
                                           const int* refratio, const int* igrid,
                                           const int* manifold, 
                                           struct fclaw2d_patch_transform_data** 
                                           transform_cptr);
  
#define FCLAW2D_CLAWPATCH5_FORT_INTERPOLATE_FACE \
                    FCLAW_F77_FUNC(fclaw2d_clawpatch5_fort_interpolate_face, \
                    FCLAW2D_CLAWPATCH5_FORT_INTERPOLATE_FACE)

void FCLAW2D_CLAWPATCH5_FORT_INTERPOLATE_FACE(const int* mx,const int* my,const int* mbc,
                                              const int* meqn,
                                              double qcoarse[],double qfine[],
                                              const int* idir, const int* iside,
                                               const int* num_neighbors,
                                               const int* refratio, const int* igrid,
                                               struct fclaw2d_patch_transform_data** 
                                               transform_cptr);

#define FCLAW2D_CLAWPATCH5_FORT_COPY_CORNER \
                          FCLAW_F77_FUNC(fclaw2d_clawpatch5_fort_copy_corner, \
                          FCLAW2D_CLAWPATCH5_FORT_COPY_CORNER)

void FCLAW2D_CLAWPATCH5_FORT_COPY_CORNER(const int* mx, const int* my, const int* mbc,
                                         const int* meqn, double this_q[], 
                                         double neighbor_q[],
                                         const int* a_corner,
                                         struct fclaw2d_patch_transform_data** 
                                         transform_cptr);

#define FCLAW2D_CLAWPATCH5_FORT_AVERAGE_CORNER \
                         FCLAW_F77_FUNC(fclaw2d_clawpatch5_fort_average_corner, \
                         FCLAW2D_CLAWPATCH5_FORT_AVERAGE_CORNER)

void FCLAW2D_CLAWPATCH5_FORT_AVERAGE_CORNER(const int* mx, const int* my, const int* mbc,
                                            const int* meqn, const int* a_refratio,
                                            double qcoarse[], double qfine[],
                                            double areacoarse[], double areafine[],
                                            const int* manifold,
                                            const int* a_corner, 
                                            struct fclaw2d_patch_transform_data** 
                                            transform_cptr);
  
#define FCLAW2D_CLAWPATCH5_FORT_INTERPOLATE_CORNER \
                    FCLAW_F77_FUNC(fclaw2d_clawpatch5_fort_interpolate_corner, \
                    FCLAW2D_CLAWPATCH5_FORT_INTERPOLATE_CORNER)

void FCLAW2D_CLAWPATCH5_FORT_INTERPOLATE_CORNER(const int* mx, const int* my, 
                                                const int* mbc,
                                                const int* meqn, const int* a_refratio, 
                                                double this_q[],
                                                double neighbor_q[], const int* a_corner,
                                                struct fclaw2d_patch_transform_data** 
                                                transform_cptr);



/* ------------------------- Pillow grid block corner routines  ----------------------- */

#define FCLAW2D_PILLOW5_COPY_BLOCK_CORNER \
               FCLAW_F77_FUNC(fclaw2d_pillow5_copy_block_corner, \
                              FCLAW2D_PILLOW5_COPY_BLOCK_CORNER)

void FCLAW2D_PILLOW5_COPY_BLOCK_CORNER(int* mx, int* my,
                                        int* mbc, int* meqn,
                                        double qthis[], 
                                        double qneighbor[], 
                                        int* icorner,
                                        int* iblock);

#define FCLAW2D_PILLOW5_AVERAGE_BLOCK_CORNER \
          FCLAW_F77_FUNC(fclaw2d_pillow5_average_block_corner,\
                         FCLAW2D_PILLOW5_AVERAGE_BLOCK_CORNER)

void  FCLAW2D_PILLOW5_AVERAGE_BLOCK_CORNER(int* mx, int* my, int* mbc,
                                            int* meqn, 
                                            int* refratio, 
                                            double qcoarse[],
                                            double qfine[], 
                                            double areacoarse[], 
                                            double areafine[],
                                            int* a_coarse_corner,
                                            int* blockno);

// Averaging at block boundaries between coarse and fine grids.
#define FCLAW2D_PILLOW5_INTERPOLATE_BLOCK_CORNER \
          FCLAW_F77_FUNC(fclaw2d_pillow5_interpolate_block_corner, \
                         FCLAW2D_PILLOW5_INTERPOLATE_BLOCK_CORNER)

void  FCLAW2D_PILLOW5_INTERPOLATE_BLOCK_CORNER(int* mx, int* my, int* mbc,
                                                int* meqn, int* refratio,
                                                double qcoarse[],
                                                double qfine[], 
                                                int* icoarse_corner,
                                                int* blockno);


/* ------------------------------------ Output functions ------------------------------ */

#define  FCLAW2D_CLAWPATCH5_FORT_OUTPUT_ASCII \
                       FCLAW_F77_FUNC(fclaw2d_clawpatch5_fort_output_ascii, \
                       FCLAW2D_CLAWPATCH5_FORT_OUTPUT_ASCII)

void FCLAW2D_CLAWPATCH5_FORT_OUTPUT_ASCII(char* matname1,
                                          int* mx,        int* my,
                                          int* meqn,      int* mbc,
                                          double* xlower, double* ylower,
                                          double* dx,     double* dy,
                                          double q[],
                                          int* patch_num, int* level,
                                          int* blockno,   int* mpirank);

#define FCLAW2D_CLAWPATCH5_FORT_HEADER_ASCII \
                   FCLAW_F77_FUNC(fclaw2d_clawpatch5_fort_header_ascii, \
                   FCLAW2D_CLAWPATCH5_FORT_HEADER_ASCII)

void FCLAW2D_CLAWPATCH5_FORT_HEADER_ASCII(const char* matname1, const char* matname2,
                                          const double* time, const int* meqn, 
                                          const int* maux, const int* ngrids);

/* ----------------------------- Diagnostics functions -------------------------------- */

#define FCLAW2D_CLAWPATCH5_FORT_CONSERVATION_CHECK \
               FCLAW_F77_FUNC(fclaw2d_clawpatch5_fort_conservation_check, \
               FCLAW2D_CLAWPATCH5_FORT_CONSERVATION_CHECK)

void FCLAW2D_CLAWPATCH5_FORT_CONSERVATION_CHECK(int *mx, int *my, int* mbc, int* meqn,
                                                double *dx, double *dy,
                                                double* area, double *q, double* sum,
                                                double *c_kahan);

#define FCLAW2D_CLAWPATCH5_FORT_COMPUTE_PATCH_AREA \
                      FCLAW_F77_FUNC(fclaw2d_clawpatch5_fort_compute_patch_area, \
                      FCLAW2D_CLAWPATCH5_FORT_COMPUTE_PATCH_AREA)

double FCLAW2D_CLAWPATCH5_FORT_COMPUTE_PATCH_AREA(int *mx, int* my, int*mbc, double* dx,
                                                  double* dy, double area[]);

#define FCLAW2D_CLAWPATCH5_FORT_COMPUTE_ERROR_NORM \
                   FCLAW_F77_FUNC(fclaw2d_clawpatch5_fort_compute_error_norm, \
                   FCLAW2D_CLAWPATCH5_FORT_COMPUTE_ERROR_NORM)

void FCLAW2D_CLAWPATCH5_FORT_COMPUTE_ERROR_NORM (int* blockno, int* mx,int* my,int* mbc,
                                                 int* meqn,double* dx,double* dy,
                                                 double area[], double error[],
                                                 double* error_norm);

/* ----------------------------- Parallel ghost patches  ------------------------------ */


#define FCLAW2D_CLAWPATCH5_FORT_LOCAL_GHOST_PACK \
                      FCLAW_F77_FUNC(fclaw2d_clawpatch5_fort_local_ghost_pack, \
                      FCLAW2D_CLAWPATCH5_FORT_LOCAL_GHOST_PACK)

void FCLAW2D_CLAWPATCH5_FORT_LOCAL_GHOST_PACK(int *mx, int *my, int *mbc,
                                              int *meqn, int *mint,
                                              double qdata[], double area[],
                                              double qpack[], int *psize,
                                              int *packmode, int *ierror);


/* ----------------------------- User convenience functions  -------------------------- */

#define CLAWPATCH5_TAG4REFINEMENT FCLAW_F77_FUNC(clawpatch5_tag4refinement, \
                                                 CLAWPATCH5_TAG4REFINEMENT)

void CLAWPATCH5_TAG4REFINEMENT(const int* mx,const int* my,
                              const int* mbc,const int* meqn,
                              const double* xlower, const double* ylower,
                              const double* dx, const double* dy,
                              const int* blockno,
                              double q[],
                              const double* tag_threshold,
                              const int* init_flag,
                              int* tag_patch);


#define CLAWPATCH5_TAG4COARSENING FCLAW_F77_FUNC(clawpatch5_tag4coarsening, \
                                                 CLAWPATCH5_TAG4COARSENING)

void CLAWPATCH5_TAG4COARSENING(const int* mx, const int* my,
                              const int* mbc, const int* meqn,
                              const double* xlower, const double* ylower,
                              const double* dx, const double* dy,
                              const int* blockno,
                              double q0[],double q1[],
                              double q2[],double q3[],
                              const double* tag_threshold,
                              const int* initflag,
                              int* tag_patch);


#ifdef __cplusplus
}
#endif

#endif
