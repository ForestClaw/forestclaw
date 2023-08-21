/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#ifndef FCLAW3DX_CLAWPATCH46_FORT_H
#define FCLAW3DX_CLAWPATCH46_FORT_H

/**
 * @file
 * C declarations of clawpatch 4.6 fortran subroutines
 */
#ifdef __cplusplus
extern "C"
{
#endif


#if 0
/* Fix syntax highlighting */
#endif

struct fclaw_patch_transform_data;  /* Should be replaced by long int?  */

/** @{ @name Time Stepping Functions *//*---------------------------------------------- */

/** Fortran subroutine name */
#define FCLAW3D_CLAWPATCH46_FORT_TIMEINTERP \
            FCLAW_F77_FUNC (fclaw3d_clawpatch46_fort_timeinterp, \
                            FCLAW3D_CLAWPATCH46_FORT_TIMEINTERP)

/** @copydoc fclaw3d_clawpatch46_fort_timeinterp() */
void FCLAW3D_CLAWPATCH46_FORT_TIMEINTERP(const int *mx, 
                                         const int* my, 
                                         const int* mz, 
                                         const int* mbc,
                                         const int *meqn, 
                                         const int* psize,
                                         double qcurr[], 
                                         double qlast[],
                                         double qinterp[],
                                         const double* alpha,
                                         const int* ierror);

/** @} */

/** @{ @name Regridding Functions *//*------------------------------------------------- */

/** Fortran subroutine name */
#define FCLAW3D_CLAWPATCH46_FORT_TAG4REFINEMENT \
           FCLAW_F77_FUNC(fclaw3d_clawpatch46_fort_tag4refinement, \
                          FCLAW3D_CLAWPATCH46_FORT_TAG4REFINEMENT)

/** @copydoc fclaw3d_clawpatch46_fort_tag4refinement() */
void FCLAW3D_CLAWPATCH46_FORT_TAG4REFINEMENT(const int* mx,
                                             const int* my, 
                                             const int* mz,
                                             const int* mbc,
                                             const int* meqn,
                                             const double* xlower, 
                                             const double* ylower, 
                                             const double* zlower,
                                             const double* dx, 
                                             const double* dy, 
                                             const double* dz,
                                             const int* blockno,
                                             double q[],
                                             const double* tag_threshold,
                                             const int* init_flag,
                                             int* tag_patch);



/** Fortran subroutine name */
#define FCLAW3DX_CLAWPATCH46_FORT_TAG4COARSENING \
           FCLAW_F77_FUNC(fclaw3dx_clawpatch46_fort_tag4coarsening, \
                          FCLAW3DX_CLAWPATCH46_FORT_TAG4COARSENING)

/** @copydoc fclaw3d_clawpatch46_fort_tag4coarsening() */
void FCLAW3DX_CLAWPATCH46_FORT_TAG4COARSENING(const int* mx, 
                                             const int* my, 
                                             const int* mz,
                                             const int* mbc, 
                                             const int* meqn,
                                             double xlower[], 
                                             double ylower[], 
                                             double zlower[],
                                             const double* dx, 
                                             const double* dy, 
                                             const double* dz,
                                             const int* blockno,
                                             double q0[],
                                             double q1[],
                                             double q2[],
                                             double q3[],
                                             const double* tag_threshold,
                                             const int* initflag,
                                             int* tag_patch);

/** Fortran subroutine name */
#define FCLAW3D_CLAWPATCH46_FORT_TAG4COARSENING \
           FCLAW_F77_FUNC(fclaw3d_clawpatch46_fort_tag4coarsening, \
                          FCLAW3D_CLAWPATCH46_FORT_TAG4COARSENING)

/** @copydoc fclaw3d_clawpatch46_fort_tag4coarsening() */
void FCLAW3D_CLAWPATCH46_FORT_TAG4COARSENING(const int* mx, 
                                             const int* my, 
                                             const int* mz,
                                             const int* mbc, 
                                             const int* meqn,
                                             double xlower[], 
                                             double ylower[], 
                                             double zlower[],
                                             const double* dx, 
                                             const double* dy, 
                                             const double* dz,
                                             const int* blockno,
                                             double q0[],
                                             double q1[],
                                             double q2[],
                                             double q3[],
                                             double q4[],
                                             double q5[],
                                             double q6[],
                                             double q7[],
                                             const double* tag_threshold,
                                             const int* initflag,
                                             int* tag_patch);


/** Fortran subroutine name */
#define FCLAW3D_CLAWPATCH46_FORT_INTERPOLATE2FINE \
           FCLAW_F77_FUNC(fclaw3d_clawpatch46_fort_interpolate2fine, \
                          FCLAW3D_CLAWPATCH46_FORT_INTERPOLATE2FINE)

/** @copydoc fclaw3d_clawpatch46_fort_interpolate2fine() */
void FCLAW3D_CLAWPATCH46_FORT_INTERPOLATE2FINE(const int* mx,
                                               const int* my, 
                                               const int* mz,
                                               const int* mbc, 
                                               const int* meqn,
                                               double qcoarse[], 
                                               double qfine[],
                                               double areacoarse[], 
                                               double areafine[],
                                               const int* igrid, 
                                               const int* manifold);
 
/** Fortran subroutine name */
#define FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE2FINE \
           FCLAW_F77_FUNC(fclaw3dx_clawpatch46_fort_interpolate2fine, \
                          FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE2FINE)

/** @copydoc fclaw3dx_clawpatch46_fort_interpolate2fine() */
void FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE2FINE(const int* mx,
                                                const int* my, 
                                                const int* mz,
                                                const int* mbc, 
                                                const int* meqn,
                                                double qcoarse[], 
                                                double qfine[],
                                                double areacoarse[], 
                                                double areafine[],
                                                const int* igrid, 
                                                const int* manifold);

 /** Fortran subroutine name */
#define FCLAW3D_CLAWPATCH46_FORT_AVERAGE2COARSE \
      FCLAW_F77_FUNC(fclaw3d_clawpatch46_fort_average2coarse, \
                     FCLAW3D_CLAWPATCH46_FORT_AVERAGE2COARSE)

/** @copydoc fclaw3d_clawpatch46_fort_average2coarse() */
void FCLAW3D_CLAWPATCH46_FORT_AVERAGE2COARSE(const int* mx, 
                                             const int* my, 
                                             const int* mz,
                                             const int* mbc, 
                                             const int* meqn,
                                             double qcoarse[],
                                             double qfine[],
                                             double areacoarse[],
                                             double areafine[],
                                             const int* igrid, 
                                             const int* manifold);

 
/** Fortran subroutine name */
#define FCLAW3DX_CLAWPATCH46_FORT_AVERAGE2COARSE \
      FCLAW_F77_FUNC(fclaw3dx_clawpatch46_fort_average2coarse, \
                     FCLAW3DX_CLAWPATCH46_FORT_AVERAGE2COARSE)

/** @copydoc fclaw3dx_clawpatch46_fort_average2coarse() */
void FCLAW3DX_CLAWPATCH46_FORT_AVERAGE2COARSE(const int* mx, 
                                              const int* my, 
                                              const int* mz,
                                              const int* mbc, 
                                              const int* meqn,
                                              double qcoarse[],
                                              double qfine[],
                                              double areacoarse[],
                                              double areafine[],
                                              const int* igrid, 
                                              const int* manifold);



/** @} */

/** @{ @name Ghost Filling *//*-------------------------------------------------------- */

/** Fortran subroutine name */
#define FCLAW3D_CLAWPATCH46_FORT_COPY_FACE \
         FCLAW_F77_FUNC(fclaw3d_clawpatch46_fort_copy_face, \
                        FCLAW3D_CLAWPATCH46_FORT_COPY_FACE)

/** @copydoc fclaw3d_clawpatch46_fort_copy_face() */
void FCLAW3D_CLAWPATCH46_FORT_COPY_FACE(const int* mx, 
                                        const int* my, 
                                        const int* mz, 
                                        const int* mbc, 
                                        const int* meqn,
                                        double qthis[],
                                        double qneighbor[], 
                                        const int* a_idir,
                                        struct fclaw_patch_transform_data** 
                                        transform_cptr);

/** Fortran subroutine name */
#define FCLAW3DX_CLAWPATCH46_FORT_COPY_FACE \
         FCLAW_F77_FUNC(fclaw3dx_clawpatch46_fort_copy_face, \
                        FCLAW3Dx_CLAWPATCH46_FORT_COPY_FACE)

/** @copydoc fclaw3d_clawpatch46_fort_copy_face() */
void FCLAW3DX_CLAWPATCH46_FORT_COPY_FACE(const int* mx, 
                                         const int* my, 
                                         const int* mz, 
                                         const int* mbc, 
                                         const int* meqn,
                                         double qthis[],
                                         double qneighbor[], 
                                         const int* a_idir,
                                         struct fclaw_patch_transform_data** 
                                         transform_cptr);


/** Fortran subroutine name */
#define FCLAW3D_CLAWPATCH46_FORT_AVERAGE_FACE \
             FCLAW_F77_FUNC(fclaw3d_clawpatch46_fort_average_face, \
                            FCLAW3D_CLAWPATCH46_FORT_AVERAGE_FACE)
/** @copydoc fclaw3d_clawpatch46_fort_average_face() */
void FCLAW3D_CLAWPATCH46_FORT_AVERAGE_FACE(const int* mx, 
                                           const int* my, 
                                           const int* mz, 
                                           const int* mbc,
                                           const int* meqn,
                                           double qcoarse[],
                                           double qfine[],
                                           double areacoarse[], 
                                           double areafine[],
                                           const int* idir, 
                                           const int* iside,
                                           const int* num_neighbors,
                                           const int* refratio, 
                                           const int* igrid,
                                           const int* manifold, 
                                           struct fclaw_patch_transform_data** 
                                           transform_cptr);
 
/** Fortran subroutine name */
#define FCLAW3DX_CLAWPATCH46_FORT_AVERAGE_FACE \
             FCLAW_F77_FUNC(fclaw3dx_clawpatch46_fort_average_face, \
                            FCLAW3DX_CLAWPATCH46_FORT_AVERAGE_FACE)
/** @copydoc fclaw3dx_clawpatch46_fort_average_face() */
void FCLAW3DX_CLAWPATCH46_FORT_AVERAGE_FACE(const int* mx, 
                                            const int* my, 
                                            const int* mz, 
                                            const int* mbc,
                                            const int* meqn,
                                            double qcoarse[],
                                            double qfine[],
                                            double areacoarse[], 
                                            double areafine[],
                                            const int* idir, 
                                            const int* iside,
                                            const int* num_neighbors,
                                            const int* refratio, 
                                            const int* igrid,
                                            const int* manifold, 
                                            struct fclaw_patch_transform_data** 
                                            transform_cptr);
 /** Fortran subroutine name */
#define FCLAW3D_CLAWPATCH46_FORT_INTERPOLATE_FACE \
              FCLAW_F77_FUNC(fclaw3d_clawpatch46_fort_interpolate_face, \
                             FCLAW3D_CLAWPATCH46_FORT_INTERPOLATE_FACE)
/** @copydoc fclaw3d_clawpatch46_fort_interpolate_face() */
void FCLAW3D_CLAWPATCH46_FORT_INTERPOLATE_FACE(const int* mx, 
                                                const int* my, 
                                                const int* mz,
                                                const int* mbc,
                                                const int* meqn,
                                                double qcoarse[],
                                                double qfine[],
                                                const int* idir, 
                                                const int* iside,
                                                const int* num_neighbors,
                                                const int* refratio, 
                                                const int* igrid,
                                                struct fclaw_patch_transform_data** 
                                                transform_cptr);

 
/** Fortran subroutine name */
#define FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE_FACE \
              FCLAW_F77_FUNC(fclaw3dx_clawpatch46_fort_interpolate_face, \
                             FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE_FACE)
/** @copydoc fclaw3dx_clawpatch46_fort_interpolate_face() */
void FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE_FACE(const int* mx, 
                                                const int* my, 
                                                const int* mz,
                                                const int* mbc,
                                                const int* meqn,
                                                double qcoarse[],
                                                double qfine[],
                                                const int* idir, 
                                                const int* iside,
                                                const int* num_neighbors,
                                                const int* refratio, 
                                                const int* igrid,
                                                struct fclaw_patch_transform_data** 
                                                transform_cptr);

/** Fortran subroutine name */
#define FCLAW3D_CLAWPATCH46_FORT_COPY_EDGE \
         FCLAW_F77_FUNC(fclaw3d_clawpatch46_fort_copy_edge, \
                        FCLAW3D_CLAWPATCH46_FORT_COPY_EDGE)

/** @copydoc fclaw3d_clawpatch46_fort_copy_face() */
void FCLAW3D_CLAWPATCH46_FORT_COPY_EDGE(const int* mx, 
                                        const int* my, 
                                        const int* mz, 
                                        const int* mbc, 
                                        const int* meqn,
                                        double qthis[],
                                        double qneighbor[], 
                                        const int* iedge,
                                        struct fclaw_patch_transform_data** 
                                        transform_cptr);

/** Fortran subroutine name */
#define FCLAW3D_CLAWPATCH46_FORT_AVERAGE_EDGE \
      FCLAW_F77_FUNC(fclaw3d_clawpatch46_fort_average_edge, \
                     FCLAW3D_CLAWPATCH46_FORT_AVERAGE_EDGE)
/** @copydoc fclaw3d_clawpatch46_fort_average_corner() */
void FCLAW3D_CLAWPATCH46_FORT_AVERAGE_EDGE(const int* mx, 
                                           const int* my, 
                                           const int* mz, 
                                           const int* mbc,
                                           const int* meqn, 
                                           const int* a_refratio,
                                           double qcoarse[], 
                                           double qfine[],
                                           double areacoarse[], 
                                           double areafine[],
                                           const int* manifold,
                                           const int* a_corner, 
                                           struct fclaw_patch_transform_data** 
                                           transform_cptr);
 /** Fortran subroutine name */
#define FCLAW3D_CLAWPATCH46_FORT_INTERPOLATE_EDGE \
      FCLAW_F77_FUNC(fclaw3d_clawpatch46_fort_interpolate_edge, \
                     FCLAW3D_CLAWPATCH46_FORT_INTERPOLATE_EDGE)
/** @copydoc fclaw3d_clawpatch46_fort_interpolate_edge() */
void FCLAW3D_CLAWPATCH46_FORT_INTERPOLATE_EDGE(const int* mx, 
                                               const int* my, 
                                               const int* mz,
                                               const int* mbc,
                                               const int* meqn, 
                                               const int* a_refratio, 
                                               double this_q[],
                                               double neighbor_q[], 
                                               const int* a_edge,
                                               struct fclaw_patch_transform_data** 
                                               transform_cptr);


/** Fortran subroutine name */
#define FCLAW3D_CLAWPATCH46_FORT_COPY_CORNER \
             FCLAW_F77_FUNC(fclaw3d_clawpatch46_fort_copy_corner, \
                            FCLAW3D_CLAWPATCH46_FORT_COPY_CORNER)
/** @copydoc fclaw3d_clawpatch46_fort_copy_corner() */
void FCLAW3D_CLAWPATCH46_FORT_COPY_CORNER(const int* mx, 
                                          const int* my, 
                                          const int* mz, 
                                          const int* mbc,
                                          const int* meqn, 
                                          double this_q[],
                                          double neighbor_q[],
                                          const int* a_corner,
                                          struct fclaw_patch_transform_data** 
                                          transform_cptr);


/** Fortran subroutine name */
#define FCLAW3DX_CLAWPATCH46_FORT_COPY_CORNER \
             FCLAW_F77_FUNC(fclaw3dx_clawpatch46_fort_copy_corner, \
                            FCLAW3DX_CLAWPATCH46_FORT_COPY_CORNER)
/** @copydoc fclaw3dx_clawpatch46_fort_copy_corner() */
void FCLAW3DX_CLAWPATCH46_FORT_COPY_CORNER(const int* mx, 
                                           const int* my, 
                                           const int* mz, 
                                           const int* mbc,
                                           const int* meqn, 
                                           double this_q[],
                                           double neighbor_q[],
                                           const int* a_corner,
                                           struct fclaw_patch_transform_data** 
                                           transform_cptr);

/** Fortran subroutine name */
#define FCLAW3D_CLAWPATCH46_FORT_AVERAGE_CORNER \
      FCLAW_F77_FUNC(fclaw3d_clawpatch46_fort_average_corner, \
                     FCLAW3D_CLAWPATCH46_FORT_AVERAGE_CORNER)
/** @copydoc fclaw3d_clawpatch46_fort_average_corner() */
void FCLAW3D_CLAWPATCH46_FORT_AVERAGE_CORNER(const int* mx, 
                                             const int* my, 
                                             const int* mz, 
                                             const int* mbc,
                                             const int* meqn, 
                                             const int* a_refratio,
                                             double qcoarse[], 
                                             double qfine[],
                                             double areacoarse[], 
                                             double areafine[],
                                             const int* manifold,
                                             const int* a_corner, 
                                             struct fclaw_patch_transform_data** 
                                             transform_cptr);
 
/** Fortran subroutine name */
#define FCLAW3DX_CLAWPATCH46_FORT_AVERAGE_CORNER \
      FCLAW_F77_FUNC(fclaw3dx_clawpatch46_fort_average_corner, \
                     FCLAW3DX_CLAWPATCH46_FORT_AVERAGE_CORNER)
/** @copydoc fclaw3dx_clawpatch46_fort_average_corner() */
void FCLAW3DX_CLAWPATCH46_FORT_AVERAGE_CORNER(const int* mx, 
                                              const int* my, 
                                              const int* mz, 
                                              const int* mbc,
                                              const int* meqn, 
                                              const int* a_refratio,
                                              double qcoarse[], 
                                              double qfine[],
                                              double areacoarse[], 
                                              double areafine[],
                                              const int* manifold,
                                              const int* a_corner, 
                                              struct fclaw_patch_transform_data** 
                                              transform_cptr);

/** Fortran subroutine name */
#define FCLAW3D_CLAWPATCH46_FORT_INTERPOLATE_CORNER \
      FCLAW_F77_FUNC(fclaw3d_clawpatch46_fort_interpolate_corner, \
                     FCLAW3D_CLAWPATCH46_FORT_INTERPOLATE_CORNER)
/** @copydoc fclaw3d_clawpatch46_fort_interpolate_corner() */
void FCLAW3D_CLAWPATCH46_FORT_INTERPOLATE_CORNER(const int* mx, 
                                                 const int* my, 
                                                 const int* mz,
                                                 const int* mbc,
                                                 const int* meqn, 
                                                 const int* a_refratio, 
                                                 double this_q[],
                                                 double neighbor_q[], 
                                                 const int* a_corner,
                                                 struct fclaw_patch_transform_data** 
                                                 transform_cptr);

 
/** Fortran subroutine name */
#define FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE_CORNER \
      FCLAW_F77_FUNC(fclaw3dx_clawpatch46_fort_interpolate_corner, \
                     FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE_CORNER)
/** @copydoc fclaw3dx_clawpatch46_fort_interpolate_corner() */
void FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE_CORNER(const int* mx, 
                                                  const int* my, 
                                                  const int* mz,
                                                  const int* mbc,
                                                  const int* meqn, 
                                                  const int* a_refratio, 
                                                  double this_q[],
                                                  double neighbor_q[], 
                                                  const int* a_corner,
                                                  struct fclaw_patch_transform_data** 
                                                  transform_cptr);

/** @} */


/* ------------------------- Pillow grid block corner routines  ----------------------- */

/** Fortran subroutine name */
#define FCLAW3DX_PILLOW46_COPY_BLOCK_CORNER \
               FCLAW_F77_FUNC(fclaw3dx_pillow46_copy_block_corner, \
                              FCLAW3DX_PILLOW46_COPY_BLOCK_CORNER)

/** @copydoc fclaw2d_pillow46_copy_block_corner() */
void FCLAW3DX_PILLOW46_COPY_BLOCK_CORNER(int* mx, 
                                         int* my, 
                                         int* mz,
                                         int* mbc, 
                                         int* meqn,
                                         double qthis[], 
                                         double qneighbor[], 
                                         int* icorner,
                                         int* iblock);

/** Fortran subroutine name */
#define FCLAW3DX_PILLOW46_AVERAGE_BLOCK_CORNER \
          FCLAW_F77_FUNC(fclaw3dx_pillow46_average_block_corner,\
                         FCLAW3DX_PILLOW46_AVERAGE_BLOCK_CORNER)

/** @copydoc fclaw2d_pillow46_average_block_corner() */
void  FCLAW3DX_PILLOW46_AVERAGE_BLOCK_CORNER(int* mx, 
                                             int* my, 
                                             int* mz, 
                                             double *dz, 
                                             int* mbc,
                                             int* meqn, 
                                             int* refratio, 
                                             double qcoarse[],
                                             double qfine[], 
                                             double areacoarse[], 
                                             double areafine[],
                                             int* coarse_corner,
                                             int* blockno);

/** Fortran subroutine name */
#define FCLAW3DX_PILLOW46_INTERPOLATE_BLOCK_CORNER \
          FCLAW_F77_FUNC(fclaw3dx_pillow46_interpolate_block_corner, \
                         FCLAW3DX_PILLOW46_INTERPOLATE_BLOCK_CORNER)

/** @copydoc fclaw2d_pillow46_interpolate_block_corner() */
void  FCLAW3DX_PILLOW46_INTERPOLATE_BLOCK_CORNER(int* mx, 
                                                 int* my, 
                                                 int* mz, 
                                                 int* mbc,
                                                 int* meqn, 
                                                 int* refratio,
                                                 double qcoarse[],
                                                 double qfine[], 
                                                 int* icoarse_corner,
                                                 int* blockno);



/* ------------------------------------ Output functions ------------------------------ */

/** Fortran subroutine name */
#define  FCLAW3D_CLAWPATCH46_FORT_OUTPUT_ASCII \
           FCLAW_F77_FUNC(fclaw3d_clawpatch46_fort_output_ascii, \
                          FCLAW3D_CLAWPATCH46_FORT_OUTPUT_ASCII)
/** @copydoc fclaw3d_clawpatch46_fort_output_ascii() */
void FCLAW3D_CLAWPATCH46_FORT_OUTPUT_ASCII(char* matname1,
                                           int* mx,        
                                           int* my, 
                                           int* mz,
                                           int* meqn,      
                                           int* mbc,
                                           double* xlower, 
                                           double* ylower, 
                                           double* zlower,
                                           double* dx,     
                                           double* dy, 
                                           double* dz,
                                           double q[],
                                           int* patch_num, 
                                           int* level,
                                           int* blockno,   
                                           int* mpirank);

/** Fortran subroutine name */
#define FCLAW3D_CLAWPATCH46_FORT_HEADER_ASCII \
         FCLAW_F77_FUNC(fclaw3d_clawpatch46_fort_header_ascii, \
                        FCLAW3D_CLAWPATCH46_FORT_HEADER_ASCII)
/** @copydoc fclaw3d_clawpatch46_fort_header_ascii() */
void FCLAW3D_CLAWPATCH46_FORT_HEADER_ASCII(const char* matname1, 
                                           const char* matname2,
                                           const double* time, 
                                           const int* meqn, 
                                           const int* maux, 
                                           const int* ngrids);

/* ----------------------------- Diagnostics functions -------------------------------- */


/** Fortran subroutine name */
#define FCLAW3D_CLAWPATCH46_FORT_CONSERVATION_CHECK \
          FCLAW_F77_FUNC(fclaw3d_clawpatch46_fort_conservation_check, \
                         FCLAW3D_CLAWPATCH46_FORT_CONSERVATION_CHECK)
/** @copydoc fclaw3d_clawpatch46_fort_conservation_check() */
void FCLAW3D_CLAWPATCH46_FORT_CONSERVATION_CHECK(int *mx, 
                                                 int *my, 
                                                 int* mz, 
                                                 int* mbc, 
                                                 int* meqn,
                                                 double *dx, 
                                                 double *dy, 
                                                 double* dz,
                                                 double* area, 
                                                 double *q, 
                                                 double* sum,
                                                 double* c_kahan);

/** Fortran subroutine name */
#define FCLAW3D_CLAWPATCH46_FORT_COMPUTE_PATCH_AREA \
          FCLAW_F77_FUNC(fclaw3d_clawpatch46_fort_compute_patch_area, \
                         FCLAW3D_CLAWPATCH46_FORT_COMPUTE_PATCH_AREA)

/** @copydoc fclaw3d_clawpatch46_fort_compute_patch_area() */
double FCLAW3D_CLAWPATCH46_FORT_COMPUTE_PATCH_AREA(int *mx, 
                                                   int* my, 
                                                   int* mz, 
                                                   int*mbc, 
                                                   double* dx, 
                                                   double* dy, 
                                                   double* dz,
                                                   double area[]);

/** Fortran subroutine name */
#define FCLAW3D_CLAWPATCH46_FORT_COMPUTE_ERROR_NORM \
           FCLAW_F77_FUNC(fclaw3d_clawpatch46_fort_compute_error_norm, \
                          FCLAW3D_CLAWPATCH46_FORT_COMPUTE_ERROR_NORM)

/** @copydoc fclaw3d_clawpatch46_fort_compute_error_norm() */
void FCLAW3D_CLAWPATCH46_FORT_COMPUTE_ERROR_NORM (int* blockno, 
                                                  int* mx,
                                                  int* my, 
                                                  int* mz,
                                                  int* mbc,
                                                  int* meqn,
                                                  double* dx, 
                                                  double* dy, 
                                                  double* dz,
                                                  double area[], 
                                                  double error[],
                                                  double* error_norm);

/* ----------------------------- Parallel ghost patches  ------------------------------ */


/** Fortran subroutine name */
#define FCLAW3D_CLAWPATCH46_FORT_LOCAL_GHOST_PACK \
          FCLAW_F77_FUNC(fclaw3d_clawpatch46_fort_local_ghost_pack, \
                         FCLAW3D_CLAWPATCH46_FORT_LOCAL_GHOST_PACK)
/** @copydoc fclaw3d_clawpatch46_fort_local_ghost_pack() */
void FCLAW3D_CLAWPATCH46_FORT_LOCAL_GHOST_PACK(const int *mx, 
                                               const int *my, 
                                               const int* mz,
                                               const int *mbc,
                                               const int *meqn, 
                                               const int *mint,
                                               double qdata[], 
                                               double area[],
                                               double qpack[], 
                                               const int *psize,
                                               const int *packmode, 
                                               int *ierror);



/* ----------------------------- User convenience headers ----------------------------- */
/** Fortran subroutine name */
#define FCLAW2DX_CLAWPATCH46_TAG4REFINEMENT \
        FCLAW_F77_FUNC(fclaw2dx_clawpatch46_tag4refinement, \
                       FCLAW2DX_CLAWPATCH46_TAG4REFINEMENT)

/** 
 * @brief @copybrief ::clawpatch_fort_tag4refinement_t 
 * 
 * For user defined clawpatch46_tag4refinment subroutine
 * 
 * @details @copydetails ::clawpatch_fort_tag4refinement_t 
 */
void FCLAW2DX_CLAWPATCH46_TAG4REFINEMENT(const int* mx,
                                         const int* my,
                                         const int* mz,
                                         const int* mbc,
                                         const int* meqn,
                                         const double* xlower, 
                                         const double* ylower,
                                         const double* zlower,
                                         const double* dx, 
                                         const double* dy,
                                         const double* dz,
                                         const int* blockno,
                                         double q[],
                                         const double* tag_threshold,
                                         const int* init_flag,
                                         int* tag_patch);



/** Fortran subroutine name */
#define FCLAW2DX_CLAWPATCH46_TAG4COARSENING \
            FCLAW_F77_FUNC(fclaw2dx_clawpatch46_tag4coarsening, \
                           FCLAW2DX_CLAWPATCH46_TAG4COARSENING)

/** 
 * @brief @copybrief ::clawpatch_fort_tag4coarsening_t 
 * 
 * For user defined clawpatch46_tag4coarsening subroutine
 * 
 * @details @copydetails ::clawpatch_fort_tag4coarsening_t 
 */
void FCLAW2DX_CLAWPATCH46_TAG4COARSENING(const int* mx, 
                                         const int* my,
                                         const int* mz,
                                         const int* mbc, 
                                         const int* meqn,
                                         const double* xlower, 
                                         const double* ylower,
                                         const double* zlower,
                                         const double* dx, 
                                         const double* dy,
                                         const double* dz,
                                         const int* blockno,
                                         double q0[],
                                         double q1[],
                                         double q2[],
                                         double q3[],
                                         const double* tag_threshold,
                                         const int* initflag,
                                         int* tag_patch);



#ifdef __cplusplus
}
#endif

#endif
