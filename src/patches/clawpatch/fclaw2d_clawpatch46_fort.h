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

#ifndef FCLAW2D_CLAWPATCH46_FORT_H
#define FCLAW2D_CLAWPATCH46_FORT_H

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

struct fclaw2d_patch_transform_data;  /* Should be replaced by long int?  */

/** @{ @name Time Stepping Functions *//*---------------------------------------------- */

/** Fortran subroutine name */
#define FCLAW2D_CLAWPATCH46_FORT_TIMEINTERP \
            FCLAW_F77_FUNC (fclaw2d_clawpatch46_fort_timeinterp, \
                            FCLAW2D_CLAWPATCH46_FORT_TIMEINTERP)

/** @copydoc fclaw2d_clawpatch46_fort_timeinterp() */
void FCLAW2D_CLAWPATCH46_FORT_TIMEINTERP(const int *mx, const int* my, const int* mbc,
                                         const int *meqn, const int* psize,
                                         double qcurr[], double qlast[],
                                         double qinterp[],const double* alpha,
                                         const int* ierror);

/** @} */

/** @{ @name Regridding Functions *//*------------------------------------------------- */

/** Fortran subroutine name */
#define FCLAW2D_CLAWPATCH46_FORT_TAG4REFINEMENT \
           FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort_tag4refinement, \
                          FCLAW2D_CLAWPATCH46_FORT_TAG4REFINEMENT)

/** @copydoc fclaw2d_clawpatch46_fort_tag4refinement() */
void FCLAW2D_CLAWPATCH46_FORT_TAG4REFINEMENT(const int* mx,const int* my,
                                             const int* mbc,const int* meqn,
                                             const double* xlower, const double* ylower,
                                             const double* dx, const double* dy,
                                             const int* blockno,
                                             double q[],
                                             const double* tag_threshold,
                                             const int* init_flag,
                                             int* tag_patch);



/** Fortran subroutine name */
#define FCLAW2D_CLAWPATCH46_FORT_TAG4COARSENING \
           FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort_tag4coarsening, \
                          FCLAW2D_CLAWPATCH46_FORT_TAG4COARSENING)

/** @copydoc fclaw2d_clawpatch46_fort_tag4coarsening() */
void FCLAW2D_CLAWPATCH46_FORT_TAG4COARSENING(const int* mx, const int* my,
                                             const int* mbc, const int* meqn,
                                             double xlower[], double ylower[],
                                             const double* dx, const double* dy,
                                             const int* blockno,
                                             double q0[],double q1[],
                                             double q2[],double q3[],
                                             const double* tag_threshold,
                                             const int* initflag,
                                             int* tag_patch);

/** Fortran subroutine name */
#define FCLAW2D_CLAWPATCH46_FORT_INTERPOLATE2FINE \
           FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort_interpolate2fine, \
                          FCLAW2D_CLAWPATCH46_FORT_INTERPOLATE2FINE)

/** @copydoc fclaw2d_clawpatch46_fort_interpolate2fine() */
void FCLAW2D_CLAWPATCH46_FORT_INTERPOLATE2FINE(const int* mx,const int* my,
                                               const int* mbc, const int* meqn,
                                               double qcoarse[], double qfine[],
                                               double areacoarse[], double areafine[],
                                               const int* igrid, const int* manifold);
  
/** Fortran subroutine name */
#define FCLAW2D_CLAWPATCH46_FORT_AVERAGE2COARSE \
      FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort_average2coarse, \
                     FCLAW2D_CLAWPATCH46_FORT_AVERAGE2COARSE)

/** @copydoc fclaw2d_clawpatch46_fort_average2coarse() */
void FCLAW2D_CLAWPATCH46_FORT_AVERAGE2COARSE(const int* mx, const int* my,
                                             const int* mbc, const int* meqn,
                                             double qcoarse[],double qfine[],
                                             double areacoarse[],double areafine[],
                                             const int* igrid, const int* manifold);



/** @} */

/** @{ @name Ghost Filling *//*-------------------------------------------------------- */

/** Fortran subroutine name */
#define FCLAW2D_CLAWPATCH46_FORT_COPY_FACE \
         FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort_copy_face, \
                        FCLAW2D_CLAWPATCH46_FORT_COPY_FACE)

/** @copydoc fclaw2d_clawpatch46_fort_copy_face() */
void FCLAW2D_CLAWPATCH46_FORT_COPY_FACE(const int* mx, const int* my, const int* mbc, 
                                        const int* meqn,
                                        double qthis[],double qneighbor[], 
                                        const int* a_idir,
                                        struct fclaw2d_patch_transform_data** 
                                        transform_cptr);


/** Fortran subroutine name */
#define FCLAW2D_CLAWPATCH46_FORT_AVERAGE_FACE \
             FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort_average_face, \
                            FCLAW2D_CLAWPATCH46_FORT_AVERAGE_FACE)
/** @copydoc fclaw2d_clawpatch46_fort_average_face() */
void FCLAW2D_CLAWPATCH46_FORT_AVERAGE_FACE(const int* mx, const int* my, const int* mbc,
                                           const int* meqn,
                                           double qcoarse[],double qfine[],
                                           double areacoarse[], double areafine[],
                                           const int* idir, const int* iside,
                                           const int* num_neighbors,
                                           const int* refratio, const int* igrid,
                                           const int* manifold, 
                                           struct fclaw2d_patch_transform_data** 
                                           transform_cptr);
  
/** Fortran subroutine name */
#define FCLAW2D_CLAWPATCH46_FORT_INTERPOLATE_FACE \
              FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort_interpolate_face, \
                             FCLAW2D_CLAWPATCH46_FORT_INTERPOLATE_FACE)
/** @copydoc fclaw2d_clawpatch46_fort_interpolate_face() */
void FCLAW2D_CLAWPATCH46_FORT_INTERPOLATE_FACE(const int* mx, const int* my, 
                                               const int* mbc,const int* meqn,
                                               double qcoarse[],double qfine[],
                                               const int* idir, const int* iside,
                                               const int* num_neighbors,
                                               const int* refratio, const int* igrid,
                                               struct fclaw2d_patch_transform_data** 
                                               transform_cptr);

/** Fortran subroutine name */
#define FCLAW2D_CLAWPATCH46_FORT_COPY_CORNER \
             FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort_copy_corner, \
                            FCLAW2D_CLAWPATCH46_FORT_COPY_CORNER)
/** @copydoc fclaw2d_clawpatch46_fort_copy_corner() */
void FCLAW2D_CLAWPATCH46_FORT_COPY_CORNER(const int* mx, const int* my, const int* mbc,
                                          const int* meqn, double this_q[],
                                          double neighbor_q[],const int* a_corner,
                                          struct fclaw2d_patch_transform_data** 
                                          transform_cptr);

/** Fortran subroutine name */
#define FCLAW2D_CLAWPATCH46_FORT_AVERAGE_CORNER \
      FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort_average_corner, \
                     FCLAW2D_CLAWPATCH46_FORT_AVERAGE_CORNER)
/** @copydoc fclaw2d_clawpatch46_fort_average_corner() */
void FCLAW2D_CLAWPATCH46_FORT_AVERAGE_CORNER(const int* mx, const int* my, const int* mbc,
                                             const int* meqn, const int* a_refratio,
                                             double qcoarse[], double qfine[],
                                             double areacoarse[], double areafine[],
                                             const int* manifold,const int* a_corner, 
                                             struct fclaw2d_patch_transform_data** 
                                             transform_cptr);
  
/** Fortran subroutine name */
#define FCLAW2D_CLAWPATCH46_FORT_INTERPOLATE_CORNER \
      FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort_interpolate_corner, \
                     FCLAW2D_CLAWPATCH46_FORT_INTERPOLATE_CORNER)
/** @copydoc fclaw2d_clawpatch46_fort_interpolate_corner() */
void FCLAW2D_CLAWPATCH46_FORT_INTERPOLATE_CORNER(const int* mx, const int* my, 
                                                 const int* mbc,const int* meqn, 
                                                 const int* a_refratio, double this_q[],
                                                 double neighbor_q[], const int* a_corner,
                                                 struct fclaw2d_patch_transform_data** 
                                                 transform_cptr);

/** @} */


/* ------------------------- Pillow grid block corner routines  ----------------------- */

/** Fortran subroutine name */
#define FCLAW2D_PILLOW46_COPY_BLOCK_CORNER \
               FCLAW_F77_FUNC(fclaw2d_pillow46_copy_block_corner, \
                              FCLAW2D_PILLOW46_COPY_BLOCK_CORNER)

/** @copydoc fclaw2d_pillow46_copy_block_corner() */
void FCLAW2D_PILLOW46_COPY_BLOCK_CORNER(int* mx, int* my,
                                        int* mbc, int* meqn,
                                        double qthis[], 
                                        double qneighbor[], 
                                        int* icorner,
                                        int* iblock);

/** Fortran subroutine name */
#define FCLAW2D_PILLOW46_AVERAGE_BLOCK_CORNER \
          FCLAW_F77_FUNC(fclaw2d_pillow46_average_block_corner,\
                         FCLAW2D_PILLOW46_AVERAGE_BLOCK_CORNER)

/** @copydoc fclaw2d_pillow46_average_block_corner() */
void  FCLAW2D_PILLOW46_AVERAGE_BLOCK_CORNER(int* mx, int* my, int* mbc,
                                            int* meqn, 
                                            int* refratio, 
                                            double qcoarse[],
                                            double qfine[], 
                                            double areacoarse[], 
                                            double areafine[],
                                            int* a_coarse_corner,
                                            int* blockno);

/** Fortran subroutine name */
#define FCLAW2D_PILLOW46_INTERPOLATE_BLOCK_CORNER \
          FCLAW_F77_FUNC(fclaw2d_pillow46_interpolate_block_corner, \
                         FCLAW2D_PILLOW46_INTERPOLATE_BLOCK_CORNER)

/** @copydoc fclaw2d_pillow46_interpolate_block_corner() */
void  FCLAW2D_PILLOW46_INTERPOLATE_BLOCK_CORNER(int* mx, int* my, int* mbc,
                                                int* meqn, int* refratio,
                                                double qcoarse[],
                                                double qfine[], 
                                                int* icoarse_corner,
                                                int* blockno);



/* ------------------------------------ Output functions ------------------------------ */

/** Fortran subroutine name */
#define  FCLAW2D_CLAWPATCH46_FORT_OUTPUT_ASCII \
           FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort_output_ascii, \
                          FCLAW2D_CLAWPATCH46_FORT_OUTPUT_ASCII)
/** @copydoc fclaw2d_clawpatch46_fort_output_ascii() */
void FCLAW2D_CLAWPATCH46_FORT_OUTPUT_ASCII(char* matname1,
                                            int* mx,        int* my,
                                            int* meqn,      int* mbc,
                                            double* xlower, double* ylower,
                                            double* dx,     double* dy,
                                            double q[],
                                            int* patch_num, int* level,
                                            int* blockno,   int* mpirank);

/** Fortran subroutine name */
#define FCLAW2D_CLAWPATCH46_FORT_HEADER_ASCII \
         FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort_header_ascii, \
                        FCLAW2D_CLAWPATCH46_FORT_HEADER_ASCII)
/** @copydoc fclaw2d_clawpatch46_fort_header_ascii() */
void FCLAW2D_CLAWPATCH46_FORT_HEADER_ASCII(const char* matname1, const char* matname2,
                                           const double* time, const int* meqn, 
                                           const int* maux, const int* ngrids);

/* ----------------------------- Diagnostics functions -------------------------------- */

/** Fortran subroutine name */
#define FCLAW2D_CLAWPATCH46_FORT_CONSERVATION_CHECK \
          FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort_conservation_check, \
                         FCLAW2D_CLAWPATCH46_FORT_CONSERVATION_CHECK)
/** @copydoc fclaw2d_clawpatch46_fort_conservation_check() */
void FCLAW2D_CLAWPATCH46_FORT_CONSERVATION_CHECK(int *mx, int *my, int* mbc, int* meqn,
                                                 double *dx, double *dy,
                                                 double* area, double *q, double* sum,
                                                 double* c_kahan);

/** Fortran subroutine name */
#define FCLAW2D_CLAWPATCH46_FORT_COMPUTE_PATCH_AREA \
          FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort_compute_patch_area, \
                         FCLAW2D_CLAWPATCH46_FORT_COMPUTE_PATCH_AREA)

/** @copydoc fclaw2d_clawpatch46_fort_compute_patch_area() */
double FCLAW2D_CLAWPATCH46_FORT_COMPUTE_PATCH_AREA(int *mx, int* my, int*mbc, double* dx,
                                                   double* dy, double area[]);

/** Fortran subroutine name */
#define FCLAW2D_CLAWPATCH46_FORT_COMPUTE_ERROR_NORM \
           FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort_compute_error_norm, \
                          FCLAW2D_CLAWPATCH46_FORT_COMPUTE_ERROR_NORM)

/** @copydoc fclaw2d_clawpatch46_fort_compute_error_norm() */
void FCLAW2D_CLAWPATCH46_FORT_COMPUTE_ERROR_NORM (int* blockno, int* mx,int* my,int* mbc,
                                                  int* meqn,double* dx,double* dy,
                                                  double area[], double error[],
                                                  double* error_norm);

/* ----------------------------- Parallel ghost patches  ------------------------------ */


/** Fortran subroutine name */
#define FCLAW2D_CLAWPATCH46_FORT_LOCAL_GHOST_PACK \
          FCLAW_F77_FUNC(fclaw2d_clawpatch46_fort_local_ghost_pack, \
                         FCLAW2D_CLAWPATCH46_FORT_LOCAL_GHOST_PACK)
/** @copydoc fclaw2d_clawpatch46_fort_local_ghost_pack() */
void FCLAW2D_CLAWPATCH46_FORT_LOCAL_GHOST_PACK(const int *mx, const int *my, 
                                               const int *mbc,
                                               const int *meqn, const int *mint,
                                               double qdata[], double area[],
                                               double qpack[], const int *psize,
                                               const int *packmode, int *ierror);



/* ----------------------------- User convenience headers ----------------------------- */
/** Fortran subroutine name */
#define CLAWPATCH46_TAG4REFINEMENT FCLAW_F77_FUNC(clawpatch46_tag4refinement, \
                                                 CLAWPATCH46_TAG4REFINEMENT)

/** 
 * @brief @copybrief ::clawpatch_fort_tag4refinement_t 
 * 
 * For user defined clawpatch46_tag4refinment subroutine
 * 
 * @details @copydetails ::clawpatch_fort_tag4refinement_t 
 */
void CLAWPATCH46_TAG4REFINEMENT(const int* mx,const int* my,
                               const int* mbc,const int* meqn,
                               const double* xlower, const double* ylower,
                               const double* dx, const double* dy,
                               const int* blockno,
                               double q[],
                               const double* tag_threshold,
                               const int* init_flag,
                               int* tag_patch);



/** Fortran subroutine name */
#define CLAWPATCH46_TAG4COARSENING FCLAW_F77_FUNC(clawpatch46_tag4coarsening, \
                                                CLAWPATCH46_TAG4COARSENING)

/** 
 * @brief @copybrief ::clawpatch_fort_tag4coarsening_t 
 * 
 * For user defined clawpatch46_tag4coarsening subroutine
 * 
 * @details @copydetails ::clawpatch_fort_tag4coarsening_t 
 */
void CLAWPATCH46_TAG4COARSENING(const int* mx, const int* my,
                               const int* mbc, const int* meqn,
                               double xlower[], double ylower[],
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
