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

#ifndef FC2D_CLAWPACK_FORT_H
#define FC2D_CLAWPACK_FORT_H

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

struct fclaw2d_transform_data;  /* Should be replaced by long int?  */

/* -------------------------------- Clawpack functions ------------------------------- */

#define CLAWPACK46_BC2_DEFAULT FCLAW_F77_FUNC(clawpack46_bc2_default,CLAWPACK46_BC2_DEFAULT)
void CLAWPACK46_BC2_DEFAULT(const int* maxmx, const int* maxmy, const int* meqn,
                     const int* mbc, const int* mx, const int* my,
                     const double* xlower, const double* ylower,
                     const double* dx, const double* dy, const double q[],
                     const int* maux, const double aux[], const double* t,
                     const double* dt, const int mthbc[]);


#define CLAWPACK46_FLUX2 FCLAW_F77_FUNC(clawpack46_flux2,CLAWPACK46_FLUX2)
void CLAWPACK46_FLUX2(const int* ixy,const int* maxm, const int* meqn,
                      const int* maux,const int* mbc,const int* mx,
                      double q1d[], double dtdx1d[],
                      double aux1[], double aux2[], double aux3[],
                      double faddm[],double faddp[], double gaddm[],
                      double gaddp[],double cfl1d[], double fwave[],
                      double s[], double amdq[],double apdq[],double cqxx[],
                      double bmasdq[], double bpasdq[],
                      fc2d_clawpack46_rpn2_t rpn2,fc2d_clawpack46_rpt2_t rpt2,
                      const int* mwaves, const int* mcapa,
                      int method[], int mthlim[]);

#define CLAWPACK46_FLUX2FW FCLAW_F77_FUNC(clawpack46_flux2fw,CLAWPACK46_FLUX2FW)
void CLAWPACK46_FLUX2FW(const int* ixy,const int* maxm, const int* meqn, //
                        const int* maux,const int* mbc,const int* mx,
                        double q1d[], double dtdx1d[],
                        double aux1[], double aux2[], double aux3[],
                        double faddm[],double faddp[], double gaddm[],
                        double gaddp[],double cfl1d[], double fwave[],
                        double s[], double amdq[],double apdq[],double cqxx[],
                        double bmasdq[], double bpasdq[],
                        fc2d_clawpack46_rpn2_t rpn2,fc2d_clawpack46_rpt2_t rpt2,
                        const int* mwaves, const int* mcapa,
                        int method[], int mthlim[]);

#define CLAWPACK46_SET_CAPACITY FCLAW_F77_FUNC(clawpack46_set_capacity,CLAWPACK46_SET_CAPACITY)
void CLAWPACK46_SET_CAPACITY(const int* mx, const int *my, const int *mbc,
                             const double *dx, const double* dy, double area[],
                             const int *mcapa, const int* maux, double aux[]);

#define CLAWPACK46_FLUX_ADD FCLAW_F77_FUNC(clawpack46_flux_add, CLAWPACK46_FLUX_ADD)
void CLAWPACK46_FLUX_ADD(const int* mx, const int* my, const int *mbc,
                         const int* meqn, const double* dx, const double *dy,
                         const double *dt, double qnew[],
                         double flux[], const int *iface,
                         double buffer[]);

/* ------------------------------ Time stepping functions ---------------------------- */

#define CLAWPACK46_STEP2_WRAP FCLAW_F77_FUNC(clawpack46_step2_wrap,CLAWPACK46_STEP2_WRAP)
void CLAWPACK46_STEP2_WRAP(const int* maxm, const int* meqn, const int* maux,
                            const int* mbc, const int method[], const int mthlim[],
                            const int* mcapa, const int* mwaves, const int* mx,
                            const int* my, double qold[], double auxold[],
                            const double* dx, const double* dy, const double* dt,
                            const double* cfl, double work[], const int* mwork,
                            const double* xlower, const double* ylower, const int* level,
                            const double* t, double fp[], double fm[], double gp[],
                            double gm[],
                            fc2d_clawpack46_rpn2_t rpn2,
                            fc2d_clawpack46_rpt2_t rpt2,
                            fc2d_clawpack46_flux2_t flux2,
                            int block_corner_count[],int* ierror);


#define FC2D_CLAWPACK46_FORT_TIMEINTERP FCLAW_F77_FUNC (fc2d_clawpack46_fort_timeinterp, \
                                                       FC2D_CLAWPACK46_FORT_TIMEINTERP)
void FC2D_CLAWPACK46_FORT_TIMEINTERP(const int *mx, const int* my, const int* mbc,
                                    const int *meqn, const int* psize,
                                    double qcurr[], double qlast[],
                                    double qinterp[],const double* alpha,
                                    const int* ierror);

/* ------------------------------ Regridding functions --------------------------- */

#define FC2D_CLAWPACK46_FORT_TAG4REFINEMENT FCLAW_F77_FUNC(fc2d_clawpack46_fort_tag4refinement, \
                                                           FC2D_CLAWPACK46_FORT_TAG4REFINEMENT)

void FC2D_CLAWPACK46_FORT_TAG4REFINEMENT(const int* mx,const int* my,
                                         const int* mbc,const int* meqn,
                                         const double* xlower, const double* ylower,
                                         const double* dx, const double* dy,
                                         const int* blockno,
                                         double q[],
                                         const double* tag_threshold,
                                         const int* init_flag,
                                         int* tag_patch);



#define FC2D_CLAWPACK46_FORT_TAG4COARSENING FCLAW_F77_FUNC(fc2d_clawpack46_fort_tag4coarsening, \
                                                          FC2D_CLAWPACK46_FORT_TAG4COARSENING)

void FC2D_CLAWPACK46_FORT_TAG4COARSENING(const int* mx, const int* my,
                                         const int* mbc, const int* meqn,
                                         const double* xlower, const double* ylower,
                                         const double* dx, const double* dy,
                                         const int* blockno,
                                         double q0[],double q1[],
                                         double q2[],double q3[],
                                         const double* tag_threshold,
                                         int* tag_patch);

#define FC2D_CLAWPACK46_FORT_INTERPOLATE2FINE FCLAW_F77_FUNC(fc2d_clawpack46_fort_interpolate2fine, \
                                                FC2D_CLAWPACK46_FORT_INTERPOLATE2FINE)
void FC2D_CLAWPACK46_FORT_INTERPOLATE2FINE(const int* mx,const int* my,
                                           const int* mbc, const int* meqn,
                                           double qcoarse[], double qfine[],
                                           double areacoarse[], double areafine[],
                                           const int* igrid, const int* manifold);

#define FC2D_CLAWPACK46_FORT_AVERAGE2COARSE FCLAW_F77_FUNC(fc2d_clawpack46_fort_average2coarse, \
                                                FC2D_CLAWPACK46_FORT_AVERAGE2COARSE)
void FC2D_CLAWPACK46_FORT_AVERAGE2COARSE(const int* mx, const int* my,
                                         const int* mbc, const int* meqn,
                                         double qcoarse[],double qfine[],
                                         double areacoarse[],double areafine[],
                                         const int* igrid, const int* manifold);



/* ---------------------------------- Ghost filling  -------------------------------------- */

#define FC2D_CLAWPACK46_FORT_COPY_FACE FCLAW_F77_FUNC(fc2d_clawpack46_fort_copy_face, \
                                                     FC2D_CLAWPACK46_FORT_COPY_FACE)

void FC2D_CLAWPACK46_FORT_COPY_FACE(const int* mx, const int* my, const int* mbc, const int* meqn,
                                   double qthis[],double qneighbor[], const int* a_idir,
                                   fclaw2d_transform_data_t** transform_cptr);


#define FC2D_CLAWPACK46_FORT_AVERAGE_FACE FCLAW_F77_FUNC(fc2d_clawpack46_fort_average_face, \
                                                        FC2D_CLAWPACK46_FORT_AVERAGE_FACE)
void FC2D_CLAWPACK46_FORT_AVERAGE_FACE(const int* mx, const int* my, const int* mbc,
                                      const int* meqn,
                                      double qcoarse[],double qfine[],
                                      double areacoarse[], double areafine[],
                                      const int* idir, const int* iside,
                                      const int* num_neighbors,
                                      const int* refratio, const int* igrid,
                                      const int* manifold, fclaw2d_transform_data_t** transform_cptr);

#define FC2D_CLAWPACK46_FORT_INTERPOLATE_FACE FCLAW_F77_FUNC(fc2d_clawpack46_fort_interpolate_face, \
                                                            FC2D_CLAWPACK46_FORT_INTERPOLATE_FACE)
void FC2D_CLAWPACK46_FORT_INTERPOLATE_FACE(const int* mx, const int* my, const int* mbc,
                                          const int* meqn,
                                          double qcoarse[],double qfine[],
                                          const int* idir, const int* iside,
                                          const int* num_neighbors,
                                          const int* refratio, const int* igrid,
                                          fclaw2d_transform_data_t** transform_cptr);

#define FC2D_CLAWPACK46_FORT_COPY_CORNER FCLAW_F77_FUNC(fc2d_clawpack46_fort_copy_corner, \
                                                       FC2D_CLAWPACK46_FORT_COPY_CORNER)
void FC2D_CLAWPACK46_FORT_COPY_CORNER(const int* mx, const int* my, const int* mbc,
                                     const int* meqn, double this_q[],double neighbor_q[],
                                     const int* a_corner,fclaw2d_transform_data_t** transform_cptr);

#define FC2D_CLAWPACK46_FORT_AVERAGE_CORNER FCLAW_F77_FUNC(fc2d_clawpack46_fort_average_corner, \
                                                          FC2D_CLAWPACK46_FORT_AVERAGE_CORNER)
void FC2D_CLAWPACK46_FORT_AVERAGE_CORNER(const int* mx, const int* my, const int* mbc,
                                        const int* meqn, const int* a_refratio,
                                        double qcoarse[], double qfine[],
                                        double areacoarse[], double areafine[],
                                        const int* manifold,
                                        const int* a_corner, fclaw2d_transform_data_t** transform_cptr);

#define FC2D_CLAWPACK46_FORT_INTERPOLATE_CORNER FCLAW_F77_FUNC(fc2d_clawpack46_fort_interpolate_corner, \
                                                             FC2D_CLAWPACK46_FORT_INTERPOLATE_CORNER)
void FC2D_CLAWPACK46_FORT_INTERPOLATE_CORNER(const int* mx, const int* my, const int* mbc,
                                            const int* meqn, const int* a_refratio, double this_q[],
                                            double neighbor_q[], const int* a_corner,
                                            fclaw2d_transform_data_t** transform_cptr);



/* --------------------------------------- Output functions ----------------------------------- */

#define  FC2D_CLAWPACK46_FORT_OUTPUT_ASCII FCLAW_F77_FUNC(fc2d_clawpack46_fort_output_ascii, \
                                                       FC2D_CLAWPACK46_FORT_OUTPUT_ASCII)
void  FC2D_CLAWPACK46_FORT_OUTPUT_ASCII(char* matname1,
                                     int* mx,        int* my,
                                     int* meqn,      int* mbc,
                                     double* xlower, double* ylower,
                                     double* dx,     double* dy,
                                     double q[],
                                     int* patch_num, int* level,
                                     int* blockno,   int* mpirank);

#define FC2D_CLAWPACK46_FORT_HEADER_ASCII FCLAW_F77_FUNC(fc2d_clawpack46_fort_header_ascii, \
                                                        FC2D_CLAWPACK46_FORT_HEADER_ASCII)
void FC2D_CLAWPACK46_FORT_HEADER_ASCII(char* matname1, char* matname2,
                                      double* time, int* meqn, int* maux, int* ngrids);



/* --------------------------------- Diagnostics functions ---------------------------------- */

#define FC2D_CLAWPACK46_FORT_CONSERVATION_CHECK FCLAW_F77_FUNC(fc2d_clawpack46_fort_conservation_check, \
                                                              FC2D_CLAWPACK46_FORT_CONSERVATION_CHECK)
void FC2D_CLAWPACK46_FORT_CONSERVATION_CHECK(int *mx, int *my, int* mbc, int* meqn,
                                            double *dx, double *dy,
                                            double* area, double *q, double* sum);

#define FC2D_CLAWPACK46_FORT_COMPUTE_PATCH_AREA FCLAW_F77_FUNC(fc2d_clawpack46_fort_compute_patch_area, \
                                                              FC2D_CLAWPACK46_FORT_COMPUTE_PATCH_AREA)

double FC2D_CLAWPACK46_FORT_COMPUTE_PATCH_AREA(int *mx, int* my, int*mbc, double* dx,
                                              double* dy, double area[]);



#if 0
#define FC2D_CLAWPACK46_FORT_COMPUTE_ERROR FCLAW_F77_FUNC(fc2d_clawpack46_fort_compute_error, \
                                                     FC2D_CLAWPACK46_FORT_COMPUTE_ERROR)

void FC2D_CLAWPACK46_FORT_COMPUTE_ERROR(int* blockno, int *mx, int *my, int* mbc, int* meqn,
                                        double *dx, double *dy, double *xlower,
                                        double *ylower, double *t, double q[],
                                        double error[]);
#endif

#define FC2D_CLAWPACK46_FORT_COMPUTE_ERROR_NORM FCLAW_F77_FUNC(fc2d_clawpack46_fort_compute_error_norm, \
                                                              FC2D_CLAWPACK46_FORT_COMPUTE_ERROR_NORM)

void FC2D_CLAWPACK46_FORT_COMPUTE_ERROR_NORM(int *mx, int *my, int *mbc, int *meqn,
                                            double *dx, double *dy, double area[],
                                            double error[], double error_norm[]);


/* ---------------------------------- Parallel ghost patches  -------------------------------- */


#define FC2D_CLAWPACK46_FORT_LOCAL_GHOST_PACK FCLAW_F77_FUNC(fc2d_clawpack46_fort_local_ghost_pack, \
                                                             FC2D_CLAWPACK46_FORT_LOCAL_GHOST_PACK)
void  FC2D_CLAWPACK46_FORT_LOCAL_GHOST_PACK(int *mx, int *my, int *mbc,
                                          int *meqn, int *mint,
                                          double qdata[], double area[],
                                          double qpack[], int *psize,
                                          int *packmode, int *ierror);

/* ----------------------------- Misc ClawPack specific functions ------------------------------ */


#define CLAWPACK46_SET_BLOCK FCLAW_F77_FUNC(clawpack46_set_block,CLAWPACK46_SET_BLOCK)
void CLAWPACK46_SET_BLOCK(int* blockno);

#define FC2D_CLAWPACK46_GET_BLOCK FCLAW_F77_FUNC(fc2d_clawpack46_get_block, \
                                                 FC2D_CLAWPACK46_GET_BLOCK)
int FC2D_CLAWPACK46_GET_BLOCK();


#define CLAWPACK46_UNSET_BLOCK FCLAW_F77_FUNC(clawpack46_unset_block, \
                                              CLAWPACK46_UNSET_BLOCK)
void CLAWPACK46_UNSET_BLOCK();




#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif

