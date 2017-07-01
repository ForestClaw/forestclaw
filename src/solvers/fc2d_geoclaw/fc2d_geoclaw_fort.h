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

#ifndef FCLAW2D_GEOCLAW_FORT_H
#define FCLAW2D_GEOCLAW_FORT_H

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

struct fclaw2d_transform_data;  /* Should be replaced by long int?  */

#define GEOCLAW_BC2 FCLAW_F77_FUNC(geoclaw_bc2,GEOCLAW_BC2)
void GEOCLAW_BC2(const int* meqn, const int* mbc,
                   const int* mx, const int* my,
                   const double* xlower, const double* ylower,
                   const double* dx, const double* dy,
                   const double q[], const int* maux,
                   const double aux[], const double* t,
                   const double* dt, const int mthbc[]);

/* --------------------------------------------------------------------
   Classic routines
   - These are provided only for convenience;  these files are not
   compiled into the library, but will be provided by the user.
   -------------------------------------------------------------------- */

/* Macros for C/Fortran portability */
#define SETPROB FCLAW_F77_FUNC(setprob,SETPROB)
#define GEOCLAW_QINIT   FCLAW_F77_FUNC(geoclaw_qinit,GEOCLAW_QINIT)
#define GEOCLAW_SETAUX  FCLAW_F77_FUNC(geoclaw_setaux,GEOCLAW_SETAUX)
#define GEOCLAW_B4STEP2 FCLAW_F77_FUNC(geoclaw_b4step2,GEOCLAW_B4STEP2)
#define GEOCLAW_SRC2    FCLAW_F77_FUNC(geoclaw_src2,GEOCLAW_SRC2)
#define BC2             FCLAW_F77_FUNC(bc2,BC2)
#define GEOCLAW_RPN2    FCLAW_F77_FUNC(geoclaw_rpn2,GEOCLAW_RPN2)
#define GEOCLAW_RPT2    FCLAW_F77_FUNC(geoclaw_rpt2,GEOCLAW_RPT2)

/* Specific to geoclaw */
#define GEOCLAW_SET_MODULES   FCLAW_F77_FUNC(geoclaw_set_modules, \
                                             GEOCLAW_SET_MODULES)
void GEOCLAW_SET_MODULES(const int* mwaves_in, const int* mcapa_in,
                         const int* meqn_in, const int* maux_in,
                         const int mthlim_in[], const int method_in[],
                         const double *ax, const double *bx, 
                         const double *ay,
            const double *by);

void SETPROB();

void GEOCLAW_QINIT(const int* meqn,const int* mbc,
                   const int* mx, const int* my,
                   const double* xlower, const double* ylower,
                   const double* dx, const double* dy,
                   double q[], const int* maux, double aux[]);

void GEOCLAW_SETAUX(const int* mbc,
                    const int* mx, const int* my,
                    const double* xlower, const double* ylower,
                    const double* dx, const double* dy,
                    const int* maux, double aux[],
                    const int* is_ghost, const int* nghost,
                    const int* mint);

void BC2(const int* meqn, const int* mbc,
         const int* mx, const int* my,
         const double* xlower, const double* ylower,
         const double* dx, const double* dy,
         const double q[], const int* maux,
         const double aux[], const double* t,
         const double* dt, const int mthbc[]);

void GEOCLAW_B4STEP2(const int* mbc,
                     const int* mx, const int* my, const int* meqn,
                     double q[], const double* xlower,
                     const double* ylower,
                     const double* dx, const double* dy,
                     const double* t, const double* dt,
                     const int* maux, double aux[]);

void GEOCLAW_SRC2(const int* meqn,
                  const int* mbc, const int* mx,const int* my,
                  const double* xlower, const double* ylower,
                  const double* dx, const double* dy, double q[],
                  const int* maux, double aux[], const double* t,
                  const double* dt);

/* Riemann solvers */
void GEOCLAW_RPN2(const int* ixy,const int* maxm, const int* meqn,
                  const int* mwaves, const int* maux,
                  const int* mbc,const int* mx,
                  double ql[], double qr[], double auxl[], double auxr[],
                  double wave[], double s[],double amdq[], double apdq[]);

void GEOCLAW_RPT2(const int* ixy, const int* imp, const int* maxm, const int* meqn,
                  const int* mwaves, const int* maux, const int* mbc,const int* mx,
                  double ql[], double qr[], double aux1[], double aux2[],
                  double aux3[],  double asdq[],
                  double bmasdq[], double bpasdq[]);

/* --------------------------------------------------------------------
   Time stepping
   -------------------------------------------------------------------- */

#define GEOCLAW_STEP2_WRAP FCLAW_F77_FUNC(geoclaw_step2_wrap,GEOCLAW_STEP2_WRAP)
void GEOCLAW_STEP2_WRAP(const int* maxm, const int* meqn, const int* maux,
                        const int* mbc, const int method[], const int mthlim[],
                        const int* mcapa, const int* mwaves, const int* mx,
                        const int* my, double qold[], double auxold[],
                        const double* dx, const double* dy, const double* dt,
                        const double* cfl, double work[], const int* mwork,
                        const double* xlower, const double* ylower, const int* level,
                        const double* t, double fp[], double fm[], double gp[],
                        double gm[],
                        fc2d_geoclaw_rpn2_t rpn2,
                        fc2d_geoclaw_rpt2_t rpt2,
                        int block_corner_count[]);

#define GEOCLAW_STEP2 FCLAW_F77_FUNC(geoclaw_step2,GEOCLAW_STEP2)
void GEOCLAW_STEP2(const int* maxm, const int* meqn, const int* maux,
                   const int* mbc, const int* mx,
                   const int* my, double qold[], double aux[],
                   const double* dx, const double* dy, const double* dt,
                   const double* cflgrid, double fm[], double fp[], double gm[],
                   double gp[],
                   fc2d_geoclaw_rpn2_t rpn2,
                   fc2d_geoclaw_rpt2_t rpt2);

#define GEOCLAW_FLUX2 FCLAW_F77_FUNC(geoclaw_flux2,GEOCLAW_FLUX2)
void GEOCLAW_FLUX2(const int* ixy,const int* maxm, const int* meqn,
                   const int* maux,const int* mbc,const int* mx,
                   double q1d[], double dtdx1d[],
                   double aux1[], double aux2[], double aux3[],
                   double faddm[],double faddp[], double gaddm[],
                   double gaddp[],double cfl1d[], double wave[],
                   double s[], double amdq[],double apdq[],double cqxx[],
                   double bmasdq[], double bpasdq[],
                   fc2d_geoclaw_rpn2_t rpn2,fc2d_geoclaw_rpt2_t rpt2);
/*
#define geoclaw_FLUX2FW FCLAW_F77_FUNC(geoclaw_flux2fw,geoclaw_FLUX2FW)
void GEOCLAW_FLUX2FW(const int* ixy,const int* maxm, const int* meqn, //
                        const int* maux,const int* mbc,const int* mx,
                        double q1d[], double dtdx1d[],
                        double aux1[], double aux2[], double aux3[],
                        double faddm[],double faddp[], double gaddm[],
                        double gaddp[],double cfl1d[], double fwave[],
                        double s[], double amdq[],double apdq[],double cqxx[],
                        double bmasdq[], double bpasdq[],
                        fc2d_geoclaw_rpn2_t rpn2,fc2d_geoclaw_rpt2_t rpt2,
                        const int* mwaves, const int* mcapa,
                        int method[], int mthlim[]);*/

#define GEOCLAW_SET_CAPACITY FCLAW_F77_FUNC(geoclaw_set_capacity,GEOCLAW_SET_CAPACITY)
void GEOCLAW_SET_CAPACITY(const int* mx, const int *my, const int *mbc,
                          const double *dx, const double* dy, double area[],
                          const int *mcapa, const int* maux, double aux[]);


#define GEOCLAW_SET_BLOCK FCLAW_F77_FUNC(geoclaw_set_block,GEOCLAW_SET_BLOCK)
void GEOCLAW_SET_BLOCK(int* blockno);

#define FC2D_GEOCLAW_GET_BLOCK FCLAW_F77_FUNC(fc2d_geoclaw_get_block, \
                                                 FC2D_GEOCLAW_GET_BLOCK)
int FC2D_GEOCLAW_GET_BLOCK();


#define GEOCLAW_UNSET_BLOCK FCLAW_F77_FUNC(geoclaw_unset_block, \
                                              GEOCLAW_UNSET_BLOCK)
void GEOCLAW_UNSET_BLOCK();

#define GEOCLAW_GAUGES_GETNUM FCLAW_F77_FUNC(geoclaw_gauges_getnum, \
                                             GEOCLAW_GAUGES_GETNUM)
int GEOCLAW_GAUGES_GETNUM(char fname[]);

#define GEOCLAW_GAUGES_INIT FCLAW_F77_FUNC(geoclaw_gauges_init,         \
                                           GEOCLAW_GAUGES_INIT)
void GEOCLAW_GAUGES_INIT(const int* restart, const int* meqn, const int* num_gauges,
                         struct geoclaw_gauge *gauges, char fname[]);

#define GEOCLAW_UPDATE_GAUGE FCLAW_F77_FUNC(geoclaw_update_gauge, \
                                            GEOCLAW_UPDATE_GAUGE)
void GEOCLAW_UPDATE_GAUGE (int* mx,int* my,int* mbc,int* meqn,double* xlower,
                           double* ylower,double* dx,double* dy,double q[],
                           int* maux,double aux[],double* xc,double* yc,double var[],
                           double* eta);

#define GEOCLAW_TOPO_UPDATE FCLAW_F77_FUNC(geoclaw_topo_update, \
                                           GEOCLAW_TOPO_UPDATE)
void GEOCLAW_TOPO_UPDATE (double* t);
/************************ Regridding ******************************/
#define FC2D_GEOCLAW_FORT_TAG4REFINEMENT FCLAW_F77_FUNC(fc2d_geoclaw_fort_tag4refinement, \
                                                        FC2D_GEOCLAW_FORT_TAG4REFINEMENT)
void FC2D_GEOCLAW_FORT_TAG4REFINEMENT(int* mx,int* my, int* mbc, int *meqn,int*maux,
                                      double * xlower,double * ylower,
                                      double* dx,double* dy, double*t, int* blockno,
                                      double q[], double aux[], int* mbathy, int* level,int* maxlevel,
                                      int* init_flag,int* tag_patch);

#define FC2D_GEOCLAW_FORT_TAG4COARSENING FCLAW_F77_FUNC(fc2d_geoclaw_fort_tag4coarsening, \
                                                        FC2D_GEOCLAW_FORT_TAG4COARSENING)
void FC2D_GEOCLAW_FORT_TAG4COARSENING(int* blockno,int* mx,int* my,int* mbc,int* meqn, int* maux,
                                      double xlower[],double ylower[],double* dx,double* dy,
                                      double* t, double q0[], double q1[], double q2[], double q3[],
                                      double aux0[], double aux1[], double aux2[], double aux3[],
                                      int* mbathy, int* level, int* maxlevel, double* dry_tolerance_c,
                                      double* wave_tolerance_c, int* speed_tolerance_entries_c,
                                      double speed_tolerance_c[], int* tag_patch);

#define FC2D_GEOCLAW_FORT_INTERPOLATE2FINE FCLAW_F77_FUNC(fc2d_geoclaw_fort_interpolate2fine, \
                                                          FC2D_GEOCLAW_FORT_INTERPOLATE2FINE)
void FC2D_GEOCLAW_FORT_INTERPOLATE2FINE(int* mx,int* my,int* mbc,int* meqn, double qcoarse[],
                                        double qfine[], int* maux, double aux_coarse[],
                                        double aux_fine[], int* mbathy, int* p4est_refineFactor,int* refratio,
                                        int* igrid);

#define FC2D_GEOCLAW_FORT_AVERAGE2COARSE FCLAW_F77_FUNC(fc2d_geoclaw_fort_average2coarse, \
                                                        FC2D_GEOCLAW_FORT_AVERAGE2COARSE)
void FC2D_GEOCLAW_FORT_AVERAGE2COARSE(int* mx,int* my,int* mbc,int* meqn, double qcoarse[],
                                      double qfine[], int* maux, double aux_coarse[],
                                      double aux_fine[], int* mcapa, int* mbathy, int* p4est_refineFactor,
                                      int* refratio, int* igrid);

#define FC2D_GEOCLAW_FORT_WRITE_HEADER FCLAW_F77_FUNC(fc2d_geoclaw_fort_write_header,\
                                                      FC2D_GEOCLAW_FORT_WRITE_HEADER)
void FC2D_GEOCLAW_FORT_WRITE_HEADER(int* iframe,double* time,int* meqn,int* maux,int* ngrids);

#define FC2D_GEOCLAW_FORT_WRITE_FILE FCLAW_F77_FUNC(fc2d_geoclaw_fort_write_file, \
                                                    FC2D_GEOCLAW_FORT_WRITE_FILE)
void FC2D_GEOCLAW_FORT_WRITE_FILE(int* mx,int* my,int* meqn,int* maux,int* mbathy,int* mbc,
                             double* xlower,double* ylower,double* dx,double* dy,
                             double q[],double aux[],int* iframe,int* patch_num,int* level,
                             int* blockno,int* mpirank);


#define FC2D_GEOCLAW_FORT_COPY_FACE FCLAW_F77_FUNC(fc2d_geoclaw_fort_copy_face, \
                                                     FC2D_GEOCLAW_FORT_COPY_FACE)

void FC2D_GEOCLAW_FORT_COPY_FACE(const int* mx, const int* my, const int* mbc, const int* meqn,
                                   double qthis[],double qneighbor[], const int* a_idir,
                                   struct fclaw2d_transform_data** transform_cptr);


#define FC2D_GEOCLAW_FORT_AVERAGE_FACE FCLAW_F77_FUNC(fc2d_geoclaw_fort_average_face, \
                                                        FC2D_GEOCLAW_FORT_AVERAGE_FACE)
void FC2D_GEOCLAW_FORT_AVERAGE_FACE(const int* mx,const int* my,const int* mbc,const int* meqn,
                                    double qcoarse[],double qfine[],const int* maux,
                                    double auxcoarse[],double auxfine[],const int* mcapa,
                                    const int* mbathy,const int* idir,const int* iface_coarse,
                                    const int* p4est_refineFactor,const int* refratio,
                                    const int* igrid,const int* manifold,
                                    struct fclaw2d_transform_data** transform_data);

#define FC2D_GEOCLAW_FORT_INTERPOLATE_FACE FCLAW_F77_FUNC(fc2d_geoclaw_fort_interpolate_face, \
                                                            FC2D_GEOCLAW_FORT_INTERPOLATE_FACE)
void FC2D_GEOCLAW_FORT_INTERPOLATE_FACE(const int* mx, const int* my, const int* mbc,
                                          const int* meqn,
                                          double qcoarse[],double qfine[],
                                          const int* maux, double aux_coarse[],
                                          double aux_fine[], const int* mbathy,
                                          const int* idir, const int* iside,
                                          const int* num_neighbors,
                                          const int* refratio, const int* igrid,
                                          struct fclaw2d_transform_data** transform_cptr);


#define FC2D_GEOCLAW_FORT_COPY_CORNER FCLAW_F77_FUNC(fc2d_geoclaw_fort_copy_corner, \
                                                       FC2D_GEOCLAW_FORT_COPY_CORNER)
void FC2D_GEOCLAW_FORT_COPY_CORNER(const int* mx, const int* my, const int* mbc,
                                   const int* meqn, double this_q[],double neighbor_q[],
                                   const int* a_corner,
                                   struct fclaw2d_transform_data** transform_cptr);

#define FC2D_GEOCLAW_FORT_AVERAGE_CORNER FCLAW_F77_FUNC(fc2d_geoclaw_fort_average_corner, \
                                                          FC2D_GEOCLAW_FORT_AVERAGE_CORNER)
void FC2D_GEOCLAW_FORT_AVERAGE_CORNER(const int* mx, const int* my, const int* mbc,
                                      const int* meqn, const int* a_refratio,
                                      double qcoarse[], double qfine[], const int* maux,
                                      double auxcoarse[], double auxfine[], const int* mcapa,
                                      const int* mbathy, const int* manifold,
                                      const int* a_corner, 
                                      struct fclaw2d_transform_data** transform_cptr);
    
#define FC2D_GEOCLAW_FORT_INTERPOLATE_CORNER FCLAW_F77_FUNC(fc2d_geoclaw_fort_interpolate_corner, \
                                                             FC2D_GEOCLAW_FORT_INTERPOLATE_CORNER)
void FC2D_GEOCLAW_FORT_INTERPOLATE_CORNER(const int* mx, const int* my, const int* mbc,
                                          const int* meqn, const int* a_refratio, double qcoarse[],
                                          double qfine[], const int* maux, double aux_coarse[],
                                          double aux_fine[], const int* mbathy, const int* a_corner,
                                          struct fclaw2d_transform_data** transform_cptr);



#define FC2D_GEOCLAW_FORT_GHOSTPACKAUX FCLAW_F77_FUNC(fc2d_geoclaw_fort_ghostpackaux, \
                                                     FC2D_GEOCLAW_FORT_GHOSTPACKAUX)
void  FC2D_GEOCLAW_FORT_GHOSTPACKAUX(int *mx, int *my, int *mbc,
                                     int *maux, int *mint,
                                     double auxdata[], double auxpack[],
                                     int *auxsize, int *packmode, int *ierror);



#define FC2D_GEOCLAW_FORT_CONSERVATION_CHECK FCLAW_F77_FUNC(fc2d_geoclaw_fort_conservation_check, \
                                                              FC2D_GEOCLAW_FORT_CONSERVATION_CHECK)

void FC2D_GEOCLAW_FORT_CONSERVATION_CHECK(int *mx, int *my, int* mbc, int* meqn,
                                            double *dx, double *dy,
                                            double* area, double *q, double* sum);

#define FC2D_GEOCLAW_FORT_COMPUTE_PATCH_AREA FCLAW_F77_FUNC(fc2d_geoclaw_fort_compute_patch_area, \
                                                              FC2D_GEOCLAW_FORT_COMPUTE_PATCH_AREA)

double FC2D_GEOCLAW_FORT_COMPUTE_PATCH_AREA(int *mx, int* my, int*mbc, double* dx,
                                              double* dy, double area[]);


#define FC2D_GEOCLAW_FORT_COMPUTE_ERROR_NORM FCLAW_F77_FUNC(fc2d_geoclaw_fort_compute_error_norm, \
                                                              FC2D_GEOCLAW_FORT_COMPUTE_ERROR_NORM)

void FC2D_GEOCLAW_FORT_COMPUTE_ERROR_NORM(int *mx, int *my, int *mbc, int *meqn,
                                            double *dx, double *dy, double area[],
                                            double error[], double error_norm[]);


#define FC2D_GEOCLAW_FORT_GHOSTPACK_QAREA FCLAW_F77_FUNC(fc2d_geoclaw_fort_ghostpack_qarea, \
                                                           FC2D_GEOCLAW_FORT_GHOSTPACK_QAREA)
void  FC2D_GEOCLAW_FORT_GHOSTPACK_QAREA(int *mx, int *my, int *mbc,
                                          int *meqn, int *mint,
                                          double qdata[], double area[],
                                          double qpack[], int *psize,
                                          int *packmode, int *ierror);

#define FC2D_GEOCLAW_FORT_TIMEINTERP FCLAW_F77_FUNC (fc2d_geoclaw_fort_timeinterp, \
                                                       FC2D_GEOCLAW_FORT_TIMEINTERP)
void FC2D_GEOCLAW_FORT_TIMEINTERP(const int *mx, const int* my, const int* mbc,
                                    const int *meqn, const int* psize,
                                    double qcurr[], double qlast[],
                                    double qinterp[],const double* alpha,
                                    const int* ierror);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif

