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

/* --------------------------------- Clawpack functions ------------------------------- */

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

/* ------------------------------- Time stepping functions ---------------------------- */

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

