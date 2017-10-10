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

#ifndef FC2D_CLAWPACK5_FORT_H
#define FC2D_CLAWPACK5_FORT_H

#include <fclaw_base.h>   /* Needed for FCLAW_F77_FUNC */

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/* -------------------------------- Clawpack functions ------------------------------- */


#define CLAWPACK5_BC2_DEFAULT FCLAW_F77_FUNC(clawpack5_bc2_default, CLAWPACK5_BC2_DEFAULT)
void CLAWPACK5_BC2_DEFAULT(const int* meqn, const int* mbc,
                           const int* mx, const int* my,
                           const double* xlower, const double* ylower,
                           const double* dx, const double* dy,
                           const double q[], const int* maux,
                           const double aux[], const double* t,
                           const double* dt, const int mthbc[]);


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
#if 0
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
                        int method[], int mthlim[]);
#endif                        

#define CLAWPACK5_SET_CAPACITY FCLAW_F77_FUNC(clawpack5_set_capacity,CLAWPACK5_SET_CAPACITY)
void CLAWPACK5_SET_CAPACITY(const int* mx, const int *my, const int *mbc,
                             const double *dx, const double* dy, double area[],
                             const int *mcapa, const int* maux, double aux[]);


/* ------------------------------ Time stepping functions ---------------------------- */


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

/* ----------------------------- Misc ClawPack specific functions ------------------------------ */
    
#define CLAWPACK5_SET_BLOCK FCLAW_F77_FUNC(clawpack5_set_block,CLAWPACK5_SET_BLOCK)
void CLAWPACK5_SET_BLOCK(int* blockno);

#define FC2D_CLAWPACK5_GET_BLOCK FCLAW_F77_FUNC(fc2d_clawpack5_get_block, \
                                                 FC2D_CLAWPACK5_GET_BLOCK)
int FC2D_CLAWPACK5_GET_BLOCK();


#define CLAWPACK5_UNSET_BLOCK FCLAW_F77_FUNC(clawpack5_unset_block, \
                                              CLAWPACK5_UNSET_BLOCK)
void CLAWPACK5_UNSET_BLOCK();

#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif /* !FC2D_CLAWPACH5_H */
