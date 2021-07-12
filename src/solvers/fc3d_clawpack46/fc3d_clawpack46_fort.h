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

#ifndef FC3D_CLAWPACK46_FORT_H
#define FC3D_CLAWPACK46_FORT_H

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/* --------------------------------- Clawpack functions ------------------------------- */

#define CLAWPACK46_BC3_DEFAULT FCLAW_F77_FUNC(clawpack46_bc3_default,CLAWPACK46_BC3_DEFAULT)
void CLAWPACK46_BC3_DEFAULT(const int* maxmx, const int* maxmy, const int* maxmz,
                            const int* meqn,
                            const int* mbc, const int* mx, const int* my,
                            const int* mz,
                            const double* xlower, const double* ylower,
                            const double* zlower, 
                            const double* dx, const double* dy, const double *dz,
                            const double q[],
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
					  clawpack46_fort_rpn2_t rpn2, clawpack46_fort_rpt2_t rpt2,
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
						clawpack46_fort_rpn2_t rpn2,clawpack46_fort_rpt2_t rpt2,
						const int* mwaves, const int* mcapa,
						int method[], int mthlim[]);

#define CLAWPACK46_SET_CAPACITY FCLAW_F77_FUNC(clawpack46_set_capacity,CLAWPACK46_SET_CAPACITY)
void CLAWPACK46_SET_CAPACITY(const int* mx, const int *my, const int *mbc,
							 const double *dx, const double* dy, double area[],
							 const int *mcapa, const int* maux, double aux[]);

/* ------------------------------------- Conservation --------------------------------- */

#define CLAWPACK46_TIME_SYNC_STORE_FLUX FCLAW_F77_FUNC(clawpack46_time_sync_store_flux, \
														 CLAWPACK46_TIME_SYNC_STORE_FLUX)

void CLAWPACK46_TIME_SYNC_STORE_FLUX(int* mx,int* my, int* mbc, int* meqn,
									   int* maux, int* blockno, int* patchno, double* dt,
									   double* el0, 
									   double* el1, 
									   double* el2, 
									   double* el3,
									   double q[], double aux[],
									   double flux0[],double flux1[], 
									   double flux2[], double flux3[],
									   clawpack46_fort_rpn2_cons_t rpn2_cons,
									   double qvec[], 
									   double auxvec_center[], double auxvec_edge[],
									   double flux[]);



#define CLAWPACK46_TIME_SYNC_ACCUMULATE_WAVES \
						FCLAW_F77_FUNC(clawpack46_time_sync_accumulate_waves, \
									   CLAWPACK46_TIME_SYNC_ACCUMULATE_WAVES)

void CLAWPACK46_TIME_SYNC_ACCUMULATE_WAVES(int* mx, int* my, int* mbc, int* meqn,
											  double* dt, double* dx, double* dy, 
											  int* patchno,
											  double el0[], double el1[], 
											  double el2[], double el3[],
											  double fp[], double fm[],
											  double gp[], double gm[],
											  double fp_left[], double fp_right[],
											  double fm_left[], double fm_right[],
											  double gp_bottom[], double gp_top[],
											  double gm_bottom[], double gm_top[]);
	

#define CLAWPACK46_FORT_TIME_SYNC_F2C FCLAW_F77_FUNC(clawpack46_fort_time_sync_f2c, \
													 CLAWPACK46_FORT_TIME_SYNC_F2C)

void  CLAWPACK46_FORT_TIME_SYNC_F2C(const int* mx,const int* my,
									const int *mbc,const int *meqn,
									const int* idir,const int* iface_coarse,
									const int* coarse_blockno, 
									const int* fine_blockno,
									const int* normal_mismatch,
									double areac0[], double areac1[],
									double areac2[], double areac3[],
									double qcoarse[], 
									double fmcoarse0[], double fpcoarse1[],
									double gmcoarse2[], double gpcoarse3[],
									double fmfine0[], double fpfine1[],
									double gmfine2[], double gpfine3[],
									double efc0[], double efc1[],
									double efc2[], double efc3[],
									double eff0[], double eff1[],
									double eff2[], double eff3[],
									double qfine_dummy[],
									struct fclaw2d_patch_transform_data** 
									transform_cptr);


#define CLAWPACK46_FORT_TIME_SYNC_SAMESIZE FCLAW_F77_FUNC(clawpack46_fort_time_sync_samesize, \
													 CLAWPACK46_FORT_TIME_SYNC_SAMESIZE)

void  CLAWPACK46_FORT_TIME_SYNC_SAMESIZE(const int* mx,const int* my,
                                         const int *mbc,const int *meqn,
                                         const int* idir,const int* iface_coarse,
                                         const int* this_blockno, 
                                         const int* neighbor_blockno,
                                         double area0[], double area1[],
                                         double area2[], double area3[],
                                         double qcoarse[], 
                                         double fmcoarse0[], double fpcoarse1[],
                                         double gmcoarse2[], double gpcoarse3[],
                                         double fmfine0[], double fpfine1[],
                                         double gmfine2[], double gpfine3[],
                                         double efc0[], double efc1[],
                                         double efc2[], double efc3[],
                                         double eff0[], double eff1[],
                                         double eff2[], double eff3[],
                                         double qfine_dummy[],
                                         struct fclaw2d_patch_transform_data** 
                                         transform_cptr);

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
							clawpack46_fort_rpn2_t rpn2,
							clawpack46_fort_rpt2_t rpt2,
							clawpack46_fort_flux2_t flux2,
							int block_corner_count[],int* ierror);

/* ----------------------------- Misc ClawPack specific functions ------------------------------ */


#define CLAWPACK46_SET_BLOCK FCLAW_F77_FUNC(clawpack46_set_block,CLAWPACK46_SET_BLOCK)
void CLAWPACK46_SET_BLOCK(int* blockno);

#define FC3D_CLAWPACK46_GET_BLOCK FCLAW_F77_FUNC(fc3d_clawpack46_get_block, \
												 FC3D_CLAWPACK46_GET_BLOCK)
int FC3D_CLAWPACK46_GET_BLOCK();


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

