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

#ifndef FC3D_CLAWPACK46_H
#define FC3D_CLAWPACK46_H

#include <fclaw_base.h>   /* Needed for FCLAW_F77_FUNC */

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

struct fclaw2d_global;
struct fclaw2d_patch;

typedef  struct fc3d_clawpack46_vtable  fc3d_clawpack46_vtable_t;


/* --------------------------- Clawpack solver functions ------------------------------ */

/* Virtualize clawpack-specific wrapper functions */
typedef void (*clawpack46_src3_t)(struct fclaw2d_global* glob,
								  struct fclaw2d_patch *this_patch,
								  int blockno,
								  int patchno,
								  double t,
								  double dt);
	
typedef void (*clawpack46_b4step3_t)(struct fclaw2d_global* glob,
									 struct fclaw2d_patch *this_patch,
									 int this_block_idx,
									 int this_patch_idx,
									 double t,
									 double dt);

/* ---------------------- Clawpack solver functions (Fortran) ------------------------- */

typedef void (*clawpack46_fort_setprob_t)(void);

typedef void (*clawpack46_fort_bc3_t)(const int* maxmx, const int* maxmy,
                                      const int* maxmz,
									  const int* meqn, const int* mbc,
									  const int* mx, const int* my,
                                      const int* mz,
									  const double* xlower, const double* ylower,
                                      const double* zlower, 
									  const double* dx, const double* dy, 
                                      const double* dz,
									  const double q[], const int* maux,
									  const double aux[], const double* t,
									  const double* dt, const int mthbc[]);

typedef  void (*clawpack46_fort_qinit_t)(const int* maxmx, const int* maxmy,
                                         const int* maxmz, 
										 const int* meqn,const int* mbc,
										 const int* mx, const int* my, const int* mz,
										 const double* xlower, const double* ylower,
                                         const double* zlower,
										 const double* dx, const double* dy,
                                         const double *dz,
										 double q[], const int* maux, double aux[]);

typedef void (*clawpack46_fort_setaux_t)(const int* maxmx, const int* maxmy, 
										 const int* mbc,
										 const int* mx, const int* my,
										 const double* xlower, const double* ylower,
										 const double* dx, const double* dy,
										 const int* maux, double aux[]);

typedef void (*clawpack46_fort_b4step3_t)(const int* maxmx, const int* maxmy,
                                          const int* maxmz,
										  const int* mbc,
										  const int* mx, const int* my, const int* mz,
                                          const int* meqn,
										  double q[], const double* xlower,
										  const double* ylower, const double *zlower,
										  const double* dx, const double* dy,
                                          const double* dz,
										  const double* t, const double* dt,
										  const int* maux, double aux[]);

typedef void (*clawpack46_fort_src3_t)(const int* maxmx, const int* maxmy, 
                                       const int* maxmz, 
									   const int* meqn,
									   const int* mbc, const int* mx,const int* my,
                                       const int* mz, 
									   const double* xlower, const double* ylower,
                                       const double* zlower, 
									   const double* dx, const double* dy, const 
                                       double *dz, double q[],
									   const int* maux, double aux[], const double* t,
									   const double* dt);

typedef void (*clawpack46_fort_rpn2_t)(const int* ixy,const int* maxm, const int* meqn,
									   const int* mwaves, const int* mbc,const int* mx,
									   double ql[], double qr[], double auxl[], 
									   double auxr[],
									   double wave[], double s[],double amdq[], 
									   double apdq[]);


typedef void (*clawpack46_fort_rpt2_t)(const int* ixy, const int* maxm, const int* meqn,
									   const int* mwaves, const int* mbc,const int* mx,
									   double ql[], double qr[], double aux1[], 
									   double aux2[],
									   double aux3[], const int* imp, double dsdq[],
									   double bmasdq[], double bpasdq[]);


typedef void (*clawpack46_fort_flux2_t)(const int* ixy,const int* maxm, const int* meqn,
										const int* maux,const int* mbc,const int* mx,
										double q1d[], double dtdx1d[],
										double aux1[], double aux2[], double aux3[],
										double faddm[],double faddp[], double gaddm[],
										double gaddp[],double cfl1d[], double fwave[],
										double s[], double amdq[],double apdq[], 
										double cqxx[],
										double bmasdq[], double bpasdq[],
										clawpack46_fort_rpn2_t rpn2,
										clawpack46_fort_rpt2_t rpt2,
										const int* mwaves, const int* mcapa,
										int method[], int mthlim[]);

typedef void (*clawpack46_fort_rpn2_cons_t)(const int* meqn, const int* maux, 
											const int *idir, const int* iface, 
                                            double q[], double auxvec_center[],
											double auxvec_edge[], double fq[]);


/* --------------------------------- Virtual table ------------------------------------ */

struct fc3d_clawpack46_vtable
{

	/* C wrappers */
	clawpack46_b4step3_t   b4step3;
	clawpack46_src3_t      src3;

	/* Fortran routines */
	clawpack46_fort_setprob_t   fort_setprob;
	clawpack46_fort_bc3_t       fort_bc3;
	clawpack46_fort_qinit_t     fort_qinit;
	clawpack46_fort_setaux_t    fort_setaux;
	clawpack46_fort_b4step3_t   fort_b4step3;
	clawpack46_fort_src3_t      fort_src3;
	
	clawpack46_fort_rpn2_t      fort_rpn2;
	clawpack46_fort_rpt2_t      fort_rpt2;
	clawpack46_fort_rpn2_cons_t   fort_rpn2_cons;

    clawpack46_fort_flux2_t     flux2;
	
	int is_set;

};

void fc3d_clawpack46_solver_initialize(void);

fc3d_clawpack46_vtable_t* fc3d_clawpack46_vt(void);


/* ----------------------------- User access to solver functions ---------------------- */

void fc3d_clawpack46_setprob(struct fclaw2d_global* glob);


void fc3d_clawpack46_setaux(struct fclaw2d_global* glob,
							struct fclaw2d_patch *this_patch,
							int this_block_idx,
							int this_patch_idx);

void fc3d_clawpack46_set_capacity(struct fclaw2d_global* glob,
								  struct fclaw2d_patch *this_patch,
								  int this_block_idx,
								  int this_patch_idx);

void fc3d_clawpack46_qinit(struct fclaw2d_global* glob,
						   struct fclaw2d_patch *this_patch,
						   int this_block_idx,
						   int this_patch_idx);

void fc3d_clawpack46_b4step3(struct fclaw2d_global* glob,
							 struct fclaw2d_patch *this_patch,
							 int this_block_idx,
							 int this_patch_idx,
							 double t,
							 double dt);

void fc3d_clawpack46_bc3(struct fclaw2d_global *glob,
						 struct fclaw2d_patch *this_patch,
						 int this_block_idx,
						 int this_patch_idx,
						 double t,
						 double dt,
						 int intersects_bc[],
						 int time_interp);

void fc3d_clawpack46_src2(struct fclaw2d_global* glob,
						  struct fclaw2d_patch *this_patch,
						  int this_block_idx,
						  int this_patch_idx,
						  double t,
						  double dt);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif /* !FC3D_CLAWPACH46_H */
