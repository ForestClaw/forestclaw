/*
  Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
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
#endif

struct fclaw2d_global;
struct fclaw_patch;

typedef  struct fc3d_clawpack46_vtable  fc3d_clawpack46_vtable_t;


/* --------------------------- Clawpack solver functions ------------------------------ */

/* Virtualize clawpack-specific wrapper functions */
typedef void (*clawpack46_src3_t)(struct fclaw2d_global* glob,
								  struct fclaw_patch *this_patch,
								  int blockno,
								  int patchno,
								  double t,
								  double dt);
	
typedef void (*clawpack46_b4step3_t)(struct fclaw2d_global* glob,
									 struct fclaw_patch *this_patch,
									 int this_block_idx,
									 int this_patch_idx,
									 double t,
									 double dt);

/* ---------------------- Clawpack solver functions (Fortran) ------------------------- */

typedef void (*clawpack46_fort_setprob_t)(void);

typedef void (*clawpack46_fort_bc3_t)(const int* meqn, const int* mbc,
									  const int* mx, const int* my, const int* mz,
									  const double* xlower, const double* ylower,
                                      const double* zlower, 
									  const double* dx, const double* dy, 
                                      const double* dz,
									  const double q[], const int* maux,
									  const double aux[], const double* t,
									  const double* dt, const int mthbc[]);

typedef  void (*clawpack46_fort_qinit_t)(const int* meqn,const int* mbc,
										 const int* mx, const int* my, const int* mz,
										 const double* xlower, const double* ylower,
                                         const double* zlower,
										 const double* dx, const double* dy,
                                         const double *dz,
										 double q[], const int* maux, double aux[]);

typedef void (*clawpack46_fort_setaux_t)(const int* mbc,
										 const int* mx, const int* my, const int *mz,
										 const double* xlower, const double* ylower,
                                         const double *zlower,
										 const double* dx, const double* dy,
                                         const double* dz,
										 const int* maux, double aux[]);

typedef void (*clawpack46_fort_b4step3_t)(const int* mbc,
										  const int* mx, const int* my, const int* mz,
                                          const int* meqn,
										  double q[], const double* xlower,
										  const double* ylower, const double *zlower,
										  const double* dx, const double* dy,
                                          const double* dz,
										  const double* t, const double* dt,
										  const int* maux, double aux[]);

typedef void (*clawpack46_fort_src3_t)(const int* meqn,
									   const int* mbc, const int* mx,const int* my,
                                       const int* mz, 
									   const double* xlower, const double* ylower,
                                       const double* zlower, 
									   const double* dx, const double* dy, const 
                                       double *dz, double q[],
									   const int* maux, double aux[], const double* t,
									   const double* dt);

typedef void (*clawpack46_fort_rpn3_t)(const int* ixyz,const int* maxm, 
                                       const int* meqn, const int* mwaves,
                                       const int* maux, const int* mbc,const int* mx, 
                                       double ql[], double qr[],
                                       double auxl[], double auxr[], double wave[],
                                       double s[], double amdq[], double apdq[]);


typedef void (*clawpack46_fort_rpt3_t)(const int* ixyz, const int* icoor, 
                                       const int* imp,
                                       const int *maxm, const int* meqn, 
                                       const int* mwaves, const int *maux,
                                       const int* mbc, const int* mx, 
                                       double ql[], double qr[],
                                       double aux1[], double aux2[], double aux3[], 
                                       double dsdq[], double bmasdq[], double bpasdq[]);

typedef void (*clawpack46_fort_rptt3_t)(const int* ixyz, const int* icoor, 
                                        const int* imp, const int* impt, 
                                        const int* maxm, const int* meqn,
                                        const int* mwaves, const int* maux,
                                        const int* mbc,const int* mx,
                                        double ql[], double qr[],
                                        double aux1[], double aux2[], double aux3[],  
                                        double bsasdq[], 
                                        double cmbsasdq[], double cpbsasdq[]);


typedef void (*clawpack46_fort_flux3_t)(const int* ixyz,const int* maxm, const int* meqn, 
        const int* maux, const int* mbc,const int* mx,
        double q1d[], double dtdx1d[], double dtdy[], double dtdz[],
        double aux1[], double aux2[], double aux3[], 
        double faddm[],    double faddp[], 
        double gadd[],     double hadd[], double *cfl1d,
        double wave[],     double s[], 
        double amdq[],     double apdq[],     double cqxx[],
        double bmamdq[],   double bmapdq[],   double bpamdq[],   double bpapdq[],
        double cmamdq[],   double cmapdq[],   double cpamdq[],   double cpapdq[],
        double cmamdq2[],  double cmapdq2[],  double cpamdq2[],  double cpapdq2[],
        double bmcqxxm[],  double bmcqxxp[],  double bpcqxxm[],  double bpcqxxp[],
        double cmcqxxm[],  double cmcqxxp[], double cpcqxxm[],  double cpcqxxp[],
        double bmcmamdq[], double bmcmapdq[], double bpcmamdq[], double bpcmapdq[],
        double bmcpamdq[], double bmcpapdq[], double bpcpamdq[], double bpcpapdq[],
        clawpack46_fort_rpn3_t  rpn3,
        clawpack46_fort_rpt3_t  rpt3, 
        clawpack46_fort_rptt3_t rptt3,
        const int* mwaves, const int* mcapa,
        int method[], int mthlim[]);

typedef void (*clawpack46_fort_rpn3_cons_t)(const int* meqn, const int* maux, 
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
	
	clawpack46_fort_rpn3_t      fort_rpn3;
	clawpack46_fort_rpt3_t      fort_rpt3;
    clawpack46_fort_rptt3_t     fort_rptt3;
	clawpack46_fort_rpn3_cons_t fort_rpn2_cons;

    clawpack46_fort_flux3_t     flux3;
	
	int is_set;

};

/**
 * @brief Initialize the clawpack46 solver
 * 
 * fclaw2d_vtables_intialize should be called before this function.
 * 
 * fclaw3dx_clawpatch_options, and fc3d_clawpack46_options should be stored in glob.
 * fc3d_clawpack46_options will be changed in this call.
 * 
 * @param glob the global context
 */
void fc3d_clawpack46_solver_initialize(struct fclaw2d_global* glob);

/**
 * @brief Get the clawpack46 vtable
 * 
 * @param glob the global context
 * @return fc3d_clawpack46_vtable_t* the vtable
 */
fc3d_clawpack46_vtable_t* fc3d_clawpack46_vt(struct fclaw2d_global* glob);


/* ----------------------------- User access to solver functions ---------------------- */

void fc3d_clawpack46_setprob(struct fclaw2d_global* glob);


void fc3d_clawpack46_setaux(struct fclaw2d_global* glob,
							struct fclaw_patch *this_patch,
							int this_block_idx,
							int this_patch_idx);

void fc3d_clawpack46_set_capacity(struct fclaw2d_global* glob,
								  struct fclaw_patch *this_patch,
								  int this_block_idx,
								  int this_patch_idx);

void fc3d_clawpack46_qinit(struct fclaw2d_global* glob,
						   struct fclaw_patch *this_patch,
						   int this_block_idx,
						   int this_patch_idx);

void fc3d_clawpack46_b4step3(struct fclaw2d_global* glob,
							 struct fclaw_patch *this_patch,
							 int this_block_idx,
							 int this_patch_idx,
							 double t,
							 double dt);

void fc3d_clawpack46_bc3(struct fclaw2d_global *glob,
						 struct fclaw_patch *this_patch,
						 int this_block_idx,
						 int this_patch_idx,
						 double t,
						 double dt,
						 int intersects_bc[],
						 int time_interp);

void fc3d_clawpack46_src2(struct fclaw2d_global* glob,
						  struct fclaw_patch *this_patch,
						  int this_block_idx,
						  int this_patch_idx,
						  double t,
						  double dt);


#ifdef __cplusplus
}
#endif


#endif /* !FC3D_CLAWPACH46_H */
