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

#ifndef FCLAW2D_CLAWPATCH_FORT3_H
#define FCLAW2D_CLAWPATCH_FORT3_H

#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix syntax highlighting */
#endif

struct fclaw2d_global;
struct fclaw2d_patch;

struct fclaw2d_patch_transform_data;  /* Should be replaced by long int?  */

/* Functions defined here are implemented in individual solvers (clawpack 4.6 and 
   clawpack 5.0) */


/* --------------------------- Ghost filling - patch specific ------------------------- */

typedef void (*clawpatch_fort_copy_face_t)(const int* mx, const int* my, const int* mz,
                                           const int* mbc, const int* meqn,
										   double qthis[],double qneighbor[], 
										   const int* a_idir,
										   struct fclaw2d_patch_transform_data** transform_cptr);

typedef void (*clawpatch_fort_average_face_t)(const int* mx, const int* my, const int* mz,
                                              const int* mbc, const int* meqn,
											  double qcoarse[],double qfine[],
											  double areacoarse[], double areafine[],
											  const int* idir, const int* iside,
											  const int* num_neighbors,
											  const int* refratio, const int* igrid,
											  const int* manifold, 
											  struct fclaw2d_patch_transform_data** transform_cptr);
	
typedef void (*clawpatch_fort_interpolate_face_t)(const int* mx, const int* my, 
                                                  const int* mz,
                                                  const int* mbc, const int* meqn,
												  double qcoarse[],double qfine[],
												  const int* idir, const int* iside,
												  const int* num_neighbors,
												  const int* refratio, const int* igrid,
												  struct fclaw2d_patch_transform_data** transform_cptr);
	
	

typedef void (*clawpatch_fort_copy_corner_t)(const int* mx, const int* my, const int* mz,
                                             const int* mbc,
											 const int* meqn, double this_q[],
                                             double neighbor_q[],
											 const int* a_corner,
											 struct fclaw2d_patch_transform_data** 
                                             transform_cptr);

typedef void (*clawpatch_fort_average_corner_t)(const int* mx, const int* my, 
                                                const int* mz,
                                                const int* mbc,
												const int* meqn, const int* a_refratio,
												double qcoarse[], double qfine[],
												double areacoarse[], double areafine[],
												const int* manifold,
												const int* a_corner, 
												struct fclaw2d_patch_transform_data** transform_cptr);

typedef void (*clawpatch_fort_interpolate_corner_t)(const int* mx, const int* my, 
                                                    const int* mz,
                                                    const int* mbc,
													const int* meqn, const int* a_refratio, 
													double this_q[],
													double neighbor_q[], const int* a_corner,
													struct fclaw2d_patch_transform_data** transform_cptr);
	

/* --------------------------------- Regridding functions ----------------------------- */

typedef void (*clawpatch_fort_tag4refinement_t)(const int* mx,const int* my,const int* mz,
												const int* mbc,const int* meqn,
												const double* xlower, const double* ylower,
                                                const double* zlower, 
												const double* dx, const double* dy, 
                                                const double* dz,
												const int* blockno,
												double q[],
												const double* tag_threshold,
												const int* init_flag,
												int* tag_patch);

typedef void (*clawpatch_fort_tag4coarsening_t)(const int* mx, 
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
												double q0[], double q1[], 
                                                double q2[], double q3[],
												const double* tag_threshold,
                                                const int* init_flag,
												int* tag_patch);


#if 0
/* Even though this is for 3d patches, we still assume that tagging can only be
   dependent on two dimensions */
typedef int (*clawpatch_fort_exceeds_threshold_t)(int *blockno,
                                                  double *qval, 
                                                  double *qmin, 
                                                  double *qmax,
                                                  double quad[], 
                                                  double *dx, 
                                                  double *dy, 
                                                  double *xc, 
                                                  double *yc, 
                                                  double *tag_threshold,
                                                  int    *init_flag,
                                                  int    *is_ghost);
#endif
    

typedef void (*clawpatch_fort_interpolate2fine_t)(const int* mx, const int* my, 
                                                  const int* mz,
												  const int* mbc, const int* meqn,
												  double qcoarse[], double qfine[],
												  double areacoarse[], double areafine[],
												  const int* igrid, const int* manifold);
	
typedef void (*clawpatch_fort_average2coarse_t)(const int* mx, const int* my, 
                                                const int* mz,
												const int* mbc, const int* meqn,
												double qcoarse[],double qfine[],
												double areacoarse[],double areafine[],
												const int* igrid, const int* manifold);
	

/* ----------------------------------- time stepping ---------------------------------- */

typedef void (*clawpatch_fort_timeinterp_t)(const int *mx, const int* my, 
                                            const int* mz, const int* mbc,
											const int *meqn, const int* psize,
											double qcurr[], double qlast[],
											double qinterp[],const double* alpha,
											const int* ierror);
	
/* ------------------------------- Parallel ghost patches ----------------------------- */

typedef void (*clawpatch_fort_local_ghost_pack_t)(const int *mx, const int *my, 
                                                  const int* mz, const int *mbc,
												  const int *meqn, const int *mint,
												  double qdata[], double area[],
												  double qpack[], int *psize,
												  const int *packmode, int *ierror);
	
/* ---------------------------------- Output functions -------------------------------- */

#if 0
typedef void  (*clawpatch_fort_header_ascii_t)(const char* matname1,const char* matname2,
											   const double* time, const int* meqn, 
                                               const int* maux, const int* ngrids);

#endif

/* Write out data */
typedef void (*clawpatch_fort_output_ascii_t)(const char* matname1,
											  const int* mx, const int* my, const int* mz,
											  const int* meqn, const int* mbc,
											  const double* xlower, const double* ylower, 
                                              const double* zlower,
											  const double* dx, const double* dy,
                                              const double* dz,
											  double q[],
											  const int* patch_num, const int* level,
											  const int* blockno,   const int* mpirank);


/* ----------------------------- Diagnostic functions --------------------------------- */

typedef void (*clawpatch_fort_error_t)(const int* blockno, const int *mx, const int *my, 
                                       const int* mz, 
                                       const int *mbc, const int *meqn,
									   const double *dx, const double *dy, const double *dz,
                                       const double *xlower, const double *ylower, 
                                       const double *zlower,
                                       double *t, double q[],
									   double error[], double soln[]);

typedef void (*clawpatch_fort_conscheck_t)(const int *mx, const int *my, const int *mz,
                                           const int* mbc, const int* meqn,
										   const double *dx, const double *dy,
                                           const double *dz, 
										   double area[], double q[], double sum[],
                                           double *c_kahan);

typedef double (*clawpatch_fort_area_t)(const int *mx, const int* my, const int *mz,
                                        const int *mbc, const double* dx,
										const double* dy, const double *dz, double area[]);

typedef void (*clawpatch_fort_norm_t)(const int* blockno, const int *mx, const int *my,
                                      const int *mz, const  int *mbc, const int *meqn,
									  const double *dx, const double *dy, const double *dz, 
                                      double area[],
									  double error[], double error_norm[]);



/* -------------------------- User convenience headers -------------------------------- */

#define TAG4REFINEMENT FCLAW_F77_FUNC(tag4refinement,TAG4REFINEMENT)
void TAG4REFINEMENT(const int* mx,const int* my, const int* mz,
					const int* mbc,const int* meqn,
					const double* xlower, const double* ylower,
                    const double* zlower, 
					const double* dx, const double* dy, const double *dz,
					const int* blockno,
					double q[],
					const double* tag_threshold,
					const int* init_flag,
					int* tag_patch);

#define TAG4COARSENING FCLAW_F77_FUNC(tag4coarsening,TAG4COARSENING)
void TAG4COARSENING(const int* mx, const int* my, const int* mz,
					const int* mbc, const int* meqn,
					const double* xlower, const double* ylower,
                    const double *zlower,
					const double* dx, const double* dy, const double *dz,
					const int* blockno,
					double q0[], double q1[],
                    double q2[], double q3[],
					const double* tag_threshold,
                    const int* initflag,
					int* tag_patch);

/* ----------------------------- interpolation/coarsening ----------------------------- */

#define INTERPOLATE2FINE FCLAW_F77_FUNC(interpolate2fine, INTERPOLATE2FINE)
void INTERPOLATE2FINE(const int* mx,const int* my, const int* mz, const int* mbc,
					  const int* meqn, double qcoarse[], double qfine[],
					  double areacoarse[], double areafine[],
					  const int* igrid,const int* manifold);

#define AVERAGE2COARSE FCLAW_F77_FUNC(average2coarse, AVERAGE2COARSE)
void AVERAGE2COARSE(const int* mx,const int* my,const int *mz,
                    const int* mbc, const int* meqn,
					double qcoarse[],double qfine[],
					double areacoarse[],double areafine[],
					const int* igrid, const int* manifold);



#ifdef __cplusplus
}
#endif

#endif
