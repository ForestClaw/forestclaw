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

#ifndef SWIRL_USER_H
#define SWIRL_USER_H

#include <amr_forestclaw.H>
#include <forestclaw2d.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


void swirl_write_tfile_(const int* iframe,const double* time,
                        const int* mfields,const int* ngrids);

void swirl_write_qfile_(const int* meqn,
                        const int* mbc, const int* mx,const int* my,
                        const double* xlower, const double* ylower,
                        const double* dx, const double* dy,
                        double q[], const int* iframe,
                        const int* patch_num,const int* level,
                        const int* blockno, const int* mpirank);

void qinit_manifold_(const int& meqn, const int& mbc,
                     const int& mx, const int& my,
                     const double& xlower, const double& ylower,
                     const double& dx, const double& dy,
                     const int& blockno,double q[],
                     const int& maux, double aux[]);

#define SETAUX_MANIFOLD FCLAW_F77_FUNC(setaux_manifold, SETAUX_MANIFOLD)

void SETAUX_MANIFOLD(const int& mbc,const int& mx, const int& my,
                     const double& xlower, const double& ylower,
                     const double& dx, const double& dy,
                     const int& maux, double aux[],
                     const int& blockno,
                     double xd[], double yd[], double zd[],
                     double area[]);

#define B4STEP2_MANIFOLD FCLAW_F77_FUNC(b4step2_manifold,B4STEP2_MANIFOLD)
void B4STEP2_MANIFOLD(const int& mbc,const int& mx,const int& my,
                      const double& dx,const double& dy,
                      const int& blockno,
                      double xd[], double yd[], double zd[],
                      const double& time,
                      const double& dt, const int& maux, double aux[]);

void swirl_link_solvers(fclaw2d_domain_t *domain);

void swirl_patch_setup(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx);

double swirl_patch_single_step_update(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int this_block_idx,
                                      int this_patch_idx,
                                      double t,
                                      double dt);

fclaw_bool swirl_patch_tag4refinement(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int this_block_idx, int this_patch_idx,
                                      int initflag);

fclaw_bool swirl_patch_tag4coarsening(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int blockno,
                                      int patchno);

void swirl_patch_b4step2(fclaw2d_domain_t *domain,
                         fclaw2d_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx,
                         double t,
                         double dt);

/* Mappings */
fclaw2d_map_context_t* fclaw2d_map_new_identity();

fclaw2d_map_context_t* fclaw2d_map_new_cart(fclaw2d_map_context_t* brick,
                                            const double scale[],
                                            const double shift[],
                                            const double rotate[]);

fclaw2d_map_context_t* fclaw2d_map_new_fivepatch(const double scale[],
                                                 const double shift[],
                                                 const double rotate[],
                                                 const double alpha);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
