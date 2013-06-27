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

#ifndef RIDGE_USER_H
#define RIDGE_USER_H

#include "amr_forestclaw.H"
#include "amr_waveprop.H"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

extern "C"
{

    void setaux_ridge_(const int& maxmx, const int& maxmy, const int& mbc,
                       const int& mx, const int& my,
                       const double& xlower, const double& ylower,
                       const double& dx, const double& dy,
                       const int& maux, double aux[],
                       double xnormals[], double ynormals[],
                       double xtangents[], double ytangents[],
                       double surfnormals[], const int& level,
                       const int& maxlevel, const int& refratio);

    void write_qfile_geo_(const int& maxmx, const int& maxmy, const int& meqn, const int& mbc,
                          const int& mx, const int& my,
                          const double& xlower, const double& ylower,
                          const double& dx, const double& dy, double q[],
                          const int& iframe,
                          const int& patch_idx,const int& level,const int& blockno,
                          const int& maux, double aux[]);

    void ridge_tag4refinement_(const int& mx,const int& my,const int& mbc,
                               const int& meqn, const double& xlower,
                               const double& ylower, const double& dx, const double& dy,
                               double q[], const int& level, const int& maux,
                               double aux[], const int& initflag, int& tag_patch);

    void ridge_tag4coarsening_(const int& mx,const int& my,const int& mbc,
                               const int& meqn,const double& xlower,
                               const double& ylower, const double& dx, const double& dy,
                               double qcoarsened[], const int& maux,
                               double aux[], int& tag_patch);

    void interpolate_geo_(const int& mx, const int& my, const int& mbc,
                          const int& meqn, double qcoarse[], double qfine[],
                          const int& maux, double aux_coarse[], double aux_fine[],
                          const int& p4est_refineFactor, const int& refratio,
                          const int& igrid);

}

void ridge_link_solvers(fclaw2d_domain_t *domain);
void ridge_link_regrid_functions(fclaw2d_domain_t* domain);


void ridge_problem_setup(fclaw2d_domain_t* domain);

void ridge_patch_setup(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx);

void ridge_patch_initialize(fclaw2d_domain_t *domain,
                            fclaw2d_patch_t *this_patch,
                            int this_block_idx,
                            int this_patch_idx);

void ridge_patch_physical_bc(fclaw2d_domain *domain,
                        fclaw2d_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx,
                        double t,
                        double dt,
                        fclaw_bool intersects_bc[]);

double ridge_patch_single_step_update(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int this_block_idx,
                                      int this_patch_idx,
                                      double t,
                                      double dt);

void ridge_patch_output(fclaw2d_domain_t *domain, fclaw2d_patch_t *this_patch,
                        int this_block_idx, int this_patch_idx,
                        int iframe,int num,int matlab_level);

fclaw_bool ridge_tag4coarsening(fclaw2d_domain_t *domain,
                                fclaw2d_patch_t *this_patch,
                                int this_blockno,
                                int this_patchno);

fclaw_bool ridge_tag4refinement(fclaw2d_domain_t *domain,
                                fclaw2d_patch_t *this_patch,
                                int this_block_idx, int this_patch_idx,
                                int initflag);

void ridge_interpolate2fine(fclaw2d_domain_t* domain, fclaw2d_patch_t *coarse_patch,
                            fclaw2d_patch_t* fine_patch,
                            int this_blockno, int coarse_patchno,
                            int fine_patchno, int igrid);



#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
