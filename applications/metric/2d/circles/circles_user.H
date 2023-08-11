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

#ifndef SPHERE_USER_H
#define SPHERE_USER_H

#include "amr_forestclaw.H"

extern "C"
{

    void setprob_();
    void qinit_manifold_(const int& maxmx, const int& maxmy, const int& meqn,
                         const int& mbc, const int& mx, const int& my,
                         const double& xlower, const double& ylower,
                         const double& dx, const double& dy,
                         double q[], const int& maux, double aux[],
                         double xp[], double yp[], double zp[]);

    void output_metric_(const int& maxmx, const int& maxmy, const int& meqn,
                        const int& mbc, const int& mx, const int& my,
                        const double& xlower, const double& ylower,
                        const double& dx, const double& dy,
                        double q[], const int& iframe, const int& patch_num,
                        const int& level, const int& blockno,
                        double xp[], double yp[]);

    void tag4refinement_(const int& mx, const int& my, const int& mbc,
                         const int& meqn,const double& xlower, const double& ylower,
                         const double& dx, const double& dy,
                         const int& init_flag,
                         const int& blockno, const int& tag_patch,
                         double xp[], double yp[]);

    void tag4coarsening_(const int& mx, const int& my, const int& mbc,
                         const int& meqn,const double& xlower, const double& ylower,
                         const double& dx, const double& dy,
                         double qcoarsened[], const int& tag_patch);
}

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

void sphere_link_patch(fclaw_domain_t *domain);

void sphere_setprob(fclaw_domain_t* domain);

void sphere_qinit(fclaw_domain_t *domain,
                      fclaw_patch_t *this_patch,
                      int this_block_idx,
                      int this_patch_idx);


int sphere_patch_tag4refinement(fclaw_domain_t *domain,
                                       fclaw_patch_t *this_patch,
                                       int this_block_idx, int this_patch_idx,
                                       int initflag);

int sphere_patch_tag4coarsening(fclaw_domain_t *domain,
                                       fclaw_patch_t *this_patch,
                                       int blockno,
                                       int patchno);

void metric_write_output(fclaw_domain_t *domain, fclaw_patch_t *this_patch,
                         int this_block_idx, int this_patch_idx,
                         int iframe,int num,int matlab_level);

void metric_write_header(fclaw_domain_t* domain, int iframe,int ngrids);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
