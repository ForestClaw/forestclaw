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

#ifndef NO_SOLVER_USER_H
#define NO_SOLVER_USER_H

#include "amr_forestclaw.H"

extern "C"
{
    void no_solver_setprob(fclaw2d_domain_t* domain);

    void initialize_(const int& mx, const int& my, const int& meqn,
                     const int& mbc,
                     const double& xlower, const double& ylower,
                     const double& dx, const double& dy,
                     double q[]);

    void initialize_ac_(const int& mx, const int& my, const int& meqn,
                        const int& mbc,
                        const double& xlower, const double& ylower,
                        const double& dx, const double& dy,
                        double q[]);

    void no_solver_tag4refinement_(const int& mx,const int& my,const int& mbc,
                                   const int& meqn, const double& xlower,
                                   const double& ylower, const double& dx, const double& dy,
                                   double q[], const int& init_flag, int& tag_patch);

    void no_solver_tag4coarsening_(const int& mx,const int& my,const int& mbc,
                                   const int& meqn,const double& xlower,
                                   const double& ylower, const double& dx, const double& dy,
                                   double qcoarsened[], int& tag_patch);

}

void no_solver_patch_initialize(fclaw2d_domain_t *domain,
                                fclaw2d_patch_t *this_patch,
                                int this_block_idx,
                                int this_patch_idx);

fclaw_bool no_solver_patch_tag4refinement(fclaw2d_domain_t *domain,
                                          fclaw2d_patch_t *this_patch,
                                          int this_block_idx, int this_patch_idx,
                                          int initflag);


fclaw_bool no_solver_patch_tag4coarsening(fclaw2d_domain_t *domain,
                                          fclaw2d_patch_t *this_patch,
                                          int blockno,
                                          int patchno);


#endif
