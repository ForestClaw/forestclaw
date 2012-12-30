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

#include "amr_forestclaw.H"
#include "amr_waveprop.H"
#include "swirl_user.H"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

void swirl_link_solvers(fclaw2d_domain_t *domain)
{
    fclaw2d_solver_functions_t* sf = get_solver_functions(domain);

    sf->f_patch_setup              = &swirl_patch_setup;
    sf->f_patch_initialize         = &swirl_patch_initialize;
    sf->f_patch_physical_bc        = &swirl_patch_physical_bc;
    sf->f_patch_single_step_update = &swirl_patch_single_step_update;

    amr_waveprop_link_to_clawpatch();
}

void swirl_problem_setup(fclaw2d_domain_t* domain)
{
    amr_waveprop_setprob(domain);
}


void swirl_patch_setup(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx)
{
    /* Set velocity data */
    amr_waveprop_setaux(domain,this_patch,this_block_idx,this_patch_idx);

    /* Set up diffusion coefficients? Read in velocity data? Material properties? */
}



void swirl_patch_initialize(fclaw2d_domain_t *domain,
                            fclaw2d_patch_t *this_patch,
                            int this_block_idx,
                            int this_patch_idx)
{
    amr_waveprop_qinit(domain,this_patch,this_block_idx,this_patch_idx);
}


void swirl_patch_physical_bc(fclaw2d_domain *domain,
                             fclaw2d_patch_t *this_patch,
                             int this_block_idx,
                             int this_patch_idx,
                             double t,
                             double dt,
                             fclaw_bool intersects_bc[])
{
    amr_waveprop_bc2(domain,this_patch,this_block_idx,this_patch_idx,
                     t,dt,intersects_bc);
}


double swirl_patch_single_step_update(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int this_block_idx,
                                      int this_patch_idx,
                                      double t,
                                      double dt)
{
    /* Should I call b4step2 in here, or make this another function? */
    amr_waveprop_b4step2(domain,this_patch,this_block_idx,this_patch_idx,t,dt);

    double maxcfl = amr_waveprop_step2(domain,this_patch,this_block_idx,
                                       this_patch_idx,t,dt);
    return maxcfl;
}


#ifdef __cplusplus
#if 0
{
#endif
}
#endif
