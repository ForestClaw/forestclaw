/*
Copyright (c) 2012-2023 Carsten Burstedde, Donna Calhoun, Scott Aiton,
Hannes Brandt
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

/* This example demonstrates the use of fclaw2d_overlap_exchange to exchange
 * interpolation data between two meshes.
 * Both meshes have a clearly assigned role:
 *  consumer (e.g. Gemini) - queries data for points - represented by swirl
 *  producer (e.g. MAGIC) - provides data - represented by filament. */

#include "filament/filament_user.h"

void filament_create_domain(fclaw_global_t *glob)
{
    fclaw_set_global_context(glob);    

    const fclaw_options_t* fclaw_opt = fclaw_get_options(glob);

    int mi = fclaw_opt->mi;
    int mj = fclaw_opt->mj;
    int a = fclaw_opt->periodic_x;
    int b = fclaw_opt->periodic_y;

    /* Square brick domain */
    fclaw_domain_t *domain =
        fclaw_domain_new_2d_brick (glob->mpicomm, mi, mj, a, b,
                                      fclaw_opt->minlevel);
    fclaw2d_map_context_t* brick = 
        fclaw2d_map_new_brick (domain, mi, mj, a, b);
        
    /* Square in [-1,1]x[-1,1], shifted by (1,1,0) */
    fclaw2d_map_context_t *cont  = 
        fclaw2d_map_new_cart(brick,
                             fclaw_opt->scale,
                             fclaw_opt->shift);

    /* Store domain in the glob */
    fclaw_global_store_domain(glob, domain);

    fclaw2d_map_store (glob, cont);

    fclaw_clear_global_context(glob);    
}

void filament_initialize(fclaw_global_t* glob)
{
    fclaw_set_global_context(glob);

    const filament_options_t *user = filament_get_options(glob);

    /* Initialize virtual table for ForestClaw */
    fclaw_vtables_initialize(glob);

    if (user->claw_version == 4)
    {
      fc2d_clawpack46_solver_initialize(glob);
    }
    else if (user->claw_version == 5)
    {
      fc2d_clawpack5_solver_initialize(glob);
    }

    filament_link_solvers(glob);
    fclaw_initialize(glob);

    fclaw_clear_global_context(glob);
}

void filament_finalize(fclaw_global_t* glob)
{
    fclaw_set_global_context(glob);

    fclaw_problem_setup(glob);
    fclaw_finalize(glob);

    fclaw_clear_global_context(glob);
}

