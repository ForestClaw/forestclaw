/*
Copyright (c) 2012-2023 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include "swirl/swirl_user.h"

void swirl_create_domain(fclaw2d_global_t *glob)
{
    fclaw2d_set_global_context(glob);    

    const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);

    fclaw2d_domain_t *domain = 
          fclaw2d_domain_new_unitsquare(glob->mpicomm, 
                                        fclaw_opt->minlevel);
    /* Create "empty" mapping */
    fclaw2d_map_context_t* cont = fclaw2d_map_new_nomap();

    /* Store domain in the glob */
    fclaw2d_global_store_domain(glob, domain);

    fclaw2d_global_store_map (glob, cont);

    fclaw2d_clear_global_context(glob);    
}

void swirl_initialize(fclaw2d_global_t* glob)
{
    fclaw2d_set_global_context(glob);

    /* ---------------------------------------------------------------
       Set domain data.
       --------------------------------------------------------------- */
    fclaw2d_domain_data_new(glob->domain);

    const swirl_options_t *swirl_opt = swirl_get_options(glob);

    /* Initialize virtual table for ForestClaw */
    fclaw2d_vtables_initialize(glob);

    /* Initialize virtual tables for solvers */
    if (swirl_opt->claw_version == 4)
    {
        fc2d_clawpack46_solver_initialize(glob);
    }
    else if (swirl_opt->claw_version == 5)
    {
        fc2d_clawpack5_solver_initialize(glob);
    }

    swirl_link_solvers(glob);

    /* ---------------------------------------------------------------
       Run
       --------------------------------------------------------------- */
    fclaw2d_initialize(glob);

    fclaw2d_clear_global_context(glob);
}


void swirl_finalize(fclaw2d_global_t* glob)
{
    fclaw2d_set_global_context(glob);

    fclaw2d_problem_setup(glob);
    fclaw2d_finalize(glob);

    fclaw2d_clear_global_context(glob);
}

