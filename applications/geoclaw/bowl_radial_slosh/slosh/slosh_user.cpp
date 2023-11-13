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

#include "slosh_user.h"
#include <fclaw_filesystem.h>

void slosh_link_solvers(fclaw_global_t *glob)
{
    fc2d_geoclaw_vtable_t* geoclaw_vt = fc2d_geoclaw_vt(glob);
    geoclaw_vt->qinit       = &FC2D_GEOCLAW_QINIT;
}

void slosh_create_domain(fclaw_global_t* glob)
{
    fclaw_options_t *fclaw_opts = fclaw_get_options(glob);

    /* Size is set by [ax,bx] x [ay, by], set in .ini file */
    fclaw_domain_t *domain = 
        fclaw_domain_new_unitsquare(glob->mpicomm, fclaw_opts->minlevel);
    fclaw_map_context_t* cont = fclaw2d_map_new_nomap();

    /* store domain and map in glob */
    fclaw_global_store_domain(glob, domain);
    fclaw_map_store(glob, cont);

    fclaw_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);
}

void slosh_run_program(fclaw_global_t* glob)
{
    fclaw_set_global_context(glob);

    /* Initialize virtual table for ForestClaw */
    fclaw_vtables_initialize(glob);

    fc2d_geoclaw_solver_initialize(glob);

    slosh_link_solvers(glob);

    /* ---------------------------------------------------------------
       Run
       --------------------------------------------------------------- */
    fc2d_geoclaw_module_setup(glob);


    fclaw_initialize(glob);
    fc2d_geoclaw_run(glob);

    fclaw_finalize(glob);

    fclaw_clear_global_context(glob);
}