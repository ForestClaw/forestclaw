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

#include <fclaw2d_global.h>

#include <fclaw_package.h>
#include <fclaw_timer.h>

#include <fclaw2d_domain.h>
#include <fclaw2d_diagnostics.h>
#include <fclaw2d_map.h>

fclaw2d_global_t* fclaw2d_global_new (void)
{
    fclaw2d_global_t *glob;

    glob = FCLAW_ALLOC (fclaw2d_global_t, 1);
    glob->pkg_container = fclaw_package_container_new ();

    glob->count_amr_advance = 0;
    glob->count_ghost_exchange = 0;
    glob->count_amr_regrid = 0;
    glob->count_amr_new_domain = 0;
    glob->count_multiproc_corner = 0;
    glob->count_grids_per_proc = 0;
    glob->count_grids_remote_boundary = 0;
    glob->count_grids_local_boundary = 0;
    glob->count_single_step = 0;
    glob->count_elliptic_grids = 0;
    glob->curr_time = 0;
    glob->cont = NULL;
    glob->acc = FCLAW_ALLOC(fclaw2d_diagnostics_accumulator_t, 1);

    return glob;
}
void
fclaw2d_global_store_domain (fclaw2d_global_t* glob, fclaw2d_domain_t* domain)
{
    glob->domain = domain;
    glob->mpicomm = domain->mpicomm;
    glob->mpisize = domain->mpisize;
    glob->mpirank = domain->mpirank;
    
    glob->cont = (fclaw2d_map_context_t*)
           fclaw2d_domain_attribute_access (glob->domain, "fclaw_map_context", NULL);
}

void
fclaw2d_global_destroy (fclaw2d_global_t * glob)
{
    FCLAW_ASSERT (glob != NULL);

    fclaw_package_container_destroy ((fclaw_package_container_t *)glob->pkg_container);
    FCLAW_FREE (glob->acc);
    FCLAW_FREE (glob);
}

void fclaw2d_global_iterate_level (fclaw2d_global_t * glob, int level,
                                   fclaw2d_patch_callback_t pcb, void *user)
{
    fclaw2d_global_iterate_t g;
    g.glob = glob;
    g.user = user;
    fclaw2d_domain_iterate_level(glob->domain, level, pcb, &g);
}

void fclaw2d_global_iterate_patches (fclaw2d_global_t * glob,
                                     fclaw2d_patch_callback_t pcb, void *user)
{
    fclaw2d_global_iterate_t g;
    g.glob = glob;
    g.user = user;
    fclaw2d_domain_iterate_patches(glob->domain, pcb, &g);
}

void fclaw2d_global_iterate_families (fclaw2d_global_t * glob,
                                      fclaw2d_patch_callback_t pcb, void *user)
{
    fclaw2d_global_iterate_t g;
    g.glob = glob;
    g.user = user;
    fclaw2d_domain_iterate_families(glob->domain, pcb, &g);
}

void fclaw2d_global_iterate_adapted (fclaw2d_global_t * glob, fclaw2d_domain_t* new_domain,
                                     fclaw2d_match_callback_t mcb, void *user)
{
    fclaw2d_global_iterate_t g;
    g.glob = glob;
    g.user = user;
    fclaw2d_domain_iterate_adapted(glob->domain, new_domain,mcb,&g);
}

void fclaw2d_global_iterate_level_mthread (fclaw2d_global_t * glob, int level,
                                           fclaw2d_patch_callback_t pcb, void *user)
{
    fclaw2d_global_iterate_t g;
    g.glob = glob;
    g.user = user;
    fclaw2d_domain_iterate_level_mthread(glob->domain, level,pcb,&g);
}

void fclaw2d_global_iterate_partitioned (fclaw2d_global_t * glob,
                                         fclaw2d_domain_t * new_domain,
                                         fclaw2d_transfer_callback_t tcb,
                                         void *user)
{
    fclaw2d_global_iterate_t g;
    g.glob = glob;
    g.user = user;
    fclaw2d_domain_iterate_partitioned (glob->domain,new_domain,tcb,&g);
}
