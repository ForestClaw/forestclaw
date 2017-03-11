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

#include <fclaw2d_clawpatch.h>

const int SpaceDim = FCLAW2D_SPACEDIM;
const int NumFaces = FCLAW2D_NUMFACES;
const int p4est_refineFactor = FCLAW2D_P4EST_REFINE_FACTOR;
const int NumCorners = FCLAW2D_NUM_CORNERS;
const int NumSiblings = FCLAW2D_NUM_SIBLINGS;

fclaw2d_global_t *
fclaw2d_global_new (fclaw_options_t * gparms, fclaw2d_domain_t* domain)
{
    fclaw2d_global_t *glob;

    glob = FCLAW_ALLOC (fclaw2d_global_t, 1);
    if (gparms == NULL)
    {
        glob->gparms_owned = 1;
        glob->gparms = FCLAW_ALLOC_ZERO (fclaw_options_t, 1);
    }
    else
    {
        glob->gparms_owned = 0;
        glob->gparms = gparms;
    }

    glob->pkgs = fclaw_package_container_new ();

    fclaw2d_clawpatch_link_global (glob);

    glob->domain = domain;
    glob->mpicomm = domain->mpicomm;
    glob->mpisize = domain->mpisize;
    glob->mpirank = domain->mpirank;

    return glob;
}

void
fclaw2d_global_destroy (fclaw2d_global_t * glob)
{
    FCLAW_ASSERT (glob != NULL);

    fclaw_package_container_destroy (glob->pkgs);

    if (glob->gparms_owned)
    {
        FCLAW_FREE (glob->gparms);
    }
    FCLAW_FREE (glob);
}

fclaw_package_container_t *
fclaw2d_global_get_container (fclaw2d_global_t * glob)
{
    FCLAW_ASSERT (glob != NULL);
    FCLAW_ASSERT (glob->pkgs != NULL);

    return glob->pkgs;
}

fclaw_app_t* fclaw2d_global_get_app(fclaw2d_global_t* glob)
{
    fclaw_app_t* app;
    app = fclaw2d_domain_get_app(glob->domain);

    return app;
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
    fclaw2d_domain_iterate_adapted(glob->domain, new_domain,
                                   mcb,
                                   &g);
}

void fclaw2d_global_iterate_level_mthread (fclaw2d_global_t * glob, int level,
                                           fclaw2d_patch_callback_t pcb, void *user)
{
    fclaw2d_global_iterate_t g;
    g.glob = glob;
    g.user = user;
    fclaw2d_domain_iterate_level_mthread(glob->domain, level,
                                         pcb, 
                                         &g);
}