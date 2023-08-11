/*
Copyright (c) 2019-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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


#include "fc2d_thunderegg_physical_bc.h"

#include "fc2d_thunderegg.h"
#include "fc2d_thunderegg_options.h"

#include <fclaw2d_elliptic_solver.h>

#include <fclaw_global.h>
#include <fclaw2d_domain.h>
#include <fclaw2d_patch.h>
#include <fclaw_clawpatch.h>
#include <fclaw2d_physical_bc.h>

void cb_fc2d_thunderegg_physical_bc(fclaw_domain_t *domain,
                                   fclaw_patch_t *patch,
                                   int blockno,
                                   int patchno,
                                   void *user)

{
    fclaw_global_iterate_t* s = (fclaw_global_iterate_t*) user;
    fc2d_thunderegg_time_info_t *tinfo = (fc2d_thunderegg_time_info_t*) s->user;

    double t = tinfo->t;


    /* Determine which faces are at the physical boundary */
    int intersects_bc[4];
    fclaw2d_physical_get_bc(s->glob,blockno,patchno,intersects_bc);

    int mx, my, mbc;
    double xlower, ylower, dx, dy;
    fclaw2d_clawpatch_grid_data(s->glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    const fc2d_thunderegg_options_t* mg_opt = fc2d_thunderegg_get_options(s->glob);
    int mfields;
    double *rhs;
    fclaw_clawpatch_rhs_data(s->glob,patch,&rhs,&mfields);

    fc2d_thunderegg_vtable_t*  mg_vt = fc2d_thunderegg_vt(s->glob);

    int cons_check = 0;
    double flux_sum[4];

    mg_vt->fort_apply_bc(&blockno, &mx, &my, &mbc, &mfields, 
                         &xlower, &ylower, &dx,&dy,&t, intersects_bc,
                         mg_opt->boundary_conditions,rhs, mg_vt->fort_eval_bc,
                         &cons_check, flux_sum);

#if 0
    fclaw2d_patch_physical_bc(s->glob,
                              patch,
                              blockno,
                              patchno,
                              intersects_bc);
#endif                              
}

/* This is needed by other routines, so we don't set it to static. */
void fc2d_thunderegg_physical_get_bc(fclaw_global_t *glob,
                                    int blockno,
                                    int patchno,
                                    int *intersects_bdry)
{
    // const int numfaces = get_faces_per_patch(domain);
    int bdry[4];
    fclaw2d_patch_boundary_type(glob->domain,blockno,patchno,bdry);
    int i;
    for(i = 0; i < 4; i++)
    {
        // Physical boundary conditions
        intersects_bdry[i] = bdry[i] == 1;
    }
}



/* -----------------------------------------------------------------------------
   Public interface : Set physical boundary conditions on a patch
   ----------------------------------------------------------------------------- */

void fc2d_thunderegg_physical_bc(fclaw_global_t *glob)
{

    fc2d_thunderegg_time_info_t tinfo;
    tinfo.t = glob->curr_time;
    fclaw_global_iterate_patches(glob,
                                   cb_fc2d_thunderegg_physical_bc,
                                   (void *) &tinfo);
}

