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

#include "bowl_user.h"
#include <fclaw2d_forestclaw.h>
#include <fclaw2d_clawpatch.h>
#include <fc2d_geoclaw.h>


static fclaw2d_vtable_t vt;
static fc2d_geoclaw_vtable_t geoclaw;

void bowl_link_solvers(fclaw2d_domain_t *domain)
{
    fclaw2d_init_vtable(&vt);

#if 0
    vt.problem_setup            = &fc2d_geoclaw_setprob;
    vt.patch_setup              = &bowl_patch_setup;
    vt.patch_initialize         = &bowl_patch_initialize;
    vt.patch_physical_bc        = &bowl_patch_physical_bc;
    vt.patch_single_step_update = &fc2d_geoclaw_update;  /* Includes b4step2 and src2 */

    vt.regrid_tag4refinement     = &bowl_patch_tag4refinement;
    vt.fort_tag4refinement      = &TAG4REFINEMENT;

    vt.regrid_tag4coarsening     = &bowl_patch_tag4coarsening;
    vt.fort_tag4coarsening      = &TAG4COARSENING;

    vt.write_header             = &fclaw2d_output_header_ascii;
    vt.fort_write_header        = &FCLAW2D_FORT_WRITE_HEADER;

    vt.patch_write_file         = &fclaw2d_output_patch_ascii;
    vt.fort_write_file          = &FCLAW2D_FORT_WRITE_FILE;
#endif

    fclaw2d_set_vtable(domain,&vt);

#if 0
    /* Needed for the geoclaw package */
    geoclaw.qinit = &QINIT;
    geoclaw.bc2 = &BC2;
    geoclaw.setaux = &SETAUX;
    geoclaw.setprob = &SETPROB;
    geoclaw.b4step2 = &B4STEP2;
    geoclaw.rpn2 = &RPN2;
    geoclaw.rpt2 = &RPT2;
#endif

    fc2d_geoclaw_set_vtable(&geoclaw);

}


void bowl_patch_setup(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx)
{
    /* Dummy setup - use multiple libraries */
    fc2d_geoclaw_setaux(domain,this_patch,this_block_idx,this_patch_idx);
}

void bowl_patch_initialize(fclaw2d_domain_t *domain,
                            fclaw2d_patch_t *this_patch,
                            int this_block_idx,
                            int this_patch_idx)
{
    /* This is an example of how to call the initialization routines explicitly
       This routine can be replaced by setting the appropriate fclaw2d_vtable_t,
       entry above, or by calling fclaw2d_geoclaw_qinit(...) from here. */

    int mx,my,mbc,meqn, maux;
    double xlower,ylower,dx,dy;
    double *q, *aux;

    vt = fclaw2d_get_vtable(domain);

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);
    fc2d_geoclaw_aux_data(domain,this_patch,&aux,&maux);

    /* Call to used defined, classic Clawpack (ver. 4.6)  'qinit' routine.
       Header is in the Clawpack package
    */
#if 0
    QINIT(&meqn,&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,q,&maux,aux);
#endif
}



void bowl_patch_physical_bc(fclaw2d_domain *domain,
                             fclaw2d_patch_t *this_patch,
                             int this_block_idx,
                             int this_patch_idx,
                             double t,
                             double dt,
                             fclaw_bool intersects_bc[],
                             fclaw_bool time_interp)
{
    /* This calls bc2 in bowl/user_4.6;  that file isn't changed but
       is included to show that both the local version of bc2.f and the
       geoclaw library code can be included */
    fc2d_geoclaw_bc2(domain,this_patch,this_block_idx,this_patch_idx,
                     t,dt,intersects_bc,time_interp);
}


int bowl_patch_tag4refinement(fclaw2d_domain_t *domain,
                               fclaw2d_patch_t *this_patch,
                               int blockno, int this_patch_idx,
                               int initflag)
{
    fclaw2d_vtable_t vt;
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *q;
    int tag_patch;
    const amr_options_t *amropt;
    double rt;

    amropt = get_domain_parms(domain);
    rt = amropt->refine_threshold;

    vt = fclaw2d_get_vtable(domain);

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);

    tag_patch = 0;
    vt.fort_tag4refinement(&mx,&my,&mbc,&meqn,&xlower,&ylower,
                           &dx,&dy,&blockno, q,&rt,&initflag,
                           &tag_patch);
    return tag_patch;
}

int bowl_patch_tag4coarsening(fclaw2d_domain_t *domain,
                               fclaw2d_patch_t *fine_patches,
                               int blockno, int patchno)

{
    fclaw2d_vtable_t vt;

    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *q[4];
    int tag_patch,igrid;
    double coarsen_threshold;
    fclaw2d_patch_t *patch0;

    patch0 = &fine_patches[0];

    const amr_options_t *amropt;
    amropt = get_domain_parms(domain);

    coarsen_threshold = amropt->coarsen_threshold;

    vt = fclaw2d_get_vtable(domain);

    fclaw2d_clawpatch_grid_data(domain,patch0,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    for (igrid = 0; igrid < 4; igrid++)
    {
        fclaw2d_clawpatch_soln_data(domain,&fine_patches[igrid],&q[igrid],&meqn);
    }
    tag_patch = 0;
    vt.fort_tag4coarsening(&mx,&my,&mbc,&meqn,&xlower,&ylower,&dx,&dy,
                           &blockno, q[0],q[1],q[2],q[3],
                           &coarsen_threshold, &tag_patch);
    return tag_patch;

}
