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

#include <fclaw2d_forestclaw.h>
#include <fclaw2d_clawpatch_output.h>


static 
void cb_clawpatch_output (fclaw2d_domain_t * domain,
                          fclaw2d_patch_t * this_patch,
                          int this_block_idx, int this_patch_idx,
                          void *user)
{
    fclaw2d_global_iterate_t* s = (fclaw2d_global_iterate_t*) user;
    fclaw2d_global_t *glob = (fclaw2d_global_t*) s->glob;
    fclaw2d_clawpatch_vt *clawpatch_vt = fclaw2d_clawpatch_vt();

    int64_t patch_num;
    int level;
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *q;
    char fname[11];

    int iframe = *((int *) s->user);

    /* Get info not readily available to user */
    fclaw2d_patch_get_info(glob->domain,this_patch,this_block,
                           this_patch_idx,this_block_idx,
                           &patch_num,&level);

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);

    sprintf(fname,"fort.q%04d",iframe);

    /* The fort routine is defined by a clawpack solver and handles 
       the layout of q in memory (i,j,m) or (m,i,j), etc */
    clawpatch_vt->fort_write_file(fname,&mx,&my,&meqn,&mbc,&xlower,&ylower,&dx,&dy,q,
                                  &patch_num,&level,&this_block_idx,
                                  &glob->domain->mpirank);
}


/* this is a separate file, neeeded for each time frame */
void fclaw2d_clawpatch_ascii_header(fclaw2d_global_t* glob, int iframe)
{
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
    fclaw2d_clawpatch_vt *clawpatch_vt = fclaw2d_clawpatch_vt();
    int meqn,ngrids;
    char matname1[11];
    char matname2[11];

    sprintf(matname1,"fort.q%04d",iframe);
    sprintf(matname2,"fort.t%04d",iframe);

    double time = glob->curr_time;

    ngrids = fclaw2d_domain_get_num_patches(glob->domain);

    meqn = clawpatch_opt->meqn;

    fclaw2d_clawpatch_vt()->fort_write_header(matname1,matname2,&time,&meqn,&ngrids);
}


void fclaw2d_clawpatch_output_ascii(fclaw2d_global_t* glob,int iframe);
{
    fclaw2d_domain_t *domain = glob->domain;
    fclaw2d_clawpatch_vt *clawpatch_vt = fclaw2d_clawpatch_vt();

    /* BEGIN NON-SCALABLE CODE */
    /* Write the file contents in serial.
       Use only for small numbers of processors. */
    fclaw2d_domain_serialization_enter (domain);

    if (domain->mpirank == 0)
    {
        fclaw2d_clawpatch_ascii_header(glob,iframe);
    }

    /* Write out each patch to fort.qXXXX */
    fclaw2d_global_iterate_patches (glob, clawpatch_vt->cb_write_ascii_file, &iframe);

    fclaw2d_domain_serialization_leave (domain);
    /* END OF NON-SCALABLE CODE */
}

