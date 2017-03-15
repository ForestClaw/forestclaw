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
#include <fclaw2d_output.h>
#include <fclaw2d_clawpatch.h>


void fclaw2d_clawpatch_output_ascii_header(fclaw2d_global_t* glob,
                                           int iframe)
{
    const amr_options_t *amropt = fclaw2d_forestclaw_get_options(glob);
    int meqn,ngrids;
    double time;
    char matname1[11];
    char matname2[11];

    sprintf(matname1,"fort.q%04d",iframe);
    sprintf(matname2,"fort.t%04d",iframe);

    time = glob->curr_time;
    ngrids = fclaw2d_domain_get_num_patches(glob->domain);

    meqn = amropt->meqn;

    fclaw2d_clawpatch_vt()->fort_write_header(matname1,matname2,&time,&meqn,&ngrids);
}


void fclaw2d_clawpatch_output_ascii(fclaw2d_global_t *glob,
                                    fclaw2d_patch_t *this_patch,
                                    int this_block_idx, int this_patch_idx,
                                    int iframe,int patch_num,int level)
{
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *q;
    char fname[11];

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);

    sprintf(fname,"fort.q%04d",iframe);
    fclaw2d_clawpatch_vt()->fort_write_file(fname,&mx,&my,&meqn,&mbc,&xlower,&ylower,&dx,&dy,q,
                                            &patch_num,&level,&this_block_idx,
                                            &glob->domain->mpirank);
}
