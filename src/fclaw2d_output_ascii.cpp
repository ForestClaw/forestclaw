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

/* #include <amr_utils.hpp> */
#include <fclaw2d_forestclaw.h>
#include <fclaw2d_output.h>

#include <fclaw2d_vtable.h>
#include <fclaw2d_clawpatch.hpp>


#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

void fclaw2d_output_header_ascii(fclaw2d_domain_t* domain,
                                 int iframe)
{
    const amr_options_t *amropt;
    fclaw2d_vtable_t vt;
    int meqn,ngrids;
    double time;

    amropt = fclaw2d_forestclaw_get_options(domain);

    time = fclaw2d_domain_get_time(domain);
    ngrids = fclaw2d_domain_get_num_patches(domain);

    meqn = amropt->meqn;

    vt = fclaw2d_get_vtable(domain);
    vt.fort_write_header(&iframe,&time,&meqn,&ngrids);

    /* Is this really necessary? */
    /* FCLAW2D_OUTPUT_NEW_QFILE(&iframe); */
}


void fclaw2d_output_patch_ascii(fclaw2d_domain_t *domain,
                               fclaw2d_patch_t *this_patch,
                               int this_block_idx, int this_patch_idx,
                               int iframe,int patch_num,int level)
{
    fclaw2d_vtable_t vt;
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *q;
    vt = fclaw2d_get_vtable(domain);

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);

    vt.fort_write_file(&mx,&my,&meqn,&mbc,&xlower,&ylower,&dx,&dy,q,
                       &iframe,&patch_num,&level,&this_block_idx,
                       &domain->mpirank);
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
