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

#ifndef REFINE_DIM
#define REFINE_DIM 2
#endif

#ifndef PATCH_DIM
#define PATCH_DIM 2
#endif

#if REFINE_DIM == 2 && PATCH_DIM == 2

#include <fclaw2d_clawpatch_output_ascii.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>

#elif REFINE_DIM ==2 && PATCH_DIM == 3

#include <fclaw3dx_clawpatch_output_ascii.h>

#include <fclaw3dx_clawpatch.h>
#include <fclaw3dx_clawpatch_options.h>

#include <_fclaw2d_to_fclaw3dx.h>

#endif
#include <fclaw2d_patch.h>
#include <fclaw2d_global.h>
#include <fclaw2d_options.h>


void cb_clawpatch_output_ascii (fclaw2d_domain_t * domain,
                                fclaw2d_patch_t * patch,
                                int blockno, int patchno,
                                void *user)
{
    fclaw2d_global_iterate_t* s = (fclaw2d_global_iterate_t*) user;
    fclaw2d_global_t *glob = (fclaw2d_global_t*) s->glob;


    /* Get info not readily available to user */
    int global_num, local_num, level;
    fclaw2d_patch_get_info(glob->domain,patch,
                           blockno,patchno,
                           &global_num,&local_num, &level);
    
    int meqn;
    double *q;
    fclaw2d_clawpatch_soln_data(glob,patch,&q,&meqn);

    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    int iframe = *((int *) s->user);
    char fname[BUFSIZ];
    snprintf (fname, BUFSIZ, "%s.q%04d", fclaw_opt->prefix, iframe);
    /* sprintf(fname,"fort.q%04d",iframe); */

    /* The fort routine is defined by a clawpack solver and handles 
       the layout of q in memory (i,j,m) or (m,i,j), etc */
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt(glob);
    FCLAW_ASSERT(clawpatch_vt->fort_output_ascii);

    int mx,my,mbc;
    double xlower,ylower,dx,dy;
#if PATCH_DIM == 2
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);
    clawpatch_vt->fort_output_ascii(fname,&mx,&my,&meqn,&mbc,
                                    &xlower,&ylower,&dx,&dy,q,
                                    &global_num,&level,&blockno,
                                    &glob->mpirank);
#else
    int mz;
    double zlower, dz;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mz,&mbc,
                                 &xlower,&ylower,&zlower,
                                 &dx,&dy,&dz);
    clawpatch_vt->fort_output_ascii(fname,&mx,&my,&mz,&meqn,&mbc,
                                    &xlower,&ylower,&zlower,
                                    &dx,&dy,&dz,q,
                                    &global_num,&level,&blockno,
                                    &glob->mpirank);
#endif

}


/* This function isn't virtualized;  should it be? */
void fclaw2d_clawpatch_time_header_ascii(fclaw2d_global_t* glob, int iframe)
{
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt(glob);
    char matname1[11];
    char matname2[11];

    sprintf(matname1,"fort.q%04d",iframe);
    sprintf(matname2,"fort.t%04d",iframe);

    double time = glob->curr_time;

    int ngrids = glob->domain->global_num_patches;

    int meqn = clawpatch_opt->meqn;
    int maux = clawpatch_opt->maux;

    clawpatch_vt->fort_header_ascii(matname1,matname2,&time,&meqn,&maux,&ngrids);
}

    /*--------------------------------------------------------------------
    Public interface
    Use this function as follows : 
           fclaw2d_vtable_t *vt = fclaw2d_vt(glob);
           vt->output_frame = &fclaw2d_clawpatch_output_ascii;
    -------------------------------------------------------------------- */

void fclaw2d_clawpatch_output_ascii(fclaw2d_global_t* glob,int iframe)
{
    fclaw2d_domain_t *domain = glob->domain;
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt(glob);

    /* BEGIN NON-SCALABLE CODE */
    /* Write the file contents in serial.
       Use only for small numbers of processors. */
    fclaw2d_domain_serialization_enter (domain);

    if (glob->mpirank == 0)
    {
        clawpatch_vt->time_header_ascii(glob,iframe);
    }

    /* Write out each patch to fort.qXXXX */
    fclaw2d_global_iterate_patches (glob, clawpatch_vt->cb_output_ascii, &iframe);

    fclaw2d_domain_serialization_leave (domain);
    /* END OF NON-SCALABLE CODE */
}

