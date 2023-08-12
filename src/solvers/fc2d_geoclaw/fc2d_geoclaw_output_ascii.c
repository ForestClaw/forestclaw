/*
Copyright (c) 2012 Carsten Burstedde, Donna Calhoun, Yu-Hsuan Shih
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

#include "fc2d_geoclaw_output_ascii.h"

#include "fc2d_geoclaw.h"
#include "fc2d_geoclaw_fort.h"
#include "fc2d_geoclaw_options.h"

#include <fclaw_clawpatch.h>  /* Include patch, domain declarations */
#include <fclaw_clawpatch_options.h>  /* Include patch, domain declarations */

#include <fclaw_patch.h>
#include <fclaw_global.h>

static
void cb_geoclaw_output_ascii(fclaw_domain_t *domain,
                             fclaw_patch_t *patch,
                             int blockno, int patchno,
                             void *user)
{
    fclaw_global_iterate_t* s = (fclaw_global_iterate_t*) user;
    fclaw_global_t *glob = (fclaw_global_t*) s->glob;

    int iframe = *((int *) s->user);    

    /* Get info not readily available to user */
    int local_num, global_num, level;
    fclaw_patch_get_info(glob->domain,patch,
                           blockno,patchno,
                           &global_num, 
                           &local_num,&level);

    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *q;
    int meqn;
    fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

    double *aux;
    int maux;
    fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

    FC2D_GEOCLAW_FORT_WRITE_FILE(&mx,&my,&meqn, &maux,&mbc,&xlower,&ylower,
                                 &dx,&dy,q,aux,&iframe,&global_num,&level,
                                 &blockno,&glob->mpirank);
}

static
void geoclaw_header_ascii(fclaw_global_t* glob,int iframe)
{
    double time = glob->curr_time;
    int ngrids = glob->domain->global_num_patches;

    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
    int meqn = clawpatch_opt->meqn;
    int maux = clawpatch_opt->maux;

    FC2D_GEOCLAW_FORT_WRITE_HEADER(&iframe,&time,&meqn,&maux,&ngrids);
}

/* --------------------------------------------------------------
	Public interface
   ------------------------------------------------------------ */

void fc2d_geoclaw_output_ascii(fclaw_global_t* glob,int iframe)
{
    fclaw_domain_t *domain = glob->domain;

    /* BEGIN NON-SCALABLE CODE */
    /* Write the file contents in serial.
       Use only for small numbers of processors. */
    fclaw_domain_serialization_enter (domain);

    if (domain->mpirank == 0)
        geoclaw_header_ascii(glob,iframe);

    /* Write out each patch to fort.qXXXX */
    fclaw_global_iterate_patches (glob, cb_geoclaw_output_ascii, &iframe);

    fclaw_domain_serialization_leave (domain);
    /* END OF NON-SCALABLE CODE */
}





