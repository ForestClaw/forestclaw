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

#include "wavetank_user.h"

/* 
#include "fc2d_geoclaw_output_ascii.h"

#include "fc2d_geoclaw.h"
#include "fc2d_geoclaw_fort.h"
#include "fc2d_geoclaw_options.h"
*/

#include <fclaw2d_clawpatch.h>  /* Include patch, domain declarations */
#include <fclaw2d_clawpatch_options.h>  /* Include patch, domain declarations */

#include <fclaw2d_patch.h>
#include <fclaw2d_global.h>

static
void cb_geoclaw_output_ascii(fclaw2d_domain_t *domain,
                             fclaw2d_patch_t *patch,
                             int blockno, int patchno,
                             void *user)
{
    fclaw2d_global_iterate_t* s = (fclaw2d_global_iterate_t*) user;
    fclaw2d_global_t *glob = (fclaw2d_global_t*) s->glob;

    int iframe = *((int *) s->user);

    /* Get info not readily available to user */
    int patch_num, level, global_num;
    fclaw2d_patch_get_info(glob->domain,patch,
                           blockno,patchno,
                           &global_num,
                           &patch_num,&level);

    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *q;
    int meqn;
    fclaw2d_clawpatch_soln_data(glob,patch,&q,&meqn);

    double *aux;
    int maux;
    fclaw2d_clawpatch_aux_data(glob,patch,&aux,&maux);

    int mbathy = 1;
    FC2D_GEOCLAW_FORT_WRITE_FILE(&mx,&my,&meqn,&maux,&mbathy,&mbc,&xlower,&ylower,
                                 &dx,&dy,q,aux,&iframe,&patch_num,&level,
                                 &blockno,&glob->mpirank);
}

static
void geoclaw_header_ascii(fclaw2d_global_t* glob,int iframe)
{
    const fclaw2d_clawpatch_options_t *clawpatch_opt;

    int meqn,maux,ngrids;
    double time;

    clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    time = glob->curr_time;
    ngrids = glob->domain->global_num_patches;

    meqn = clawpatch_opt->meqn;
    maux = clawpatch_opt->maux;

    FC2D_GEOCLAW_FORT_WRITE_HEADER(&iframe,&time,&meqn,&maux,&ngrids);

}

/* --------------------------------------------------------------
	Public interface
   ------------------------------------------------------------ */

void fc2d_geoclaw_output_ascii(fclaw2d_global_t* glob,int iframe)
{
    fclaw2d_domain_t *domain = glob->domain;

    /* BEGIN NON-SCALABLE CODE */
    /* Write the file contents in serial.
       Use only for small numbers of processors. */
    fclaw2d_domain_serialization_enter (domain);

    if (domain->mpirank == 0)
    {
        geoclaw_header_ascii(glob,iframe);
    }

    /* Write out each patch to fort.qXXXX */
    fclaw2d_global_iterate_patches (glob, cb_geoclaw_output_ascii, &iframe);

    fclaw2d_domain_serialization_leave (domain);
    /* END OF NON-SCALABLE CODE */
}





