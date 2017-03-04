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

#include "slosh_user.h"
#include <fclaw2d_forestclaw.h>
#include <fclaw2d_clawpatch.h>
#include <fc2d_geoclaw.h>


// static fclaw2d_vtable_t vt;
// static fc2d_geoclaw_vtable_t geoclaw;

void slosh_link_solvers(fclaw2d_domain_t *domain)
{

    /* These are set by GeoClaw for convenience, but the user
       can set these with customized functions, if desired. */
    // fc2d_geoclaw_init_vtables(&vt, &geoclaw);
    // fc2d_geoclaw_init_vtables();
    
    fclaw2d_patch_vtable_t* patch_vt = fclaw2d_patch_vt();
    patch_vt->patch_initialize       = &slosh_patch_initialize;

    // fc2d_geoclaw_set_vtables(domain,&vt,&geoclaw);
    fclaw2d_set_vtable();

}

void slosh_patch_initialize(fclaw2d_domain_t *domain,
                            fclaw2d_patch_t *this_patch,
                            int this_block_idx,
                            int this_patch_idx)
{
    /* This is an example of how to call the initialization routines explicitly
       This routine can be replaced by setting the appropriate fclaw2d_vtable_t,
       entry above, or by calling fclaw2d_clawpack5_qinit(...) from here. */

    int mx,my,mbc,meqn, maux;
    double xlower,ylower,dx,dy;
    double *q, *aux;

    // int minlevel, maxlevel;
    // minlevel = 4;
    // maxlevel = 4;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);
    fc2d_geoclaw_aux_data(domain,this_patch,&aux,&maux);

    /* Call to used defined, classic Clawpack (ver. 4.6)  'qinit' routine.
       Header is in the Clawpack package
    */
    // PREQINIT(&minlevel,&maxlevel);
    QINIT(&meqn,&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,q,&maux,aux);
}
