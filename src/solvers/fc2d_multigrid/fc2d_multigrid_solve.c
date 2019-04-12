/*
Copyright (c) 2019 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#include "fc2d_multigrid.h"
#include "fc2d_multigrid_options.h"

#include <fclaw2d_elliptic_solver.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_clawpatch_output_ascii.h> 
#include <fclaw2d_clawpatch_output_vtk.h>

#include <fclaw2d_patch.h>
#include <fclaw2d_global.h>
#include <fclaw2d_vtable.h>

#include <p4est_bits.h>
#include <p4est_wrap.h>

typedef struct mg_copy_info
{
    double *sdata;
    int copy_mode;    /* 0 = copy to sdata;  1 = copy from sdata */
} mg_copy_info_t;

void cb_multigrid_copy_data(fclaw2d_domain_t *domain,
                            fclaw2d_patch_t *patch,
                            int blockno,
                            int patchno,
                            void* user)
{
    fclaw2d_global_iterate_t* g = (fclaw2d_global_iterate_t*) user;

    int mx,my,mbc,meqn;

    int global_num, level;
    double *q; 

    fclaw2d_clawpatch_options_t *clawpatch_opt;

    clawpatch_opt = fclaw2d_clawpatch_get_options(g->glob);
    mx = clawpatch_opt->mx;
    my = clawpatch_opt->my;
    mbc = clawpatch_opt->mbc;

    mg_copy_info_t *cpinfo = (mg_copy_info_t*) g->user; 

    /* We assume that q stores right hand side. This is fine when solving a single 
       elliptic problem, but more generally, we should store the rhs in some other 
       location */
    fclaw2d_clawpatch_soln_data(g->glob,patch,&q,&meqn);
    FCLAW_ASSERT(meqn == 1);

    fclaw2d_patch_get_info(domain,patch, blockno, patchno,
                           &global_num, &level);

    size_t size = (mx+2*mbc)*(my+2*mbc);

    if (cpinfo->copy_mode == 0)
    {
        /* Copy to cpinfo->sdata; */
        memcpy(&cpinfo->sdata[global_num*size],q,size);        
    }
    else
    {
        /* Copy from cpinfo->sdata; */
        memcpy(q,&cpinfo->sdata[global_num*size],size);                
    }

}

static 
void multigrid_copy_data(fclaw2d_global_t *glob,double *solver_data,int copy_mode)
{
    /* Copy patch data to new memory location */
    mg_copy_info_t cpinfo;
    cpinfo.sdata = solver_data;
    cpinfo.copy_mode = copy_mode;

    fclaw2d_global_iterate_patches (glob, cb_multigrid_copy_data, &cpinfo);
}


void fc2d_multigrid_solve(fclaw2d_global_t *glob)
{
    fclaw2d_domain_t *domain = glob->domain;

    int mx, my, mbc;

    fclaw2d_clawpatch_options_t *clawpatch_opt;

    clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
    mx = clawpatch_opt->mx;
    my = clawpatch_opt->my;
    mbc = clawpatch_opt->mbc;

    size_t size = (mx+2*mbc)*(my+2*mbc)*domain->local_num_patches;

    double *qsoln = FCLAW_ALLOC(double, size);

    /* Copy data to location that can be used by solver. */
    int copy_mode = 0;
    multigrid_copy_data(glob,qsoln,copy_mode);

#if 0
    /* A start ... */

    p4est_wrap_t *wrap = (p4est_wrap_t *) domain->pp;
    p4est_t *p4est = wrap->p4est;
    p4est_mesh_t *mesh = wrap->match_aux ? wrap->mesh_aux : wrap->mesh;

    // See forestclaw2d.c for more examples on how to access p4est elements.

    // .....

    // Solve the problem here.  Store solution in qsoln

#endif


    /* copy solution back to patches  */
    copy_mode = 1;
    multigrid_copy_data(glob,qsoln,copy_mode);

    FCLAW_FREE(qsoln);
}