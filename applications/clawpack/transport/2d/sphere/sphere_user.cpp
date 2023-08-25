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

#include "sphere_user.h"

static
void sphere_problem_setup(fclaw_global_t* glob)
{
    const user_options_t* user = sphere_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,  "%-24d   %s",user->example,"\% example\n");
        fprintf(f,  "%-24d   %s",user->mapping,"\% mapping\n");
        fprintf(f,  "%-24d   %s",user->initial_condition,"\% initial_condition\n");
        fprintf(f,  "%-24.16f   %s",user->omega[0],"\% omega[0]\n");
        fprintf(f,  "%-24.16f   %s",user->omega[1],"\% omega[1]\n");
        fprintf(f,  "%-24.16f   %s",user->omega[2],"\% omega[2]\n");
        fprintf(f,  "%-24d   %s",user->refine_pattern,"\% refinement_pattern\n");
        fclose(f);
    }
    fclaw_domain_barrier (glob->domain);
    SETPROB();
}

void sphere_patch_setup_manifold(fclaw_global_t *glob,
                                    fclaw_patch_t *patch,
                                    int blockno,
                                    int patchno)
{
    const user_options_t *user = sphere_get_options(glob);
    transport_patch_setup_manifold(glob,patch,blockno,patchno,
                                   user->claw_version);
}


static
void sphere_b4step2(fclaw_global_t *glob,
                    fclaw_patch_t *patch,
                    int blockno,
                    int patchno,
                    double t, double dt)

{
    const user_options_t* user = sphere_get_options(glob);
    transport_b4step2_manifold(glob,patch,blockno,patchno,t, dt,
                               user->claw_version);
}


static
void cb_sphere_output_ascii (fclaw_domain_t * domain,
                            fclaw_patch_t * patch,
                            int blockno, int patchno,
                            void *user)
{
    fclaw_global_iterate_t* s = (fclaw_global_iterate_t*) user;
    fclaw_global_t      *glob = (fclaw_global_t*) s->glob;

    //fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);

    int iframe = *((int *) s->user);

    /* Get info not readily available to user */
    int global_patch_num, local_patch_num, level;
    fclaw_patch_get_info(glob->domain,patch,
                           blockno,patchno,
                           &global_patch_num,&local_patch_num,&level);
    
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw_clawpatch_grid_data_2d(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int meqn;
    double *q;
    fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

    double* error = fclaw_clawpatch_get_error(glob,patch);
    double* soln = fclaw_clawpatch_get_exactsoln(glob,patch);

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    char fname[BUFSIZ];
    snprintf (fname, BUFSIZ, "%s.q%04d", fclaw_opt->prefix, iframe);


    /* Here, we pass in q and the error, so need special headers and files */
    double time = glob->curr_time;

    const user_options_t* user_opt = sphere_get_options(glob);
    if (user_opt->claw_version == 4)
    {
        SPHERE46_FORT_WRITE_FILE(fname, &mx,&my,&meqn,&mbc,
                                 &xlower,&ylower,
                                 &dx,&dy, q,error,soln, &time,
                                 &global_patch_num,&level,
                                 &blockno, &glob->mpirank);
    }
    else
    {
        SPHERE5_FORT_WRITE_FILE(fname, &mx,&my,&meqn,&mbc,
                               &xlower,&ylower,
                               &dx,&dy, q,error,soln, &time,
                               &global_patch_num,&level,
                               &blockno, &glob->mpirank);
    }
}



void sphere_link_solvers(fclaw_global_t *glob)
{
    /* ForestClaw core functions */
    fclaw_vtable_t *vt = fclaw_vt(glob);
    vt->problem_setup = &sphere_problem_setup;  /* Version-independent */

    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    patch_vt->setup   = &sphere_patch_setup_manifold;


    const user_options_t* user_opt = sphere_get_options(glob);
    if (user_opt->mapping == 1)
        fclaw_clawpatch_use_pillowsphere(glob);

    if (user_opt->claw_version == 4)
    {
        fc2d_clawpack46_vtable_t  *clawpack46_vt = fc2d_clawpack46_vt(glob);
        clawpack46_vt->b4step2        = sphere_b4step2;
        clawpack46_vt->fort_qinit     = CLAWPACK46_QINIT;
        clawpack46_vt->fort_rpn2fw    = &CLAWPACK46_RPN2CONS_FW_MANIFOLD; 
        clawpack46_vt->fort_rpt2fw    = &CLAWPACK46_RPT2CONS_MANIFOLD;      
        clawpack46_vt->fort_rpn2_cons = &RPN2CONS_UPDATE_MANIFOLD;

        /* Clawpatch functions */    
        fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);

        /* Include error in output files */
        const fclaw_options_t* fclaw_opt = fclaw_get_options(glob);
        if (fclaw_opt->compute_error)
        {
            clawpatch_vt->d2->fort_compute_patch_error = &SPHERE46_COMPUTE_ERROR;
            clawpatch_vt->fort_header_ascii   = &SPHERE46_FORT_HEADER_ASCII;
            clawpatch_vt->cb_output_ascii     = &cb_sphere_output_ascii;                
        }
    }
    else
    {
        fc2d_clawpack5_vtable_t  *clawpack5_vt = fc2d_clawpack5_vt(glob);
        clawpack5_vt->b4step2        = sphere_b4step2;
        clawpack5_vt->fort_qinit     = CLAWPACK5_QINIT;
        clawpack5_vt->fort_rpn2      = &CLAWPACK5_RPN2CONS_MANIFOLD; 
        clawpack5_vt->fort_rpt2      = &CLAWPACK5_RPT2CONS_MANIFOLD;      
        clawpack5_vt->fort_rpn2_cons = &RPN2CONS_UPDATE_MANIFOLD;

        /* Clawpatch functions */    
        fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);

        /* Include error in output files */
        const fclaw_options_t* fclaw_opt = fclaw_get_options(glob);
        if (fclaw_opt->compute_error)
        {
            clawpatch_vt->d2->fort_compute_patch_error = &SPHERE5_COMPUTE_ERROR;
            clawpatch_vt->fort_header_ascii   = &SPHERE5_FORT_HEADER_ASCII;
            clawpatch_vt->cb_output_ascii     = &cb_sphere_output_ascii;                
        }
    }
 }










