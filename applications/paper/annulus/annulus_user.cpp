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

#include "annulus_user.h"

#include <fclaw2d_include_all.h>

#include "fclaw2d_clawpatch.h"

#include <fc2d_clawpack46.h>

#include "fclaw2d_clawpatch_options.h"
#include <fc2d_clawpack46_options.h>


#include "clawpack46_advection_user_fort.h"


static
void annulus_problem_setup(fclaw2d_global_t *glob)
{
    const user_options_t *user = annulus_get_options(glob);

    if (glob->mpirank == 0)
    {

        FILE *f = fopen("setprob.data","w");
        fprintf(f,  "%-24d   %s",user->example,       "\% example\n");  
        fprintf(f,"%-24.4f   %s",user->revs_per_s,    "\% revs_per_s\n"); 
        fprintf(f,"%-24.4f   %s",user->cart_speed,    "\% cart_speed\n"); 
        fprintf(f,"%-24.4f   %s",user->amplitude,     "\% amplitude\n");    
        fprintf(f,"%-24.4f   %s",user->freq,          "\% freq\n");         
        fprintf(f,"%-24.4f   %s",user->beta,          "\% beta\n");    
        fprintf(f,"%-12.4f%-12.4f    %s",user->theta[0],user->theta[1],"\% beta\n");    
        fprintf(f,  "%-24d   %s",user->refine_pattern,"\% refine_pattern\n");    
        fprintf(f,"%-24.4f   %s",user->init_radius,   "\% init_radius\n");    
        fclose(f);
    }

    /* We want to make sure node 0 gets here before proceeding */
#ifdef FCLAW_ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    SETPROB_ANNULUS();  /* Reads file created above */
 
}


static
void annulus_patch_setup(fclaw2d_global_t *glob,
                         fclaw2d_patch_t *this_patch,
                         int blockno,
                         int patchno)
{
    int mx,my,mbc,maux;
    double xlower,ylower,dx,dy;
    double *aux,*xd,*yd,*zd,*area;
    double *xp,*yp,*zp;
    double *curvature, *surfnormals;
    double *edgelengths;
    double *xnormals, *ynormals, *xtangents, *ytangents;

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_metric_data(glob,this_patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    fclaw2d_clawpatch_metric_scalar(glob, this_patch,&area,&edgelengths,
                                    &curvature);

    fclaw2d_clawpatch_metric_vector(glob,this_patch,
                                    &xnormals, &ynormals,
                                    &xtangents, &ytangents,
                                    &surfnormals);

    ANNULUS46_SETAUX(&mbc,&mx,&my,&xlower,&ylower,
                     &dx,&dy,&maux,aux,&blockno,
                     area,xd,yd,zd,
                     edgelengths,xnormals,ynormals,
                     xtangents, ytangents, surfnormals);
}


void annulus_link_solvers(fclaw2d_global_t *glob)
{
    /* ForestClaw virtual functions */
    fclaw2d_vtable_t  *vt = fclaw2d_vt();
    vt->problem_setup = &annulus_problem_setup;

    /* Patch virtual functions */
    fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
    patch_vt->setup = &annulus_patch_setup;

    /* Clawpatch virtual table */
    fclaw2d_clawpatch_vtable_t *clawpatch_vt  = fclaw2d_clawpatch_vt();
    clawpatch_vt->fort_tag4refinement = &CLAWPACK46_TAG4REFINEMENT;
    clawpatch_vt->fort_tag4coarsening = &CLAWPACK46_TAG4COARSENING;

    /* Clawpack options */
    fc2d_clawpack46_options_t  *clawopt       = fc2d_clawpack46_get_options(glob);
    clawopt->use_fwaves = 1;

    /* Clawpack virtual functions */
    fc2d_clawpack46_vtable_t   *clawpack46_vt = fc2d_clawpack46_vt();
    clawpack46_vt->fort_qinit = &CLAWPACK46_QINIT;


    /* Get different refinement patterns */
    clawpack46_vt->fort_rpn2 = &RPN2CONS_FW_MANIFOLD;   
    clawpack46_vt->fort_rpt2 = &RPT2CONS_MANIFOLD;    


    /* Time dependent velocity field */
    const user_options_t *user = annulus_get_options(glob);
    if (user->example == 2 || user->example == 3) 
        clawpack46_vt->fort_b4step2 = &CLAWPACK46_B4STEP2;        


    /* Forestclaw equivalent of QAD */
    if (user->flux_correction) 
        clawpack46_vt->fort_rpn2_cons = &RPN2QAD_FLUX;        
    else
        clawpack46_vt->fort_rpn2_cons = &RPN2QAD_ZERO;

}


