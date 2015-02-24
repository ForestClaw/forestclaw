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

#include "amr_includes.H"
#include "fc2d_clawpack46.H"
#include "torus_user.H"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

static fc2d_clawpack46_vtable_t classic_user;

void torus_link_solvers(fclaw2d_domain_t *domain)
{
    const amr_options_t *gparms = get_domain_parms(domain);

    fclaw2d_solver_functions_t* sf = get_solver_functions(domain);
    sf->use_single_step_update = fclaw_true;
    sf->use_mol_update = fclaw_false;

    if (!gparms->manifold)
    {
        sf->f_patch_setup              = &fc2d_clawpack46_setaux;
    }
    else
    {
        sf->f_patch_setup              = &torus_patch_manifold_setup;
    }

    sf->f_patch_initialize         = &fc2d_clawpack46_qinit;
    sf->f_patch_physical_bc        = &fc2d_clawpack46_bc2;  /* Needed for lat-long grid */
    sf->f_patch_single_step_update = &fc2d_clawpack46_update;


    fclaw2d_regrid_functions_t *rf = get_regrid_functions(domain);
    rf->f_patch_tag4refinement     = &torus_patch_tag4refinement;
    rf->f_patch_tag4coarsening     = &torus_patch_tag4coarsening;

    fclaw2d_output_functions_t *of = get_output_functions(domain);
    of->f_patch_write_header       = &torus_parallel_write_header;
    of->f_patch_write_output       = &torus_parallel_write_output;

    classic_user.setprob = &SETPROB;
    classic_user.qinit = &QINIT;
    classic_user.setaux = &SETAUX;
    classic_user.rpn2 = &RPN2;
    classic_user.rpt2 = &RPT2;

    fc2d_clawpack46_set_vtable(&classic_user);

}

void torus_patch_manifold_setup(fclaw2d_domain_t *domain,
                                fclaw2d_patch_t *this_patch,
                                int this_block_idx,
                                int this_patch_idx)
{

    int mx,my,mbc,maux;
    double xlower,ylower,dx,dy;
    double *xd,*yd,*zd,*area;
    double *xp,*yp,*zp;
    double *aux;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_metric_data(domain,this_patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    fc2d_clawpack46_define_auxarray2(domain,this_patch);
    fc2d_clawpack46_aux_data(domain,this_patch,&aux,&maux);

    setaux_manifold_(mbc,mx,my,this_block_idx,xlower,ylower,dx,dy,
                     maux,aux, area);
}


/* -----------------------------------------------------------------
   Default routine for tagging patches for refinement and coarsening
   ----------------------------------------------------------------- */
fclaw_bool torus_patch_tag4refinement(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int this_block_idx, int this_patch_idx,
                                      int initflag)
{
    int mx,my,mbc, meqn;
    double xlower,ylower,dx,dy;
    double *q;
    int tag_patch;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q, &meqn);

    tag_patch = 0;
    torus_tag4refinement_(mx,my,mbc,meqn,xlower,ylower,dx,dy,q,initflag,
                           this_block_idx,tag_patch);
    return tag_patch == 1;
}

fclaw_bool torus_patch_tag4coarsening(fclaw2d_domain_t *domain,
                                             fclaw2d_patch_t *this_patch,
                                             int blockno_idx,
                                             int patchno)
{
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *q;
    int tag_patch;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);

    torus_tag4coarsening_(mx,my,mbc,meqn,xlower,ylower,dx,dy,q,tag_patch);
    return tag_patch == 0;
}

void torus_parallel_write_header(fclaw2d_domain_t* domain, int iframe, int ngrids)
{
    int meqn, maux;
    double time;

    time = get_domain_time(domain);

    fclaw_global_essentialf("Matlab output Frame %d  at time %16.8e\n\n",iframe,time);

    fclaw2d_clawpatch_meqn(domain, &meqn);
    fc2d_clawpack46_maux(domain, &maux);

    torus_write_tfile_(iframe,time,meqn,ngrids,maux);

    new_qfile_(iframe);
}


void torus_parallel_write_output(fclaw2d_domain_t *domain, fclaw2d_patch_t *this_patch,
                                  int this_block_idx, int this_patch_idx,
                                  int iframe,int num,int level)
{
    int mx,my,mbc,meqn,maxmx,maxmy;
    double xlower,ylower,dx,dy;
    double *q;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);

    maxmx = mx;
    maxmy = my;
    /* This opens a file for append and writes in the 'clawout' style. */
    torus_write_qfile_(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,
                        iframe,num,level,this_block_idx,domain->mpirank);
}



#ifdef __cplusplus
#if 0
{
#endif
}
#endif
