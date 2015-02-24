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
#include "fc2d_dummy.H"
#include "swirl_user.H"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

static const fc2d_clawpack46_vtable_t classic_user =
{
    SETPROB,
    BC2,  /* bc2 - added here just as a compiler check */
    QINIT,
    SETAUX,
    B4STEP2,
    NULL,  /* src2 */
    RPN2,
    RPT2
};


void swirl_link_solvers(fclaw2d_domain_t *domain)
{
    fclaw2d_solver_functions_t* sf = get_solver_functions(domain);

    sf->use_single_step_update = fclaw_true;
    sf->use_mol_update = fclaw_false;

    sf->f_patch_setup              = &swirl_patch_setup;
    sf->f_patch_initialize         = &fc2d_clawpack46_qinit;
    sf->f_patch_physical_bc        = &swirl_patch_physical_bc;
    sf->f_patch_single_step_update = &fc2d_clawpack46_update;

    fclaw2d_output_functions_t* of = get_output_functions(domain);
    of->f_patch_write_header = &swirl_parallel_write_header;
    of->f_patch_write_output = &swirl_parallel_write_output;

    fc2d_clawpack46_set_vtable(&classic_user);

}


void swirl_patch_setup(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx)
{
    /* Dummy setup - use multiple libraries */
    fc2d_clawpack46_setaux(domain,this_patch,this_block_idx,this_patch_idx);
    fc2d_dummy_setup_patch(domain,this_patch,this_block_idx,this_patch_idx);

}



void swirl_patch_physical_bc(fclaw2d_domain *domain,
                             fclaw2d_patch_t *this_patch,
                             int this_block_idx,
                             int this_patch_idx,
                             double t,
                             double dt,
                             fclaw_bool intersects_bc[],
                             fclaw_bool time_interp)
{
    /* This calls bc2 in swirl/user_4.6;  that file isn't changed but
       is included to show that both the local version of bc2.f and the
       clawpack46 library code can be included */
    fc2d_clawpack46_bc2(domain,this_patch,this_block_idx,this_patch_idx,
                     t,dt,intersects_bc,time_interp);
}


/* -----------------------------------------------------------------
   Default routine for tagging patches for refinement and coarsening
   ----------------------------------------------------------------- */
fclaw_bool swirl_patch_tag4refinement(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int this_block_idx, int this_patch_idx,
                                      int initflag)
{
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *q;
    int tag_patch;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);

    tag_patch = 0;
    swirl_tag4refinement_(mx,my,mbc,meqn,xlower,ylower,dx,dy,q,initflag,
                          tag_patch);
    return tag_patch == 1;
}

fclaw_bool swirl_patch_tag4coarsening(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int blockno,
                                      int patchno)
{
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *qcoarse;
    int tag_patch;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&qcoarse,&meqn);

    tag_patch = 1;
    swirl_tag4coarsening_(mx,my,mbc,meqn,xlower,ylower,dx,dy,qcoarse,tag_patch);
    return tag_patch == 0;
}

void swirl_parallel_write_header(fclaw2d_domain_t* domain, int iframe, int ngrids)
{
    int maux, meqn;
    double time;

    time = get_domain_time(domain);
    fc2d_clawpack46_maux(domain,&maux);
    fclaw2d_clawpatch_meqn(domain,&meqn);

    fclaw_global_essentialf("Matlab output Frame %d  at time %16.8e\n\n",iframe,time);

    swirl_write_tfile_(iframe,time,meqn,ngrids,maux);

    /* This opens file 'fort.qXXXX' for replace and closes the file. */
    new_qfile_(iframe);
}


void swirl_parallel_write_output(fclaw2d_domain_t *domain, fclaw2d_patch_t *this_patch,
                                  int this_block_idx, int this_patch_idx,
                                  int iframe,int global_patch_num,int level)
{
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *q;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);

    swirl_write_qfile_(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,
                       iframe,global_patch_num,level,this_block_idx,
                       domain->mpirank);
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
