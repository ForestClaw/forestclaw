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

#include "amr_forestclaw.H"
#include "fc2d_clawpack46.H"
#include "pwconst_user.H"
#include <fclaw2d_vtable.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

static fclaw2d_vtable_t vt;
static fc2d_clawpack46_vtable_t classic_user;

void pwconst_link_solvers(fclaw2d_domain_t *domain)
{
    vt.problem_setup              = NULL;
    vt.patch_setup                = NULL;
    vt.patch_initialize           = &fc2d_clawpack46_qinit;
    vt.patch_physical_bc          = &fc2d_clawpack46_bc2;
    vt.patch_single_step_update   = &fc2d_clawpack46_update;

    vt.write_header               = &pwconst_parallel_write_header;
    vt.patch_write_output         = &pwconst_parallel_write_output;

    vt.patch_tag4refinement       = &pwconst_patch_tag4refinement;
    vt.patch_tag4coarsening       = &pwconst_patch_tag4coarsening;

    fclaw2d_set_vtable(domain,&vt);

    classic_user.qinit = &QINIT;
    classic_user.rpn2 = &RPN2;
    classic_user.rpt2 = &RPT2;

    fc2d_clawpack46_set_vtable(&classic_user);

}

/* -----------------------------------------------------------------
   Default routine for tagging patches for refinement and coarsening
   ----------------------------------------------------------------- */
fclaw_bool pwconst_patch_tag4refinement(fclaw2d_domain_t *domain,
                                          fclaw2d_patch_t *this_patch,
                                          int this_block_idx, int this_patch_idx,
                                          int initflag)
{
    /* ----------------------------------------------------------- */
    // Global parameters
    const amr_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;

    /* ----------------------------------------------------------- */
    // Patch specific parameters
    ClawPatch *cp = get_clawpatch(this_patch);
    double xlower = cp->xlower();
    double ylower = cp->ylower();
    double dx = cp->dx();
    double dy = cp->dy();

    /* ------------------------------------------------------------ */
    // Pointers needed to pass to Fortran
    double* q = cp->q();

    int tag_patch = 0;
    pwconst_tag4refinement_(mx,my,mbc,meqn,xlower,ylower,dx,dy,q,initflag,tag_patch);
    return tag_patch == 1;
}

fclaw_bool pwconst_patch_tag4coarsening(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int blockno,
                                      int patchno)
{
    /* ----------------------------------------------------------- */
    // Global parameters
    const amr_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;

    /* ----------------------------------------------------------- */
    // Patch specific parameters
    ClawPatch *cp = get_clawpatch(this_patch);
    double xlower = cp->xlower();
    double ylower = cp->ylower();
    double dx = cp->dx();
    double dy = cp->dy();

    /* ------------------------------------------------------------ */
    // Pointers needed to pass to Fortran
    double* qcoarse = cp->q();

    int tag_patch = 1;  // == 0 or 1
    pwconst_tag4coarsening_(mx,my,mbc,meqn,xlower,ylower,dx,dy,qcoarse,tag_patch);
    return tag_patch == 0;
}


void pwconst_parallel_write_header(fclaw2d_domain_t* domain, int iframe, int ngrids)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    double time = get_domain_time(domain);

    fclaw_global_essentialf("Matlab output Frame %d  at time %16.8e\n\n",iframe,time);

    int mfields = gparms->meqn;
    int maux = 0;
    pwconst_write_tfile_(iframe,time,mfields,ngrids,maux);

    new_qfile_(iframe);
}


void pwconst_parallel_write_output(fclaw2d_domain_t *domain, fclaw2d_patch_t *this_patch,
                                     int this_block_idx, int this_patch_idx,
                                     int iframe,int num,int level)
{
    /* ----------------------------------------------------------- */
    // Global parameters
    const amr_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;

    /* ----------------------------------------------------------- */
    // Patch specific parameters
    ClawPatch *cp = get_clawpatch(this_patch);
    double xlower = cp->xlower();
    double ylower = cp->ylower();
    double dx = cp->dx();
    double dy = cp->dy();

    /* ------------------------------------------------------------ */
    // Pointers needed to pass to Fortran
    double* q = cp->q();

    /* ------------------------------------------------------------- */

    int mpirank = domain->mpirank;
    pwconst_write_qfile_(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,
                           iframe,num,level,this_block_idx,mpirank);
}


#ifdef __cplusplus
#if 0
{
#endif
}
#endif
