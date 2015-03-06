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
#include "ClawPatch.H"
#include "fclaw2d_map_query.h"
#include "fclaw2d_clawpatch.h"
#include "fclaw2d_vtable.h"

#include "metric_user.H"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

static fclaw2d_vtable_t vt;

void metric_link_patch(fclaw2d_domain_t *domain)
{
    fclaw2d_init_vtable(&vt);

    vt.problem_setup = &metric_problem_setup;

    vt.patch_initialize = &metric_patch_initialize;
    vt.patch_physical_bc = &fclaw2d_physbc_default;

    vt.run_diagnostics = &metric_diagnostics;

    fclaw2d_set_vtable(domain,&vt);

}

void metric_problem_setup(fclaw2d_domain_t* domain)
{
    /* Any general problem set up here */
    SETPROB();  /* Set value of pi */
}

#if 0
void metric_link_patch(fclaw2d_domain_t *domain)
{
    fclaw2d_solver_functions_t* sf = get_solver_functions(domain);
    sf->f_patch_initialize = &metric_patch_initialize;

    fclaw2d_regrid_functions_t *rf = get_regrid_functions(domain);
    rf->f_patch_tag4refinement = &metric_patch_tag4refinement;
    rf->f_patch_tag4coarsening = &metric_patch_tag4coarsening;

    fclaw2d_output_functions_t *of = get_output_functions(domain);
    of->f_patch_write_header = &metric_parallel_write_header;
    of->f_patch_write_output = &metric_parallel_write_output;
}
#endif

void metric_patch_initialize(fclaw2d_domain_t *domain,
                             fclaw2d_patch_t *this_patch,
                             int this_block_idx,
                             int this_patch_idx)
{
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *q, *area, *curvature, *error_ptr;
    ClawPatch *cp;
    int blockno;

#if 0
    cp = fclaw2d_clawpatch_cp(domain,this_patch);
#endif
    cp = get_clawpatch(this_patch);

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);

    area = cp->area();
    curvature = cp->curvature();

    /* Create an array with same dimensions as q, and one field */
    FArrayBox error;
    error.define(cp->dataBox(),1);
    error_ptr = error.dataPtr();

    fclaw2d_map_context_t* cont = get_map_context(domain);
    blockno = this_block_idx;
    compute_error(meqn,mbc,mx,my,&cont,blockno,xlower,ylower,dx,dy,
                  curvature,error_ptr);
    initialize(mx,my,meqn,mbc,xlower,ylower,dx,dy,q,
               error_ptr,curvature,area);
}


#if 0
/* -----------------------------------------------------------------
   Default routine for tagging patches for refinement and coarsening
   ----------------------------------------------------------------- */
fclaw_bool metric_patch_tag4refinement(fclaw2d_domain_t *domain,
                                       fclaw2d_patch_t *this_patch,
                                       int blockno, int this_patch_idx,
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

    int tag_patch = 0;  // == 0 or 1
    metric_tag4refinement_(mx,my,mbc,meqn,xlower,ylower,dx,dy,q,initflag,
                           blockno, tag_patch);
    return tag_patch == 1;
}

fclaw_bool metric_patch_tag4coarsening(fclaw2d_domain_t *domain,
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
    double* q = cp->q();

    int tag_patch = 1;  // == 0 or 1
    metric_tag4coarsening_(mx,my,mbc,meqn,xlower,ylower,dx,dy,q,tag_patch);
    return tag_patch == 0;
}

void metric_parallel_write_header(fclaw2d_domain_t* domain, int iframe,int ngrids)
{
    double time = get_domain_time(domain);

    fclaw_global_essentialf("Matlab output Frame %d  at time %16.8e\n\n",iframe,time);

    /* Write out header file containing global information for 'iframe' */
    int mfields = 3;  // Write out an extra fields
    int maux = 0;
    write_tfile_(iframe,time,mfields,ngrids,maux);

    /* Initialize fort.qXXXX */
    new_qfile_(iframe);
}


void metric_parallel_write_output(fclaw2d_domain_t *domain, fclaw2d_patch_t *this_patch,
                                  int this_block_idx, int this_patch_idx,
                                  int iframe,int patch_num,int level)
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
    double *q = cp->q();
    int blockno = this_block_idx;
    int mpirank = domain->mpirank;

    /* ------------------------------------------------------------- */
    metric_output(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,
                  iframe,patch_num,level,blockno,mpirank);
}
#endif

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
