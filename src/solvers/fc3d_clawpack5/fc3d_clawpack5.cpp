/*
Copyright (c) 2012-2023 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include <fclaw_global.h>

#include <fclaw_forestclaw.h>
#include <fclaw_clawpatch3.hpp>
#include <fclaw_clawpatch.hpp>

#include "fc3d_clawpack5.h"
#include "fc3d_clawpack5_options.h"

static fc3d_clawpack5_vtable_t classic_vt;
// static fclaw_clawpatch_vtable_t clawpatch_vt;

fc3d_clawpack5_vtable_t* fc3d_clawpack5_vt()
{
    return &classic_vt;
}

static
int fc3d_clawpack5_get_package_id (void)
{
    return s_clawpack5_package_id;
}

#if 0
void fc3d_clawpack5_set_vtable(const fc3d_clawpack5_vtable_t user_vt)
{
    classic_vt = user_vt;
}
#endif

/* This is called from the user application. */
void fc3d_clawpack5_set_vtable_defaults()
{
    fclaw2d_vtable_t*               fclaw_vt = fclaw2d_vt(glob);
    fclaw2d_patch_vtable_t*         patch_vt = fclaw2d_patch_vt(glob);
    fclaw_clawpatch3_vtable_t*  clawpatch_vt = fclaw_clawpatch3_vt();
    fc3d_clawpack5_vtable_t*        claw5_vt = fc3d_clawpack5_vt();

    fclaw2d_init_vtable();
    fclaw_clawpatch3_init_vtable_defaults();

    /* Required functions  - error if NULL*/
    // claw5_vt->bc2 = CLAWPACK5_BC2_DEFAULT;
    claw5_vt->qinit = NULL;
    claw5_vt->rpn2 = NULL;
    claw5_vt->rpt2 = NULL;

    /* Optional functions - call only if non-NULL */
    claw5_vt->setprob = NULL;
    claw5_vt->setaux = NULL;
    claw5_vt->b4step2 = NULL;
    claw5_vt->src2 = NULL;

    if (fclaw_vt->problem_setup == NULL)
    {
        /* This call shouldn't override a version-independent setting
           for this function */
        fclaw_vt->problem_setup        = &fc3d_clawpack5_setprob;
    }

    /* Default patch functions */
    patch_vt->initialize           = &fc3d_clawpack5_qinit;
    patch_vt->setup                = &fc3d_clawpack5_setaux;
    patch_vt->physical_bc          = &fc3d_clawpack5_bc2;
    patch_vt->single_step_update   = &fc3d_clawpack5_update;

    /* Forestclaw functions */
    clawpatch_vt->fort_average2coarse    = &FC3D_CLAWPACK5_FORT_AVERAGE2COARSE;
    clawpatch_vt->fort_interpolate2fine  = &FC3D_CLAWPACK5_FORT_INTERPOLATE2FINE;

    clawpatch_vt->fort_tag4refinement    = &FC3D_CLAWPACK5_FORT_TAG4REFINEMENT;
    clawpatch_vt->fort_tag4coarsening    = &FC3D_CLAWPACK5_FORT_TAG4COARSENING;

    /* output functions */
    patch_vt->write_header                = &fc3d_clawpack5_output_header_ascii;
    // fclaw_vt->fort_write_header        = &FC3D_CLAWPACK5_FORT_WRITE_HEADER;
    clawpatch_vt->fort_write_file         = &FC3D_CLAWPACK5_FORT_WRITE_FILE;
#if 0
    /* diagnostic functions */
    clawpatch_vt->fort_compute_patch_error    = NULL;  /* User defined */
    clawpatch_vt->fort_compute_error_norm     = &FC3D_CLAWPACK5_FORT_COMPUTE_ERROR_NORM;
    clawpatch_vt->fort_compute_patch_area     = &FC3D_CLAWPACK5_FORT_COMPUTE_PATCH_AREA;
    clawpatch_vt->fort_conservation_check     = &FC3D_CLAWPACK5_FORT_CONSERVATION_CHECK;

    /* Patch functions */
    clawpatch_vt->fort_copy_face          = &FC3D_CLAWPACK5_FORT_COPY_FACE;
    clawpatch_vt->fort_average_face       = &FC3D_CLAWPACK5_FORT_AVERAGE_FACE;
    clawpatch_vt->fort_interpolate_face   = &FC3D_CLAWPACK5_FORT_INTERPOLATE_FACE;

    clawpatch_vt->fort_copy_corner        = &FC3D_CLAWPACK5_FORT_COPY_CORNER;
    clawpatch_vt->fort_average_corner     = &FC3D_CLAWPACK5_FORT_AVERAGE_CORNER;
    clawpatch_vt->fort_interpolate_corner = &FC3D_CLAWPACK5_FORT_INTERPOLATE_CORNER;

    clawpatch_vt->ghostpack_extra         = NULL;
    clawpatch_vt->fort_ghostpack_qarea    = &FC3D_CLAWPACK5_FORT_GHOSTPACK_QAREA;

    clawpatch_vt->fort_timeinterp         = &FC3D_CLAWPACK5_FORT_TIMEINTERP;
#endif
}

fc3d_clawpack5_options_t* fc3d_clawpack5_get_options(fclaw_global_t *glob)
{
    return (fc3d_clawpack5_options_t*) fclaw_global_get_options(glob,
                                                                  "fc3d_clawpack5");
}

/* -----------------------------------------------------------
   Public interface to routines in this file
   ----------------------------------------------------------- */
#if 0
void fc3d_clawpack5_register (fclaw_app_t* app, const char *configfile, fclaw_global_t* glob)
{
    fc3d_clawpack5_options_t* clawopt;
    int id;

    /* Register the options */
    clawopt = fc3d_clawpack5_options_register(app,configfile);

    /* And the package */

    /* Don't register a package more than once */
    FCLAW_ASSERT(s_clawpack5_package_id == -1);

    id = fclaw_package_container_add_pkg(glob,clawopt);

    s_clawpack5_package_id = id;
}
#endif

void fc3d_clawpack5_set_options (fclaw_global_t* glob, fc3d_clawpack5_options_t* clawopt)
{
    fclaw_global_options_store(glob, "fc3d_clawpack5", clawopt);
}

void fc3d_clawpack5_aux_data(fclaw_global_t *glob,
                              fclaw_patch_t *this_patch,
                              double **aux, int* maux)
{
    fclaw2d_clawpatch_aux_data(glob, this_patch, aux, maux);
}

void fc3d_clawpack5_setprob(fclaw_global_t *glob)
{
    if (classic_vt.setprob != NULL)
    {
        classic_vt.setprob();
    }
}

/* This should only be called when a new fclaw_clawpatch_t is created. */
void fc3d_clawpack5_setaux(fclaw_global_t *glob,
                            fclaw_patch_t *this_patch,
                            int this_block_idx,
                            int this_patch_idx)
{
    if (classic_vt.setaux == NULL)
    {
        /* For consistency, class_vt.patch_setup is set to
           fc3d_clawpack5_setaux by default.  We
           check here that the user has actually set a
           SETAUX routine. If they haven't, we assume that
           they didn't really want a setaux routine. */
        return;
    }

    if (fclaw2d_patch_is_ghost(this_patch))
    {
        /* This is going to be removed at some point */
        return;
    }


    int mx,my,mz,mbc,maux;
    double xlower,ylower,zlower,dx,dy,dz;
    double *aux;


    fclaw_clawpatch3_grid_data(glob,this_patch, &mx,&my,&mz,&mbc,
                                &xlower,&ylower,&zlower,&dx,&dy,&dz);
    fclaw_clawpatch3_aux_data(glob,this_patch,&aux,&maux);

    CLAWPACK5_SET_BLOCK(&this_block_idx);
    classic_vt.setaux(&mbc,&mx,&my,&mz,&xlower,&ylower,&zlower,&dx,&dy,
                      &dz,&maux,aux);
    CLAWPACK5_UNSET_BLOCK();
}

/* This should only be called when a new fclaw_clawpatch_t is created. */
void fc3d_clawpack5_set_capacity(fclaw_global_t *glob,
                                  fclaw_patch_t *this_patch,
                                  int this_block_idx,
                                  int this_patch_idx)
{
    int mx,my,mbc,maux,mcapa;
    double dx,dy,xlower,ylower;
    double *aux, *area;
    fc3d_clawpack5_options_t *clawopt;

    clawopt = fc3d_clawpack5_get_options(glob);
    mcapa = clawopt->mcapa;

    fclaw_clawpatch_grid_data_2d(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    area = fclaw_clawpatch_get_area_2d(glob,this_patch);

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);
    FCLAW_ASSERT(maux >= mcapa && mcapa > 0);

    CLAWPACK5_SET_CAPACITY(&mx,&my,&mbc,&dx,&dy,area,&mcapa,
                            &maux,aux);
}


void fc3d_clawpack5_qinit(fclaw_global_t *glob,
                           fclaw_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx)
{
    FCLAW_ASSERT(classic_vt.qinit != NULL); /* Must initialized */
    int mx,my,mz,mbc,meqn,maux;
    double dx,dy,dz,xlower,ylower,zlower;
    double *q, *aux;

    fclaw_clawpatch3_grid_data(glob,this_patch,&mx,&my,&mz,&mbc,
                                &xlower,&ylower,&zlower,&dx,&dy,&dz);

    fclaw_clawpatch3_soln_data(glob,this_patch,&q,&meqn);
    fclaw_clawpatch3_aux_data(glob,this_patch,&aux,&maux);

    /* Call to classic Clawpack 'qinit' routine.  This must be user defined */
    CLAWPACK5_SET_BLOCK(&this_block_idx);
    classic_vt.qinit(&meqn,&mbc,&mx,&my,&mz,&xlower,&ylower,&zlower,&dx,&dy,&dz,q,
                     &maux,aux);
    CLAWPACK5_UNSET_BLOCK();
}

void fc3d_clawpack5_b4step2(fclaw_global_t *glob,
                             fclaw_patch_t *this_patch,
                             int this_block_idx,
                             int this_patch_idx,
                             double t, double dt)

{
    FCLAW_ASSERT(classic_vt.b4step2 != NULL);

    int mx,my,mbc,meqn, maux;
    double xlower,ylower,dx,dy;
    double *aux,*q;

    fclaw_clawpatch_grid_data_2d(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);
    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    CLAWPACK5_SET_BLOCK(&this_block_idx);
    classic_vt.b4step2(&mbc,&mx,&my,&meqn,q,&xlower,&ylower,
                       &dx,&dy,&t,&dt,&maux,aux);
    CLAWPACK5_UNSET_BLOCK();
}

void fc3d_clawpack5_src2(fclaw_global_t *glob,
                          fclaw_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt)
{
    FCLAW_ASSERT(classic_vt.src2 != NULL);

    int mx,my,mbc,meqn, maux;
    double xlower,ylower,dx,dy;
    double *aux,*q;

    fclaw_clawpatch_grid_data_2d(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);
    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    CLAWPACK5_SET_BLOCK(&this_block_idx);
    classic_vt.src2(&meqn,&mbc,&mx,&my,&xlower,&ylower,
                    &dx,&dy,q,&maux,aux,&t,&dt);
    CLAWPACK5_UNSET_BLOCK();
}



/* Use this to return only the right hand side of the clawpack algorithm */
double fc3d_clawpack5_step2_rhs(fclaw_global_t *glob,
                                 fclaw_patch_t *this_patch,
                                 int this_block_idx,
                                 int this_patch_idx,
                                 double t,
                                 double *rhs)
{
    /* This should evaluate the right hand side, but not actually do the update.
       This will be useful in cases where we want to use something other than
       a single step method.  For example, in a RK scheme, one might want to
       call the right hand side to evaluate stages. */
    return 0;
}


void fc3d_clawpack5_bc2(fclaw_global_t *glob,
                         fclaw_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx,
                         double t,
                         double dt,
                         fclaw_bool intersects_phys_bdry[],
                         fclaw_bool time_interp)
{
    fc3d_clawpack5_options_t *clawopt = fc3d_clawpack5_get_options(glob);

    FCLAW_ASSERT(classic_vt.bc2 != NULL);

    int mx,my,mbc,meqn, maux;
    double xlower,ylower,dx,dy;
    double *aux,*q;

    fclaw_clawpatch_grid_data_2d(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    int *block_mthbc = clawopt->mthbc;

    /* Set a local copy of mthbc that can be used for a patch. */
    int mthbc[FCLAW2D_NUMFACES];
    for(int i = 0; i < FCLAW2D_NUMFACES; i++)
    {
        if (intersects_phys_bdry[i])
        {
            mthbc[i] = block_mthbc[i];
        }
        else
        {
            mthbc[i] = -1;
        }
    }

    /*
      We may be imposing boundary conditions on time-interpolated data;
      and is being done just to help with fine grid interpolation.
      In this case, this boundary condition won't be used to update
      anything
    */
    fclaw2d_clawpatch_timesync_data(glob,this_patch,time_interp,&q,&meqn);

    CLAWPACK5_SET_BLOCK(&this_block_idx);
    classic_vt.bc2(&meqn,&mbc,&mx,&my,&xlower,&ylower,
                   &dx,&dy,q,&maux,aux,&t,&dt,mthbc);
    CLAWPACK5_UNSET_BLOCK();
}


/* This is called from the single_step callback. and is of type 'flaw_single_step_t' */
double fc3d_clawpack5_step2(fclaw_global_t *glob,
                             fclaw_patch_t *this_patch,
                             int this_block_idx,
                             int this_patch_idx,
                             double t,
                             double dt)
{
    fc3d_clawpack5_options_t* clawpack_options;
    int level;
    double *qold, *aux;
    int mx, my, meqn, maux, mbc;
    double xlower, ylower, dx,dy;

    FCLAW_ASSERT(classic_vt.rpn2 != NULL);
    FCLAW_ASSERT(classic_vt.rpt2 != NULL);

    clawpack_options = fc3d_clawpack5_get_options(glob);

    // cp = fclaw2d_clawpatch_get_cp(this_patch);

    level = this_patch->level;

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    fclaw2d_clawpatch_save_current_step(glob, this_patch);

    fclaw_clawpatch_grid_data_2d(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&qold,&meqn);

    int mwaves = clawpack_options->mwaves;

    int maxm = fmax(mx,my);

    double cflgrid = 0.0;

    int mwork = (maxm+2*mbc)*(12*meqn + (meqn+1)*mwaves + 3*maux + 2);
    double* work = new double[mwork];

    int size = meqn*(mx+2*mbc)*(my+2*mbc);
    double* fp = new double[size];
    double* fm = new double[size];
    double* gp = new double[size];
    double* gm = new double[size];


    int ierror = 0;
    //fc3d_clawpack5_flux2_t flux2 = clawpack_options->use_fwaves ?
    //                                CLAWPACK5_FLUX2FW : CLAWPACK5_FLUX2;
    fc3d_clawpack5_flux2_t flux2 = CLAWPACK5_FLUX2;
    int* block_corner_count = fclaw2d_patch_block_corner_count(glob,this_patch);
    CLAWPACK5_STEP2_WRAP(&maxm, &meqn, &maux, &mbc, clawpack_options->method,
                          clawpack_options->mthlim, &clawpack_options->mcapa,
                          &mwaves,&mx, &my, qold, aux, &dx, &dy, &dt, &cflgrid,
                          work, &mwork, &xlower, &ylower, &level,&t, fp, fm, gp, gm,
                          classic_vt.rpn2, classic_vt.rpt2,flux2,
                          block_corner_count, &ierror);

    FCLAW_ASSERT(ierror == 0);

    //CLAWPACK5_STEP2(&maxm, &meqn, &maux, &mbc, &mx, &my, qold, aux, &dx, &dy, &dt, &cflgrid,
    //                fm, fp, gm, gp, classic_vt.rpn2, classic_vt.rpt2);
    /* Accumulate fluxes needed for conservative fix-up */
    if (classic_vt.fluxfun != NULL)
    {
        /* Accumulate fluxes */
    }


    delete [] fp;
    delete [] fm;
    delete [] gp;
    delete [] gm;

    delete [] work;

    return cflgrid;
}

double fc3d_clawpack5_update(fclaw_global_t *glob,
                              fclaw_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx,
                              double t,
                              double dt)
{    
    const fc3d_clawpack5_options_t* clawpack_options;
    clawpack_options = fc3d_clawpack5_get_options(glob);
    if (classic_vt.b4step2 != NULL)
    {
        fc3d_clawpack5_b4step2(glob,
                               this_patch,
                               this_block_idx,
                               this_patch_idx,t,dt);
    }
    double maxcfl = fc3d_clawpack5_step2(glob,
                                         this_patch,
                                         this_block_idx,
                                         this_patch_idx,t,dt);
    if (clawpack_options->src_term > 0)
    {
        fc3d_clawpack5_src2(glob,
                             this_patch,
                             this_block_idx,
                             this_patch_idx,t,dt);
    }
    return maxcfl;
}

void fc3d_clawpack5_output_header_ascii(fclaw_global_t* glob,
                                        int iframe)
{
    const fclaw_clawpatch_options_t *clawpatch_opt;

    int meqn,maux,ngrids;
    double time;

    clawpatch_opt = fclaw_clawpatch_get_options(glob);

    time = glob->curr_time;
    ngrids = fclaw2d_domain_get_num_patches(glob->domain);

    meqn = clawpatch_opt->meqn;
    maux = clawpatch_opt->maux;

    FC3D_CLAWPACK5_FORT_WRITE_HEADER(&iframe,&time,&meqn,&maux,&ngrids);
}
