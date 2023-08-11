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

#include "fc2d_clawpack5.h"
#include "fc2d_clawpack5_fort.h"
#include "fc2d_clawpack5_options.h"


#include <fclaw_pointer_map.h>

#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>
#include <fclaw_clawpatch.hpp>

#include <fclaw_clawpatch_output_ascii.h>
#include <fclaw_clawpatch_output_vtk.h>
#include <fclaw2d_clawpatch_fort.h>

#include <fclaw2d_clawpatch_conservation.h>

#include <fclaw2d_patch.h>
#include <fclaw2d_global.h>
#include <fclaw2d_vtable.h>
#include <fclaw2d_options.h>
#include <fclaw2d_defs.h>


/* -------------------------- Clawpack solver functions ------------------------------ */

static
void clawpack5_setprob(fclaw2d_global_t *glob)
{
    fc2d_clawpack5_vtable_t*  claw5_vt = fc2d_clawpack5_vt(glob);
    if (claw5_vt->fort_setprob != NULL)
    {
        claw5_vt->fort_setprob();
    }
}

/* This should only be called when a new fclaw_clawpatch_t is created. */
static
void clawpack5_setaux(fclaw2d_global_t *glob,
                      fclaw_patch_t *this_patch,
                      int this_block_idx,
                      int this_patch_idx)
{
    fc2d_clawpack5_vtable_t*  claw5_vt = fc2d_clawpack5_vt(glob);
    if (claw5_vt->fort_setaux == NULL)
    {
        return;
    }

    if (fclaw2d_patch_is_ghost(this_patch))
    {
        /* This is going to be removed at some point */
        return;
    }


    int mx,my,mbc,maux;
    double xlower,ylower,dx,dy;
    double *aux;


    fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);
    fclaw_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    CLAWPACK5_SET_BLOCK(&this_block_idx);
    claw5_vt->fort_setaux(&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,
                          &maux,aux);
    CLAWPACK5_UNSET_BLOCK();
}

static
void clawpack5_qinit(fclaw2d_global_t *glob,
                     fclaw_patch_t *this_patch,
                     int this_block_idx,
                     int this_patch_idx)
{
    fc2d_clawpack5_vtable_t*  claw5_vt = fc2d_clawpack5_vt(glob);
    
    FCLAW_ASSERT(claw5_vt->fort_qinit != NULL); /* Must initialized */
    int mx,my,mbc,meqn,maux;
    double dx,dy,xlower,ylower;
    double *q, *aux;

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw_clawpatch_soln_data(glob,this_patch,&q,&meqn);
    fclaw_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    /* Call to classic Clawpack 'qinit' routine.  This must be user defined */
    CLAWPACK5_SET_BLOCK(&this_block_idx);
    claw5_vt->fort_qinit(&meqn,&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,q,
                         &maux,aux);
    CLAWPACK5_UNSET_BLOCK();
}

static
void clawpack5_b4step2(fclaw2d_global_t *glob,
                       fclaw_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       double t, double dt)

{
    fc2d_clawpack5_vtable_t*  claw5_vt = fc2d_clawpack5_vt(glob);
    
    if (claw5_vt->fort_b4step2 == NULL)
    {
        return;
    }

    int mx,my,mbc,meqn, maux;
    double xlower,ylower,dx,dy;
    double *aux,*q;

    fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw_clawpatch_soln_data(glob,this_patch,&q,&meqn);
    fclaw_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    CLAWPACK5_SET_BLOCK(&this_block_idx);
    claw5_vt->fort_b4step2(&mbc,&mx,&my,&meqn,q,&xlower,&ylower,
                           &dx,&dy,&t,&dt,&maux,aux);
    CLAWPACK5_UNSET_BLOCK();
}

static
void clawpack5_src2(fclaw2d_global_t *glob,
                    fclaw_patch_t *this_patch,
                    int this_block_idx,
                    int this_patch_idx,
                    double t,
                    double dt)
{
    fc2d_clawpack5_vtable_t*  claw5_vt = fc2d_clawpack5_vt(glob);
    
    if (claw5_vt->fort_src2 == NULL)
    {
        return;
    }

    int mx,my,mbc,meqn, maux;
    double xlower,ylower,dx,dy;
    double *aux,*q;

    fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw_clawpatch_soln_data(glob,this_patch,&q,&meqn);
    fclaw_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    CLAWPACK5_SET_BLOCK(&this_block_idx);
    claw5_vt->fort_src2(&meqn,&mbc,&mx,&my,&xlower,&ylower,
                        &dx,&dy,q,&maux,aux,&t,&dt);
    CLAWPACK5_UNSET_BLOCK();
}

static
void clawpack5_bc2(fclaw2d_global_t *glob,
                   fclaw_patch_t *this_patch,
                   int this_block_idx,
                   int this_patch_idx,
                   double t,
                   double dt,
                   int intersects_phys_bdry[],
                   int time_interp)
{
    fc2d_clawpack5_vtable_t*  claw5_vt = fc2d_clawpack5_vt(glob);

    fc2d_clawpack5_options_t *clawopt = fc2d_clawpack5_get_options(glob);
    
    FCLAW_ASSERT(claw5_vt->fort_bc2 != NULL);

    int mx,my,mbc,meqn, maux;
    double xlower,ylower,dx,dy;
    double *aux,*q;

    fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    int *block_mthbc = clawopt->mthbc;

    /* Set a local copy of mthbc that can be used for a patch. */
    int mthbc[4];
    for(int i = 0; i < 4; i++)
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
    fclaw_clawpatch_timesync_data(glob,this_patch,time_interp,&q,&meqn);

    CLAWPACK5_SET_BLOCK(&this_block_idx);
    claw5_vt->fort_bc2(&meqn,&mbc,&mx,&my,&xlower,&ylower,
                   &dx,&dy,q,&maux,aux,&t,&dt,mthbc);
    CLAWPACK5_UNSET_BLOCK();

}


/* This is called from the single_step callback. and is of type 'flaw_single_step_t' */
static
double clawpack5_step2(fclaw2d_global_t *glob,
                       fclaw_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       double t,
                       double dt)
{
    fc2d_clawpack5_vtable_t*  claw5_vt = fc2d_clawpack5_vt(glob);

    fc2d_clawpack5_options_t* clawpack_options;
    int level;
    double *qold, *aux;
    int mx, my, meqn, maux, mbc;
    double xlower, ylower, dx,dy;

    FCLAW_ASSERT(claw5_vt->fort_rpn2 != NULL);
    FCLAW_ASSERT(claw5_vt->fort_rpt2 != NULL);

    clawpack_options = fc2d_clawpack5_get_options(glob);

    level = this_patch->level;

    fclaw_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    fclaw_clawpatch_save_current_step(glob, this_patch);

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw_clawpatch_soln_data(glob,this_patch,&qold,&meqn);

    int mwaves = clawpack_options->mwaves;

    int maxm = fmax(mx,my);

    double cflgrid = 0.0;

    fclaw2d_clawpatch_registers_t* cr = 
          fclaw2d_clawpatch_get_registers(glob,this_patch);

    /* Evaluate fluxes needed in correction terms */
    const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
    if (fclaw_opt->time_sync && fclaw_opt->flux_correction)
    {
        FCLAW_ASSERT(claw5_vt->fort_rpn2_cons != NULL);
        double *qvec          = FCLAW_ALLOC(double, meqn);
        double *auxvec_center = FCLAW_ALLOC(double, maux);
        double *auxvec_edge   = FCLAW_ALLOC(double, maux);
        double *flux          = FCLAW_ALLOC(double, meqn);     /* f(qr) - f(ql) = amdq+apdq */

#if 1
        CLAWPACK5_TIME_SYNC_STORE_FLUX(&mx,&my,&mbc,&meqn,&maux,
                                       &this_block_idx,&this_patch_idx, &dt,
                                       cr->edgelengths[0], 
                                       cr->edgelengths[1], 
                                       cr->edgelengths[2], 
                                       cr->edgelengths[3],
                                       qold,aux,
                                       cr->edge_fluxes[0],cr->edge_fluxes[1],
                                       cr->edge_fluxes[2],cr->edge_fluxes[3],
                                       claw5_vt->fort_rpn2_cons,
                                       qvec,auxvec_center,auxvec_edge,flux);
#endif
        FCLAW_FREE(qvec);
        FCLAW_FREE(auxvec_center);
        FCLAW_FREE(auxvec_edge);
        FCLAW_FREE(flux);
    }

    int mwork = (maxm+2*mbc)*(12*meqn + (meqn+1)*mwaves + 3*maux + 2);
    double* work = new double[mwork];

    int size = meqn*(mx+2*mbc)*(my+2*mbc);
    double* fp = new double[size];
    double* fm = new double[size];
    double* gp = new double[size];
    double* gm = new double[size];


    int ierror = 0;
    //fc2d_clawpack5_flux2_t flux2 = clawpack_options->use_fwaves ?
    //                                CLAWPACK5_FLUX2FW : CLAWPACK5_FLUX2;
    clawpack5_fort_flux2_t flux2 = CLAWPACK5_FLUX2;

    int* block_corner_count = fclaw2d_patch_block_corner_count(glob,this_patch);
    CLAWPACK5_SET_BLOCK(&this_block_idx);

    CLAWPACK5_STEP2_WRAP(&maxm, &meqn, &maux, &mbc, clawpack_options->method,
                          clawpack_options->mthlim, &clawpack_options->mcapa,
                          &mwaves,&mx, &my, qold, aux, &dx, &dy, &dt, &cflgrid,
                          work, &mwork, &xlower, &ylower, &level,&t, fp, fm, gp, gm,
                          claw5_vt->fort_rpn2, claw5_vt->fort_rpt2,flux2,
                          block_corner_count, &ierror);
    CLAWPACK5_UNSET_BLOCK();

    FCLAW_ASSERT(ierror == 0);

    if (fclaw_opt->time_sync && fclaw_opt->fluctuation_correction)
    {
        CLAWPACK5_TIME_SYNC_ACCUMULATE_WAVES(&mx,&my,&mbc,&meqn, &dt, &dx, 
                                              &dy, &this_patch_idx,
                                              cr->edgelengths[0],
                                              cr->edgelengths[1],
                                              cr->edgelengths[2],
                                              cr->edgelengths[3],
                                              fp,fm,gp,gm,
                                              cr->fp[0],cr->fp[1],
                                              cr->fm[0],cr->fm[1],
                                              cr->gp[0],cr->gp[1],
                                              cr->gm[0],cr->gm[1]);
    }       

    delete [] fp;
    delete [] fm;
    delete [] gp;
    delete [] gm;

    delete [] work;

    return cflgrid;
}

static
double clawpack5_update(fclaw2d_global_t *glob,
                        fclaw_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx,
                        double t,
                        double dt,
                        void* user)
{    
    fc2d_clawpack5_vtable_t*  claw5_vt = fc2d_clawpack5_vt(glob);

    const fc2d_clawpack5_options_t* clawpack_options;
    clawpack_options = fc2d_clawpack5_get_options(glob);

    if (claw5_vt->b4step2 != NULL)
    {
        fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_ADVANCE_B4STEP2]);       
        claw5_vt->b4step2(glob,
                          this_patch,
                          this_block_idx,
                          this_patch_idx,t,dt);
        fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_ADVANCE_B4STEP2]);       
    }

    fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_ADVANCE_STEP2]);       
    double maxcfl = clawpack5_step2(glob,
                                    this_patch,
                                    this_block_idx,
                                    this_patch_idx,t,dt);
    fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_ADVANCE_STEP2]);       

    if (clawpack_options->src_term > 0 && claw5_vt->src2 != NULL)
    {
        claw5_vt->src2(glob,
                       this_patch,
                       this_block_idx,
                       this_patch_idx,t,dt);
    }
    return maxcfl;
}

/* ---------------------------------- Output functions -------------------------------- */

static
void clawpack5_output(fclaw2d_global_t *glob, int iframe)
{
    const fc2d_clawpack5_options_t* clawpack_options;
    clawpack_options = fc2d_clawpack5_get_options(glob);

    if (clawpack_options->ascii_out != 0)
    {
        fclaw_clawpatch_output_ascii(glob,iframe);
    }

    if (clawpack_options->vtk_out != 0)
    {
        fclaw_clawpatch_output_vtk(glob,iframe);
    }

}

/* ---------------------------------- Virtual table  ------------------------------------- */

static
fc2d_clawpack5_vtable_t* fc2d_clawpack5_vt_new()
{
    return (fc2d_clawpack5_vtable_t*) FCLAW_ALLOC_ZERO (fc2d_clawpack5_vtable_t, 1);
}

static
void fc2d_clawpack5_vt_destroy(void* vt)
{
    FCLAW_FREE (vt);
}

/* This is called from the user application. */
void fc2d_clawpack5_solver_initialize(fclaw2d_global_t* glob)
{
	fclaw_clawpatch_options_t* clawpatch_opt = fclaw_clawpatch_get_options(glob);
	fc2d_clawpack5_options_t* clawopt = fc2d_clawpack5_get_options(glob);

    clawopt->method[6] = clawpatch_opt->maux;

    if (clawpatch_opt->maux == 0 && clawopt->mcapa > 0)
    {
        fclaw_global_essentialf("clawpack5 : mcapa > 0 but maux == 0.\n");
        exit(FCLAW_EXIT_ERROR);
    }


    int claw_version = 5;
    fclaw2d_clawpatch_vtable_initialize(glob, claw_version);

    fclaw2d_vtable_t*          fclaw_vt = fclaw2d_vt(glob);
    fclaw2d_patch_vtable_t*    patch_vt = fclaw2d_patch_vt(glob);
    fclaw_clawpatch_vtable_t*      clawpatch_vt = fclaw_clawpatch_vt(glob);

    fc2d_clawpack5_vtable_t*   claw5_vt = fc2d_clawpack5_vt_new();

    fclaw_vt->output_frame      = clawpack5_output;
    fclaw_vt->problem_setup     = clawpack5_setprob;    

    /* Default patch functions */
    patch_vt->initialize            = clawpack5_qinit;
    patch_vt->setup                 = clawpack5_setaux;
    patch_vt->physical_bc           = clawpack5_bc2;
    patch_vt->single_step_update    = clawpack5_update;

    /* Conservation updates (based on Clawpack updates) */
    clawpatch_vt->d2->fort_time_sync_f2c         = CLAWPACK5_FORT_TIME_SYNC_F2C;
    clawpatch_vt->d2->fort_time_sync_samesize    = CLAWPACK5_FORT_TIME_SYNC_SAMESIZE;


    claw5_vt->b4step2   = clawpack5_b4step2;
    claw5_vt->src2      = clawpack5_src2;

    /* Required functions  - error if NULL */
    claw5_vt->fort_bc2       = CLAWPACK5_BC2_DEFAULT;
    claw5_vt->fort_qinit     = NULL;
    claw5_vt->fort_rpn2      = NULL;
    claw5_vt->fort_rpt2      = NULL;

    /* Optional functions - call only if non-NULL */
    claw5_vt->fort_setprob   = NULL;
    claw5_vt->fort_setaux    = NULL;
    claw5_vt->fort_b4step2   = NULL;
    claw5_vt->fort_src2      = NULL;

    claw5_vt->is_set = 1;

	FCLAW_ASSERT(fclaw_pointer_map_get(glob->vtables,"fc2d_clawpack5") == NULL);
	fclaw_pointer_map_insert(glob->vtables, "fc2d_clawpack5", claw5_vt, fc2d_clawpack5_vt_destroy);
}


/* ----------------------------- User access to solver functions --------------------------- */


/* These are here in case the user wants to call Clawpack routines directly */

fc2d_clawpack5_vtable_t* fc2d_clawpack5_vt(fclaw2d_global_t* glob)
{
    fc2d_clawpack5_vtable_t* claw5_vt = (fc2d_clawpack5_vtable_t*) 
	   							fclaw_pointer_map_get(glob->vtables, "fc2d_clawpack5");
	FCLAW_ASSERT(claw5_vt != NULL);
	FCLAW_ASSERT(claw5_vt->is_set != 0);
	return claw5_vt;
}


void fc2d_clawpack5_setprob(fclaw2d_global_t *glob)
{
    clawpack5_setprob(glob);
}

void fc2d_clawpack5_setaux(fclaw2d_global_t *glob,
                            fclaw_patch_t *this_patch,
                            int this_block_idx,
                            int this_patch_idx)
{
    clawpack5_setaux(glob,this_patch,this_block_idx,this_patch_idx);
}


/* This should only be called when a new fclaw_clawpatch_t is created. */
void fc2d_clawpack5_set_capacity(fclaw2d_global_t *glob,
                                  fclaw_patch_t *this_patch,
                                  int this_block_idx,
                                  int this_patch_idx)
{
    int mx,my,mbc,maux,mcapa;
    double dx,dy,xlower,ylower;
    double *aux, *area;
    fc2d_clawpack5_options_t *clawopt;

    clawopt = fc2d_clawpack5_get_options(glob);
    mcapa = clawopt->mcapa;

    fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    area = fclaw2d_clawpatch_get_area(glob,this_patch);

    fclaw_clawpatch_aux_data(glob,this_patch,&aux,&maux);
    FCLAW_ASSERT(maux >= mcapa && mcapa > 0);

    CLAWPACK5_SET_CAPACITY(&mx,&my,&mbc,&dx,&dy,area,&mcapa,
                            &maux,aux);
}

void fc2d_clawpack5_qinit(fclaw2d_global_t *glob,
                           fclaw_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx)
{
    clawpack5_qinit(glob,this_patch,this_block_idx,this_patch_idx);
}

void fc2d_clawpack5_b4step2(fclaw2d_global_t* glob,
                             fclaw_patch_t *this_patch,
                             int this_block_idx,
                             int this_patch_idx,
                             double t,
                             double dt)
{
    clawpack5_b4step2(glob,this_patch,this_block_idx,this_patch_idx,t,dt);
}

void fc2d_clawpack5_bc2(fclaw2d_global_t *glob,
                         fclaw_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx,
                         double t,
                         double dt,
                         int intersects_bc[],
                         int time_interp)
{
    clawpack5_bc2(glob,this_patch,this_block_idx,this_block_idx,t,dt,
                   intersects_bc,time_interp);
}

void fc2d_clawpack5_src2(fclaw2d_global_t* glob,
                          fclaw_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt)
{
    clawpack5_src2(glob,this_patch,this_block_idx,this_block_idx,t,dt);
}






