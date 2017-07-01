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

#include "fc2d_geoclaw.h"
#include "fc2d_geoclaw_options.h"

#include <fclaw2d_clawpatch.hpp>
#include <fclaw2d_clawpatch.h>

/* Basic forestclaw functions */
#include <fclaw2d_convenience.h>  /* Needed to get search function for gauges */
#include "fclaw2d_options.h"
#include <fclaw2d_global.h>
#include <fclaw2d_vtable.h>
#include <fclaw2d_diagnostics.h>
#include <fclaw2d_defs.h>
#include <fclaw2d_transform.h>

/* Some mapping functions */
#include <fclaw2d_map_brick.h>
#include <fclaw2d_map.h>
#include <fclaw2d_map_query.h>


/* Needed for debugging */
#include "types.h"

struct region_type region_type_for_debug;

static fc2d_geoclaw_vtable_t s_geoclaw_vt;

static
fc2d_geoclaw_vtable_t* fc2d_geoclaw_vt_init()
{
    FCLAW_ASSERT(s_geoclaw_vt.is_set == 0);
    return &s_geoclaw_vt;
}

fc2d_geoclaw_vtable_t* fc2d_geoclaw_vt()
{
    FCLAW_ASSERT(s_geoclaw_vt.is_set != 0);
    return &s_geoclaw_vt;
}

typedef struct fc2d_geoclaw_gauge_info
{
    sc_array_t *block_offsets;
    sc_array_t *coordinates;
} fc2d_geoclaw_gauge_info_t;

static fc2d_geoclaw_gauge_info_t gauge_info;

typedef struct fc2d_geoclaw_gauge_acc
{
    int num_gauges;
    int is_latest_domain;
    geoclaw_gauge_t *gauges;
} fc2d_geoclaw_gauge_acc_t;



/* This is provided as a convencience to the user, and is
   called by the user app */
void fc2d_geoclaw_vtable_initialize()
{
    fclaw2d_clawpatch_vtable_initialize();

    fclaw2d_vtable_t*                fclaw_vt = fclaw2d_vt();
    fclaw2d_patch_vtable_t*          patch_vt = fclaw2d_patch_vt();
    fclaw2d_clawpatch_vtable_t*  clawpatch_vt = fclaw2d_clawpatch_vt();
    fclaw2d_diagnostics_vtable_t *    diag_vt = fclaw2d_diagnostics_vt();

    fc2d_geoclaw_vtable_t*  geoclaw_vt = fc2d_geoclaw_vt_init();

    /* ForestClaw virtual tables */
    fclaw_vt->problem_setup      = &fc2d_geoclaw_setprob;  
    fclaw_vt->output_frame       = &fc2d_geoclaw_output;
    fclaw_vt->after_regrid       = &fc2d_geoclaw_after_regrid;  /* handles gauges for new domain */

    /* Set basic patch operations */
    patch_vt->setup              = &fc2d_geoclaw_patch_setup;
    patch_vt->initialize         = &fc2d_geoclaw_qinit;
    patch_vt->physical_bc        = &fc2d_geoclaw_bc2;
    patch_vt->single_step_update = &fc2d_geoclaw_update;  /* Includes b4step2 and src2 */

    /* Regridding */
    patch_vt->tag4refinement     = &fc2d_geoclaw_patch_tag4refinement;
    patch_vt->tag4coarsening     = &fc2d_geoclaw_patch_tag4coarsening;
    patch_vt->interpolate2fine   = &fc2d_geoclaw_interpolate2fine;
    patch_vt->average2coarse     = &fc2d_geoclaw_average2coarse;

    /* Ghost filling */
    patch_vt->average_face       = &fc2d_geoclaw_average_face;
    patch_vt->interpolate_face   = &fc2d_geoclaw_interpolate_face;      
    patch_vt->average_corner     = &fc2d_geoclaw_average_corner;
    patch_vt->interpolate_corner = &fc2d_geoclaw_interpolate_corner;

    patch_vt->remote_ghost_setup = &fc2d_geoclaw_patch_setup;
  
    /* ClawPatch functions (mostly Fortran functions which implement details of the 
       above patch functions */
    clawpatch_vt->fort_ghostpack_extra    = &fc2d_geoclaw_ghostpack_aux;
    clawpatch_vt->fort_compute_error_norm = &FC2D_GEOCLAW_FORT_COMPUTE_ERROR_NORM;
    clawpatch_vt->fort_compute_patch_area = &FC2D_GEOCLAW_FORT_COMPUTE_PATCH_AREA;
    clawpatch_vt->fort_conservation_check = &FC2D_GEOCLAW_FORT_CONSERVATION_CHECK;
    clawpatch_vt->fort_copy_face          = &FC2D_GEOCLAW_FORT_COPY_FACE;
    clawpatch_vt->fort_copy_corner        = &FC2D_GEOCLAW_FORT_COPY_CORNER;
    clawpatch_vt->fort_ghostpack_qarea    = &FC2D_GEOCLAW_FORT_GHOSTPACK_QAREA;
    clawpatch_vt->fort_timeinterp         = &FC2D_GEOCLAW_FORT_TIMEINTERP;

    geoclaw_vt->setprob          = NULL;                   
    geoclaw_vt->setaux           = &GEOCLAW_SETAUX;
    geoclaw_vt->qinit            = &GEOCLAW_QINIT;
    geoclaw_vt->bc2              = &GEOCLAW_BC2;
    geoclaw_vt->b4step2          = &GEOCLAW_B4STEP2;
    geoclaw_vt->src2             = &GEOCLAW_SRC2;
    geoclaw_vt->rpn2             = &GEOCLAW_RPN2;
    geoclaw_vt->rpt2             = &GEOCLAW_RPT2;

    /* Update gauges */
    diag_vt->solver_init_diagnostics     = &fc2d_geoclaw_gauge_initialize;
    diag_vt->solver_compute_diagnostics  = &fc2d_geoclaw_gauge_update;
    diag_vt->solver_finalize_diagnostics = &fc2d_geoclaw_gauge_finalize;

    geoclaw_vt->is_set = 1;
}

/* -----------------------------------------------------------
   Public interface to routines in this file
   ----------------------------------------------------------- */
void fc2d_geoclaw_aux_data(fclaw2d_global_t* glob,
                           fclaw2d_patch_t *this_patch,
                           double **aux, int* maux)
{
    fclaw2d_clawpatch_aux_data(glob, this_patch, aux, maux);
}

void fc2d_geoclaw_setup(fclaw2d_global_t *glob)
{
    const fclaw_options_t* gparms = fclaw2d_get_options(glob);
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
    const fc2d_geoclaw_options_t *geoclaw_options = fc2d_geoclaw_get_options(glob);

    GEOCLAW_SET_MODULES(&geoclaw_options->mwaves, &geoclaw_options->mcapa,
                        &clawpatch_opt->meqn, &clawpatch_opt->maux,
                        geoclaw_options->mthlim, geoclaw_options->method,
                        &gparms->ax, &gparms->bx, &gparms->ay, &gparms->by);
}

void fc2d_geoclaw_gauge_initialize(fclaw2d_global_t* glob, void** acc)
{
    fc2d_geoclaw_gauge_acc_t* gauge_acc;
    gauge_acc = FCLAW_ALLOC(fc2d_geoclaw_gauge_acc_t,1);
    *acc = gauge_acc;

    const fclaw_options_t * gparms = fclaw2d_get_options(glob);
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    /* --------------------------------------------------------
       Read gauges files 'gauges.data' to get number of gauges
       -------------------------------------------------------- */
    char fname[] = "gauges.data";
    int num = GEOCLAW_GAUGES_GETNUM(fname);
    int restart = 0;

    gauge_acc->num_gauges = num;
    gauge_acc->gauges = FCLAW_ALLOC(geoclaw_gauge_t,num);
    
    if (num == 0)
    {
        return;
    }

    /* Read gauges file for the locations, etc. of all gauges */
    GEOCLAW_GAUGES_INIT(&restart, &clawpatch_opt->meqn, &num,  gauge_acc->gauges, fname);

    /* -----------------------------------------------------
       Open gauge files and add header information
       ----------------------------------------------------- */
    char filename[14];    /* gaugexxxxx.txt */
    FILE *fp;

    for (int i = 0; i < gauge_acc->num_gauges; ++i)
    {
        geoclaw_gauge_t g = gauge_acc->gauges[i];
        sprintf(filename,"gauge%05d.txt",g.num);
        fp = fopen(filename, "w");
        fprintf(fp, "# gauge_id= %5d location=( %15.7e %15.7e ) num_eqn= %2d\n",
                g.num, g.xc, g.yc, clawpatch_opt->meqn+1);
        fprintf(fp, "# Columns: level time h    hu    hv    eta\n");
        fclose(fp);
    }

    /* -----------------------------------------------------
       Set up block offsets and coordinate list for p4est
       search function
       ----------------------------------------------------- */

    fclaw2d_map_context_t* cont = glob->cont;

    int is_brick = FCLAW2D_MAP_IS_BRICK(&cont);

    gauge_info.block_offsets = sc_array_new_size(sizeof(int), glob->domain->num_blocks+1);
    gauge_info.coordinates = sc_array_new_size(2*sizeof(double), num);

    int *block_offsets = (int*) sc_array_index_int(gauge_info.block_offsets, 0);
    double *coordinates = (double*) sc_array_index_int(gauge_info.coordinates, 0);

    geoclaw_gauge_t *gauges = gauge_acc->gauges;

    if (is_brick)
    {
        /* We don't know how the blocks are arranged in the brick domain
           so we reverse engineer this information
        */
        int nb,mi,mj;
        double x,y;
        double z;
        double xll,yll;
        double xur,yur;
        int number_of_gauges_set;


        double x0,y0,x1,y1;
        x0 = 0;
        y0 = 0;
        x1 = 1;
        y1 = 1;

        FCLAW2D_MAP_BRICK_GET_DIM(&cont,&mi,&mj);

        number_of_gauges_set = 0;
        block_offsets[0] = 0;
        for (nb = 0; nb < glob->domain->num_blocks; ++nb)
        {
            /* Scale to [0,1]x[0,1], based on blockno */
            fclaw2d_map_c2m_nomap_brick(cont,nb,x0,y0,&xll,&yll,&z);
            fclaw2d_map_c2m_nomap_brick(cont,nb,x1,y1,&xur,&yur,&z);
            for(int i = 0; i < num; i++)
            {
                /* Map gauge to global [0,1]x[0,1] space */
                x = (gauges[i].xc - gparms->ax)/(gparms->bx-gparms->ax);
                y = (gauges[i].yc - gparms->ay)/(gparms->by-gparms->ay);
                if (xll <= x && x <= xur && yll <= y && y <= yur)
                {
                    int ng = number_of_gauges_set;
                    gauges[i].blockno = nb;
                    coordinates[2*ng] = mi*(x - xll);
                    coordinates[2*ng+1] = mj*(y - yll);
                    gauges[i].location_in_results = ng;
                    number_of_gauges_set++;
                }
            }
            block_offsets[nb+1] = number_of_gauges_set;
        }
    }
    else
    {
        for(int i = 0; i < num; i++)
        {
            gauges[i].blockno = 0;
            gauges[i].location_in_results = i;
        }

        block_offsets[0] = 0;
        block_offsets[1] = num;

        for (int i = 0; i < num; ++i)
        {
            coordinates[2*i] = (gauges[i].xc - gparms->ax)/(gparms->bx-gparms->ax);
            coordinates[2*i+1] = (gauges[i].yc - gparms->ay)/(gparms->by-gparms->ay);
        }
    }
}

void fc2d_geoclaw_gauge_update(fclaw2d_global_t *glob, void* solver_acc)
{
    fc2d_geoclaw_gauge_acc_t* gauge_acc = (fc2d_geoclaw_gauge_acc_t*) solver_acc;

    const double tcurr = glob->curr_time;
    int mx,my,mbc,meqn,maux;
    double dx,dy,xlower,ylower,eta;
    double *q, *aux;
    double *var;
    char filename[14];  /* gaugexxxxx.txt */
    FILE *fp;

    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    if (gauge_acc->num_gauges == 0)
    {
        return;
    }

    fclaw2d_block_t *block;
    fclaw2d_patch_t *patch;

    var = FCLAW_ALLOC(double, clawpatch_opt->meqn);
    geoclaw_gauge_t *gauges = gauge_acc->gauges;
    geoclaw_gauge_t g;
    for (int i = 0; i < gauge_acc->num_gauges; ++i)
    {
        g = gauges[i];
        block = &glob->domain->blocks[g.blockno];
        if (g.patchno >= 0)
        {
            patch = &block->patches[g.patchno];
            FCLAW_ASSERT(patch != NULL);
            fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                        &xlower,&ylower,&dx,&dy);

            fclaw2d_clawpatch_soln_data(glob,patch,&q,&meqn);
            fc2d_geoclaw_aux_data(glob,patch,&aux,&maux);

            FCLAW_ASSERT(g.xc >= xlower && g.xc <= xlower+mx*dx);
            FCLAW_ASSERT(g.yc >= ylower && g.yc <= ylower+my*dy);
            if (tcurr >= g.t1 && tcurr <= g.t2)
            {
                GEOCLAW_UPDATE_GAUGE(&mx,&my,&mbc,&meqn,&xlower,&ylower,&dx,&dy,
                                      q,&maux,aux,&g.xc,&g.yc,var,&eta);
                sprintf(filename,"gauge%05d.txt",g.num);
                fp = fopen(filename, "a");
                fprintf(fp, "%5d %15.7e %15.7e %15.7e %15.7e %15.7e\n",
                        patch->level,tcurr,var[0],var[1],var[2],eta);
                fclose(fp);
            }
        }
    }
    FCLAW_FREE(var);
}

void fc2d_geoclaw_gauge_finalize(fclaw2d_global_t *glob, void** acc)
{
    fc2d_geoclaw_gauge_acc_t* gauge_acc = *((fc2d_geoclaw_gauge_acc_t**) acc);
    FCLAW_FREE(gauge_acc->gauges);
    FCLAW_FREE(gauge_acc);
    *acc = NULL;
}


void fc2d_geoclaw_finalize(fclaw2d_global_t *glob)
{
    if (gauge_info.block_offsets != NULL)
    {
        sc_array_destroy(gauge_info.block_offsets);
    }
    if (gauge_info.coordinates != NULL)
    {
        sc_array_destroy(gauge_info.coordinates);
    }
}


void fc2d_geoclaw_after_regrid(fclaw2d_global_t *glob)
{
    int i,index;
    fc2d_geoclaw_gauge_acc_t* gauge_acc = (fc2d_geoclaw_gauge_acc_t*) glob->acc->solver_accumulator;

    /* Locate each gauge in the new mesh */
    int num = gauge_acc->num_gauges;

    if (num == 0)
    {
        return;
    }

    sc_array_t *results = sc_array_new_size(sizeof(int), num);
    fclaw2d_domain_search_points(glob->domain, gauge_info.block_offsets,
                                 gauge_info.coordinates, results);

    for (i = 0; i < gauge_acc->num_gauges; ++i)
    {
        index = gauge_acc->gauges[i].location_in_results;
        FCLAW_ASSERT(index >= 0 && index < num);

        /* patchno == -1  : Patch is not on this processor
           patchno >= 0   : Patch number is local patch list.
        */

        int patchno = *((int *) sc_array_index_int(results, index));
        gauge_acc->gauges[i].patchno = patchno;
    }
    sc_array_destroy(results);
}


void fc2d_geoclaw_setprob(fclaw2d_global_t *glob)
{
    fc2d_geoclaw_vtable_t *geoclaw_vt = fc2d_geoclaw_vt();
    if (geoclaw_vt->setprob != NULL)
    {
        geoclaw_vt->setprob();
    }
}


void fc2d_geoclaw_patch_setup(fclaw2d_global_t *glob,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx)
{
    fc2d_geoclaw_setaux(glob,this_patch,this_block_idx,this_patch_idx);
}

void fc2d_geoclaw_ghostpatch_setup(fclaw2d_global_t *glob,
                                   fclaw2d_patch_t *this_patch,
                                   int this_block_idx,
                                   int this_patch_idx)
{
    fclaw2d_clawpatch_options_t* clawpatch_options;
    clawpatch_options = fclaw2d_clawpatch_get_options(glob);

    if (!clawpatch_options->ghost_patch_pack_aux)
    {
        fc2d_geoclaw_setaux(glob,this_patch,this_block_idx,this_patch_idx);
    }
    else
    {
        /* the aux array data has been packed and transferred as MPI messages */
    }
}

/* This should only be called when a new fclaw2d_clawpatch_t is created. */
void fc2d_geoclaw_setaux(fclaw2d_global_t *glob,
                         fclaw2d_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx)
{
    fc2d_geoclaw_vtable_t *geoclaw_vt = fc2d_geoclaw_vt();
    int mx,my,mbc,maux;
    double xlower,ylower,dx,dy;
    double *aux;

    int is_ghost = fclaw2d_patch_is_ghost(this_patch);

    fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int mint = 2*mbc;
    int nghost = mbc;

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);
    GEOCLAW_SET_BLOCK(&this_block_idx);
    geoclaw_vt->setaux(&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,
                      &maux,aux,&is_ghost,&nghost,&mint);
    GEOCLAW_UNSET_BLOCK();
}


void fc2d_geoclaw_qinit(fclaw2d_global_t *glob,
                        fclaw2d_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx)
{
    fc2d_geoclaw_vtable_t *geoclaw_vt = fc2d_geoclaw_vt();

    FCLAW_ASSERT(geoclaw_vt->qinit != NULL); /* Must initialized */
    int mx,my,mbc,meqn,maux;
    double dx,dy,xlower,ylower;
    double *q, *aux;

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);
    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    /* Call to classic Clawpack 'qinit' routine.  This must be user defined */
    GEOCLAW_SET_BLOCK(&this_block_idx);
    geoclaw_vt->qinit(&meqn,&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,q,
                     &maux,aux);
    GEOCLAW_UNSET_BLOCK();
}

void fc2d_geoclaw_b4step2(fclaw2d_global_t *glob,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t, double dt)

{
    fc2d_geoclaw_vtable_t *geoclaw_vt = fc2d_geoclaw_vt();

    FCLAW_ASSERT(geoclaw_vt->b4step2 != NULL);

    int mx,my,mbc,meqn, maux;
    double xlower,ylower,dx,dy;
    double *aux,*q;

    fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);
    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    GEOCLAW_SET_BLOCK(&this_block_idx);
    geoclaw_vt->b4step2(&mbc,&mx,&my,&meqn,q,&xlower,&ylower,
                       &dx,&dy,&t,&dt,&maux,aux);
    GEOCLAW_UNSET_BLOCK();
}

void fc2d_geoclaw_src2(fclaw2d_global_t *glob,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       double t,
                       double dt)
{
    fc2d_geoclaw_vtable_t *geoclaw_vt = fc2d_geoclaw_vt();

    FCLAW_ASSERT(geoclaw_vt->src2 != NULL);

    int mx,my,mbc,meqn, maux;
    double xlower,ylower,dx,dy;
    double *aux,*q;

    fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);
    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    GEOCLAW_SET_BLOCK(&this_block_idx);
    geoclaw_vt->src2(&meqn,&mbc,&mx,&my,&xlower,&ylower,
                    &dx,&dy,q,&maux,aux,&t,&dt);
    GEOCLAW_UNSET_BLOCK();

}



/* Use this to return only the right hand side of the clawpack algorithm */
double fc2d_geoclaw_step2_rhs(fclaw2d_global_t *glob,
                              fclaw2d_patch_t *this_patch,
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


void fc2d_geoclaw_bc2(fclaw2d_global_t *glob,
                      fclaw2d_patch_t *this_patch,
                      int this_block_idx,
                      int this_patch_idx,
                      double t,
                      double dt,
                      int intersects_phys_bdry[],
                      int time_interp)
{
    fc2d_geoclaw_options_t *geo_opt = fc2d_geoclaw_get_options(glob);
    fc2d_geoclaw_vtable_t *geoclaw_vt = fc2d_geoclaw_vt();

    FCLAW_ASSERT(geoclaw_vt->bc2 != NULL);

    int mx,my,mbc,meqn, maux;
    double xlower,ylower,dx,dy;
    double *aux,*q;

    fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    int *block_mthbc = geo_opt->mthbc;

    /* Set a local copy of mthbc that can be used for a patch. */
    int mthbc[NumFaces];
    for(int i = 0; i < NumFaces; i++)
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

    GEOCLAW_SET_BLOCK(&this_block_idx);
    geoclaw_vt->bc2(&meqn,&mbc,&mx,&my,&xlower,&ylower,
                   &dx,&dy,q,&maux,aux,&t,&dt,mthbc);
    GEOCLAW_UNSET_BLOCK();

}


/* This is called from the single_step callback. and is of type 'flaw_single_step_t' */
double fc2d_geoclaw_step2(fclaw2d_global_t *glob,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt)
{
    fc2d_geoclaw_vtable_t *geoclaw_vt = fc2d_geoclaw_vt();
    fc2d_geoclaw_options_t* geoclaw_options;
    int level;
    double *qold, *aux;
    int mx, my, meqn, maux, mbc;
    double xlower, ylower, dx,dy;

    FCLAW_ASSERT(geoclaw_vt->rpn2 != NULL);
    FCLAW_ASSERT(geoclaw_vt->rpt2 != NULL);

    geoclaw_options = fc2d_geoclaw_get_options(glob);

    level = this_patch->level;

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    fclaw2d_clawpatch_save_current_step(glob, this_patch);

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&qold,&meqn);

    int mwaves = geoclaw_options->mwaves;

    int maxm = fmax(mx,my);

    double cflgrid = 0.0;

    /* This is no longer needed */
    int mwork = (maxm+2*mbc)*(12*meqn + (meqn+1)*mwaves + 3*maux + 2);
    double* work = new double[mwork];

    int size = meqn*(mx+2*mbc)*(my+2*mbc);
    double* fp = new double[size];
    double* fm = new double[size];
    double* gp = new double[size];
    double* gm = new double[size];

    int* block_corner_count = fclaw2d_patch_block_corner_count(glob,this_patch);

    GEOCLAW_STEP2_WRAP(&maxm, &meqn, &maux, &mbc, geoclaw_options->method,
                       geoclaw_options->mthlim, &geoclaw_options->mcapa,
                       &mwaves,&mx, &my, qold, aux, &dx, &dy, &dt, &cflgrid,
                       work, &mwork, &xlower, &ylower, &level,&t, fp, fm, gp, gm,
                       geoclaw_vt->rpn2, geoclaw_vt->rpt2,
                       block_corner_count);

    /* Accumulate fluxes needed for conservative fix-up */
    if (geoclaw_vt->fluxfun != NULL)
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

double fc2d_geoclaw_update(fclaw2d_global_t *glob,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx,
                           double t,
                           double dt)
{
    fc2d_geoclaw_vtable_t *geoclaw_vt = fc2d_geoclaw_vt();

    const fc2d_geoclaw_options_t* geoclaw_options;
    geoclaw_options = fc2d_geoclaw_get_options(glob);

    GEOCLAW_TOPO_UPDATE(&t);

    if (geoclaw_vt->b4step2 != NULL)
    {
        fc2d_geoclaw_b4step2(glob,
                             this_patch,
                             this_block_idx,
                             this_patch_idx,t,dt);
    }

    double maxcfl = fc2d_geoclaw_step2(glob,
                                       this_patch,
                                       this_block_idx,
                                       this_patch_idx,t,dt);

    if (geoclaw_options->src_term > 0)
    {
        fc2d_geoclaw_src2(glob,
                          this_patch,
                          this_block_idx,
                          this_patch_idx,t,dt);
    }

    return maxcfl;
}

int fc2d_geoclaw_patch_tag4refinement(fclaw2d_global_t *glob,
                                      fclaw2d_patch_t *this_patch,
                                      int blockno, int this_patch_idx,
                                      int initflag)
{
    int mx,my,mbc,meqn,maux,mbathy;
    double xlower,ylower,dx,dy, t;
    double *q, *aux;
    int tag_patch;
    int level,maxlevel;
    const fclaw_options_t * gparms = fclaw2d_get_options(glob);
    const fc2d_geoclaw_options_t* geoclaw_options;

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);
    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);
    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    level = this_patch->level;
    maxlevel = gparms->maxlevel;
    t = glob->curr_time;

    geoclaw_options = fc2d_geoclaw_get_options(glob);
    mbathy = geoclaw_options->mbathy;

    tag_patch = 0;

    FC2D_GEOCLAW_FORT_TAG4REFINEMENT(&mx,&my,&mbc,&meqn,&maux,&xlower,&ylower,
                                     &dx,&dy,&t,&blockno,q,aux,&mbathy,&level,&maxlevel,
                                     &initflag,&tag_patch);

    return tag_patch;
}


int fc2d_geoclaw_patch_tag4coarsening(fclaw2d_global_t *glob,
                                      fclaw2d_patch_t *fine_patches,
                                      int blockno, int patchno)

{
    const fclaw_options_t *amropt;
    fc2d_geoclaw_options_t* geoclaw_options;

    int mx,my,mbc,meqn,maux,mbathy;
    double xlower,ylower,dx,dy,t;
    double *q[4], *aux[4];
    int tag_patch,igrid;
    //double coarsen_threshold;
    int level, maxlevel;

    fclaw2d_patch_t *patch0;
    patch0 = &fine_patches[0];


    amropt = fclaw2d_get_options(glob);
    geoclaw_options = fc2d_geoclaw_get_options(glob);
    mbathy = geoclaw_options->mbathy;
    //coarsen_threshold = amropt->coarsen_threshold;
    t = glob->curr_time;

    fclaw2d_clawpatch_grid_data(glob,patch0,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    for (igrid = 0; igrid < 4; igrid++)
    {
        fclaw2d_clawpatch_soln_data(glob,&fine_patches[igrid],&q[igrid],&meqn);
        fclaw2d_clawpatch_aux_data(glob,&fine_patches[igrid],&aux[igrid],&maux);
    }

    level = fine_patches[0].level;
    maxlevel = amropt->maxlevel;

    tag_patch = 0;
    FC2D_GEOCLAW_FORT_TAG4COARSENING(&patchno,&mx,&my,&mbc,&meqn,&maux,&xlower,&ylower,&dx,&dy,
                                     &t,q[0],q[1],q[2],q[3],aux[0],aux[1],aux[2],aux[3],
                                     &mbathy,&level,&maxlevel, &geoclaw_options->dry_tolerance_c,
                                     &geoclaw_options->wave_tolerance_c,
                                     &geoclaw_options->speed_tolerance_entries_c,
                                     geoclaw_options->speed_tolerance_c, &tag_patch);
    return tag_patch;

}

void fc2d_geoclaw_interpolate2fine(fclaw2d_global_t *glob,
                                   fclaw2d_patch_t *coarse_patch,
                                   fclaw2d_patch_t *fine_patches,
                                   int this_blockno, int coarse_patchno,
                                   int fine0_patchno)

{
    int mx,my,mbc,meqn,maux,refratio,p4est_refineFactor,mbathy;
    double *qcoarse,*qfine,*auxcoarse,*auxfine;
    int igrid;

    const fclaw_options_t* gparms;
    const fc2d_geoclaw_options_t*  geoclaw_options;
    const fclaw2d_clawpatch_options_t *clawpatch_opt;
    fclaw2d_patch_t* fine_patch;

    gparms = fclaw2d_get_options(glob);
    geoclaw_options = fc2d_geoclaw_get_options(glob);
    clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    mbathy = geoclaw_options->mbathy;
    mx = clawpatch_opt->mx;
    my = clawpatch_opt->my;
    mbc = clawpatch_opt->mbc;
    refratio = gparms->refratio;
    p4est_refineFactor = RefineFactor;

    // fclaw2d_clawpatch_metric_data(domain,coarse_patch,&xp,&yp,&zp,
    //                               &xd,&yd,&zd,&areacoarse);

    fclaw2d_clawpatch_soln_data(glob,coarse_patch,&qcoarse,&meqn);
    fclaw2d_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);

    /* Loop over four siblings (z-ordering) */
    for (igrid = 0; igrid < 4; igrid++)
    {
        fine_patch = &fine_patches[igrid];

        fclaw2d_clawpatch_soln_data(glob,fine_patch,&qfine,&meqn);
        fclaw2d_clawpatch_aux_data(glob,fine_patch,&auxfine,&maux);

        FC2D_GEOCLAW_FORT_INTERPOLATE2FINE(&mx,&my,&mbc,&meqn,qcoarse,qfine,
                                           &maux,auxcoarse,auxfine,&mbathy,
                                           &p4est_refineFactor,&refratio,
                                           &igrid);
    }
}

void fc2d_geoclaw_average2coarse(fclaw2d_global_t *glob,
                                 fclaw2d_patch_t *fine_patches,
                                 fclaw2d_patch_t *coarse_patch,
                                 int blockno, int fine0_patchno,
                                 int coarse_patchno)

{
    int mx,my,mbc,meqn,maux,refratio,p4est_refineFactor,mbathy;
    double *qcoarse,*qfine,*auxcoarse,*auxfine;
    // double *areacoarse,*areafine;
    // double *xp,*yp,*zp,*xd,*yd,*zd;
    int igrid,mcapa;

    const fclaw_options_t* gparms;
    const fc2d_geoclaw_options_t*  geoclaw_options;
    const fclaw2d_clawpatch_options_t *clawpatch_opt;

    fclaw2d_patch_t* fine_patch;

    gparms = fclaw2d_get_options(glob);
    geoclaw_options = fc2d_geoclaw_get_options(glob);
    clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    mcapa = geoclaw_options->mcapa;
    mx  = clawpatch_opt->mx;
    my  = clawpatch_opt->my;
    mbc = clawpatch_opt->mbc;
    refratio = gparms->refratio;
    p4est_refineFactor = RefineFactor;
    mbathy = geoclaw_options->mbathy;

    // fclaw2d_clawpatch_metric_data(domain,coarse_patch,&xp,&yp,&zp,
    //                               &xd,&yd,&zd,&areacoarse);

    fclaw2d_clawpatch_soln_data(glob,coarse_patch,&qcoarse,&meqn);
    fclaw2d_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);

    /* Loop over four siblings (z-ordering) */
    for (igrid = 0; igrid < 4; igrid++)
    {
        fine_patch = &fine_patches[igrid];

        fclaw2d_clawpatch_soln_data(glob,fine_patch,&qfine,&meqn);
        fclaw2d_clawpatch_aux_data(glob,fine_patch,&auxfine,&maux);

        // if (gparms->manifold)
        // {
        //     fclaw2d_clawpatch_metric_data(domain,fine_patch,&xp,&yp,&zp,
        //                                   &xd,&yd,&zd,&areafine);
        // }
        FC2D_GEOCLAW_FORT_AVERAGE2COARSE(&mx,&my,&mbc,&meqn,qcoarse,qfine,
                                         &maux,auxcoarse,auxfine,&mcapa,&mbathy,
                                         &p4est_refineFactor,&refratio,
                                         &igrid);

    }
}

void fc2d_geoclaw_average_face(fclaw2d_global_t *glob,
                               fclaw2d_patch_t *coarse_patch,
                               fclaw2d_patch_t *fine_patch,
                               int idir,
                               int iface_coarse,
                               int p4est_refineFactor,
                               int refratio,
                               int time_interp,
                               int igrid,
                               fclaw2d_transform_data_t* transform_data)
{
    int meqn,mx,my,mbc;
    int mcapa, mbathy, maux;
    double *qcoarse, *qfine;
    double *auxcoarse, *auxfine;
    
    const fclaw_options_t *gparms = fclaw2d_get_options(glob);
    const fc2d_geoclaw_options_t *geoclaw_options;
    const fclaw2d_clawpatch_options_t *clawpatch_opt;

    geoclaw_options = fc2d_geoclaw_get_options(glob);
    clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    fclaw2d_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
    qfine = fclaw2d_clawpatch_get_q(glob,fine_patch);

    /* These will be empty for non-manifolds cases */
    fclaw2d_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);
    fclaw2d_clawpatch_aux_data(glob,fine_patch,&auxfine,&maux);

    mcapa = geoclaw_options->mcapa;
    mbathy = geoclaw_options->mbathy;

    mx = clawpatch_opt->mx;
    my = clawpatch_opt->my;
    mbc = clawpatch_opt->mbc;

    int manifold = gparms->manifold;
    FC2D_GEOCLAW_FORT_AVERAGE_FACE(&mx,&my,&mbc,&meqn,qcoarse,qfine,
                                   &maux,auxcoarse,auxfine,&mcapa,&mbathy,
                                   &idir,&iface_coarse,&p4est_refineFactor,&refratio,
                                   &igrid,&manifold,&transform_data);
}

void fc2d_geoclaw_interpolate_face(fclaw2d_global_t *glob,
                                   fclaw2d_patch_t *coarse_patch,
                                   fclaw2d_patch_t *fine_patch,
                                   int idir,
                                   int iside,
                                   int p4est_refineFactor,
                                   int refratio,
                                   int time_interp,
                                   int igrid,
                                   fclaw2d_transform_data_t* transform_data)
{

    int meqn,mx,my,mbc,maux,mbathy;
    double *qcoarse, *qfine;
    double *auxcoarse, *auxfine;
    
    const fc2d_geoclaw_options_t *geoclaw_options;
    const fclaw2d_clawpatch_options_t *clawpatch_opt;
    
    geoclaw_options = fc2d_geoclaw_get_options(glob);
    clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    fclaw2d_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
    qfine = fclaw2d_clawpatch_get_q(glob,fine_patch);

    fclaw2d_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);
    fclaw2d_clawpatch_aux_data(glob,fine_patch,&auxfine,&maux);

    mbathy = geoclaw_options->mbathy;
    mx = clawpatch_opt->mx;
    my = clawpatch_opt->my;
    mbc = clawpatch_opt->mbc;

    FC2D_GEOCLAW_FORT_INTERPOLATE_FACE(&mx,&my,&mbc,&meqn,qcoarse,qfine,&maux,
                                       auxcoarse,auxfine,&mbathy,
                                       &idir, &iside,
                                       &p4est_refineFactor,&refratio,&igrid,
                                       &transform_data);
}

void fc2d_geoclaw_average_corner(fclaw2d_global_t *glob,
                                 fclaw2d_patch_t *coarse_patch,
                                 fclaw2d_patch_t *fine_patch,
                                 int coarse_corner,
                                 int refratio,
                                 int time_interp,
                                 fclaw2d_transform_data_t* transform_data)
{
    int meqn,mx,my,mbc;
    int maux,mbathy,mcapa;
    double *qcoarse, *qfine;
    double *auxcoarse, *auxfine;
    
    const fclaw_options_t *gparms = fclaw2d_get_options(glob);
    const fc2d_geoclaw_options_t *geoclaw_options;
    const fclaw2d_clawpatch_options_t *clawpatch_opt;

    geoclaw_options = fc2d_geoclaw_get_options(glob);
    clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    fclaw2d_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
    qfine = fclaw2d_clawpatch_get_q(glob,fine_patch);

    fclaw2d_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);
    fclaw2d_clawpatch_aux_data(glob,fine_patch,&auxfine,&maux);

    mx = clawpatch_opt->mx;
    my = clawpatch_opt->my;
    mbc = clawpatch_opt->mbc;

    mcapa=geoclaw_options->mcapa;
    mbathy=geoclaw_options->mbathy;

    int manifold = gparms->manifold;
    FC2D_GEOCLAW_FORT_AVERAGE_CORNER(&mx,&my,&mbc,&meqn,&refratio,
                                     qcoarse,qfine,&maux,auxcoarse,auxfine,
                                     &mcapa,&mbathy,&manifold,&coarse_corner,
                                     &transform_data);

}

void fc2d_geoclaw_interpolate_corner(fclaw2d_global_t* glob,
                                     fclaw2d_patch_t* coarse_patch,
                                     fclaw2d_patch_t* fine_patch,
                                     int coarse_corner,
                                     int refratio,
                                     int time_interp,
                                     fclaw2d_transform_data_t* transform_data)

{
    int meqn,mx,my,mbc;
    int mbathy,maux;
    double *qcoarse, *qfine;
    double *auxcoarse, *auxfine;

    const fclaw2d_clawpatch_options_t *clawpatch_opt;
    const fc2d_geoclaw_options_t *geoclaw_options;

    geoclaw_options = fc2d_geoclaw_get_options(glob);
    clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    mx = clawpatch_opt->mx;
    my = clawpatch_opt->my;
    mbc = clawpatch_opt->mbc;
    meqn = clawpatch_opt->meqn;

    mbathy = geoclaw_options->mbathy;

    fclaw2d_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
    qfine = fclaw2d_clawpatch_get_q(glob,fine_patch);

    fclaw2d_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);
    fclaw2d_clawpatch_aux_data(glob,fine_patch,&auxfine,&maux);

    FC2D_GEOCLAW_FORT_INTERPOLATE_CORNER(&mx,&my,&mbc,&meqn,
                                         &refratio,qcoarse,qfine,&maux,
                                         auxcoarse,auxfine,&mbathy,
                                         &coarse_corner,&transform_data);

}

void fc2d_geoclaw_output(fclaw2d_global_t *glob, int iframe)
{
    const fc2d_geoclaw_options_t* geo_opt;
    geo_opt = fc2d_geoclaw_get_options(glob);

    if (geo_opt->ascii_out != 0)
    {
        fc2d_geoclaw_output_ascii(glob,iframe);        
    }
}


void fc2d_geoclaw_ghostpack_aux(fclaw2d_global_t *glob,
                                fclaw2d_patch_t *this_patch,
                                int mint,
                                double *auxpack,
                                int auxsize, int packmode,
                                int* ierror)
{
    int mx,my,mbc,maux;
    double xlower,ylower,dx,dy;
    double *aux;

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);
    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);
    FC2D_GEOCLAW_FORT_GHOSTPACKAUX(&mx,&my,&mbc,&maux,
                                   &mint,aux,auxpack,&auxsize,
                                   &packmode,ierror);
}
