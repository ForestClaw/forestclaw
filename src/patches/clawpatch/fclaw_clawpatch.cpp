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

#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch.hpp>

#include <fclaw_clawpatch_diagnostics.h>
#include <fclaw_clawpatch_options.h>
#include <fclaw_clawpatch_output_ascii.h> 
#include <fclaw_clawpatch_output_vtk.h>
#include <fclaw_clawpatch_pillow.h>  

#include <fclaw2d_clawpatch_fort.h>
#include <fclaw2d_clawpatch_conservation.h>
#include <fclaw2d_clawpatch_conservation_fort.h>
#include <fclaw2d_clawpatch_transform.h>

#include <fclaw2d_clawpatch46_fort.h>
#include <fclaw2d_clawpatch5_fort.h>

#include <fclaw2d_metric.h>
#include <fclaw2d_metric.hpp>

#include <fclaw3d_clawpatch_fort.h>
#include <fclaw2d_clawpatch_transform.h>

#include <fclaw3d_clawpatch_transform.h>

#include <fclaw3d_clawpatch46_fort.h>

#include <fclaw3d_metric.h>
#include <fclaw3d_metric.hpp>


#include <fclaw_patch.h>  /* Needed to get enum for build modes */

#include <fclaw2d_defs.h>
#include <fclaw_global.h>
#include <fclaw_vtable.h>
#include <fclaw_options.h>

#include <fclaw_timeinterp.h>
#include <fclaw_diagnostics.h>

#include <fclaw2d_map_query.h>

#include <fclaw_pointer_map.h>



/* ------------------------------- Static function defs ------------------------------- */

/* Added to turn off time_interp */
static int fill_ghost(fclaw_global_t* glob, int time_interp)
{
    const fclaw_options_t* fclaw_opt = fclaw_get_options(glob);
    if (fclaw_opt->timeinterp2fillghost)
        /* This will always fill ghost cells using data from "qsync", which is either 
           coarse grid data at time step, or time interpolated data */
        return 1;
    else
        /* Only fill ghost cells with neighboring data if not doing time interpolation. 
          If doing time interpolation, then fine grid data at intermediate time steps 
          will be filled in some other way (.e.g. by advancing the solution in ghost 
          cells.) */
        return !time_interp;
}


/* Store virtual table for retrieval from anywhere */

static
fclaw_clawpatch_t* get_clawpatch(fclaw_patch_t *patch)
{
    fclaw_clawpatch_t *cp = (fclaw_clawpatch_t*) 
                     fclaw_patch_get_user_patch(patch);
    return cp;
}

/* Needed for virtual patch function so that the metric class can be independent 
of a clawpatch object */
static 
void* clawpatch_get_metric_patch_2d(fclaw_patch_t* patch)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    return cp->d2->mp;
}

static 
void* clawpatch_get_metric_patch_3d(fclaw_patch_t* patch)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    return cp->d3->mp;
}

static
fclaw2d_metric_patch_t* get_metric_patch_2d(fclaw_patch_t *patch)
{
    return (fclaw2d_metric_patch_t*) clawpatch_get_metric_patch_2d(patch);
}

static
fclaw3d_metric_patch_t* get_metric_patch_3d(fclaw_patch_t *patch)
{
    return (fclaw3d_metric_patch_t*) clawpatch_get_metric_patch_3d(patch);
}

/* Return a pointer to either time interpolated data or regular grid data */
static 
double* q_time_sync(fclaw_patch_t* patch, int time_interp)
{
    fclaw_clawpatch_t* cp = get_clawpatch(patch);
    if (time_interp)
        return cp->griddata_time_interpolated.dataPtr();
    else
        return cp->griddata.dataPtr();
}

static 
double* clawpatch_get_area(struct fclaw_global* glob,
                           fclaw_patch_t* patch)
{
    return fclaw2d_metric_patch_get_area(glob, patch);
}

static 
double* clawpatch_get_volume(struct fclaw_global* glob,
                           fclaw_patch_t* patch)
{
    return fclaw3d_metric_patch_get_volume(glob, patch);
}

/* ----------------------------- Creating/deleting patches ---------------------------- */

static
void* clawpatch_new(int dim)
{
    fclaw_clawpatch_t *cp = new fclaw_clawpatch_t;    
    if(dim == 2)
    {
        cp->patch_dim = 2;
        cp->d2 = new fclaw_clawpatch_2d_t;
        /* This patch will only be defined if we are on a manifold. */
        cp->d2->mp = fclaw2d_metric_patch_new();
    }
    else 
    {
        cp->patch_dim = 3;
        cp->d3 = new fclaw_clawpatch_3d_t;
        /* This patch will only be defined if we are on a manifold. */
        cp->d3->mp = fclaw3d_metric_patch_new();
    }

    return (void*) cp;
}

static
void* clawpatch_new_2d()
{
    return clawpatch_new(2);
}

static
void* clawpatch_new_3d()
{
    return clawpatch_new(3);
}

static
void clawpatch_delete(void *patchcp)
{
    FCLAW_ASSERT(patchcp != NULL);
    fclaw_clawpatch_t* cp = (fclaw_clawpatch_t*) patchcp;

    if(cp->patch_dim == 2)
    {
        fclaw2d_clawpatch_time_sync_delete(&cp->d2->registers);

        FCLAW_ASSERT(cp->d2->mp != NULL);
        fclaw2d_metric_patch_delete(&cp->d2->mp);

        delete cp->d2;
    }
    else 
    {
        FCLAW_ASSERT(cp->d3->mp != NULL);
        fclaw3d_metric_patch_delete(&cp->d3->mp);

        delete cp->d3;
    }

    delete cp;
    patchcp = NULL;
}

/* Maybe this should just be a 'build' function? */
static
void clawpatch_define(fclaw_global_t* glob,
                      fclaw_patch_t *patch,
                      int blockno, int patchno,
                      fclaw_build_mode_t build_mode)
{
    /* We are getting closer to getting rid the class fclaw_clawpatch_t */
    fclaw_clawpatch_t *cp = get_clawpatch(patch);

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    const fclaw_clawpatch_options_t *clawpatch_opt = 
                         fclaw_clawpatch_get_options(glob);

    cp->mx = clawpatch_opt->mx;
    cp->my = clawpatch_opt->my;
    cp->mz = clawpatch_opt->mz;

    cp->mbc = clawpatch_opt->mbc;
    cp->blockno = blockno;
    cp->meqn = clawpatch_opt->meqn;
    cp->maux = clawpatch_opt->maux;
    cp->mfields = clawpatch_opt->rhs_fields;

    for (int icorner=0; icorner < 4; icorner++)
    {
        fclaw_patch_set_block_corner_count(glob,patch,icorner,0);
    }

    fclaw_map_context_t* cont = fclaw_map_get(glob);

    int is_brick;
    if(glob->domain->refine_dim == 2)
    {
        is_brick = FCLAW2D_MAP_IS_BRICK(&cont);
    }else {
        //TODO mappings
        is_brick = (fclaw_opt->mi > 1) || (fclaw_opt->mj > 1) || (fclaw_opt->mk > 1);
    }

    cp->manifold = fclaw_opt->manifold;
    if (cp->manifold)
    {
        cp->xlower = patch->xlower;
        cp->ylower = patch->ylower;
        cp->zlower = patch->zlower;
        cp->xupper = patch->xupper;
        cp->yupper = patch->yupper;
        cp->zupper = patch->zupper;
    }
    else
    {
        double ax = fclaw_opt->ax;
        double bx = fclaw_opt->bx;
        double ay = fclaw_opt->ay;
        double by = fclaw_opt->by;
        double az = fclaw_opt->az;
        double bz = fclaw_opt->bz;

        double xl = patch->xlower;
        double yl = patch->ylower;
        double zl = patch->zlower;
        double xu = patch->xupper;
        double yu = patch->yupper;
        double zu = patch->zupper;

        double xlower, ylower, zlower, xupper, yupper, zupper;

        if (is_brick)
        {
            if(glob->domain->refine_dim == 2)
            {
                /* Scale to [0,1]x[0,1], based on blockno */
                double z;
                fclaw_map_2d_c2m_nomap_brick(cont,cp->blockno,xl,yl,&xlower,&ylower,&z);
                fclaw_map_2d_c2m_nomap_brick(cont,cp->blockno,xu,yu,&xupper,&yupper,&z);
                zlower = 0;
                zupper = 1;
            }
            else 
            {
                //TODO fclaw3d_map_c2m_nomap_brick
                //hardcoded for now
                fclaw_block_t *block = &glob->domain->blocks[blockno];

                //map in [0,1] for entire brick
                double block_xlower = block->vertices[0];
                double block_ylower = block->vertices[1];
                double block_zlower = block->vertices[2];

                xlower = (block_xlower + xl)/(double)fclaw_opt->mi;
                xupper = (block_xlower + xu)/(double)fclaw_opt->mi;
                ylower = (block_ylower + yl)/(double)fclaw_opt->mj;
                yupper = (block_ylower + yu)/(double)fclaw_opt->mj;
                zlower = (block_zlower + zl)/(double)fclaw_opt->mk;
                zupper = (block_zlower + zu)/(double)fclaw_opt->mk;
            }
        }
        else
        {
            xlower = xl;
            ylower = yl;
            zlower = zl;
            xupper = xu;
            yupper = yu;
            zupper = zu;
        }

        cp->xlower = ax + (bx - ax)*xlower;
        cp->xupper = ax + (bx - ax)*xupper;
        cp->ylower = ay + (by - ay)*ylower;
        cp->yupper = ay + (by - ay)*yupper;
        cp->zlower = az + (bz - az)*zlower;
        cp->zupper = az + (bz - az)*zupper;
    }

    cp->dx = (cp->xupper - cp->xlower)/cp->mx;
    cp->dy = (cp->yupper - cp->ylower)/cp->my;
    cp->dz = (cp->zupper - cp->zlower)/cp->mz;


    // bounds for 3d, even if 2d, z bounds will be ignored by Box constructor in 2d
    int ll[3];
    int ur[3];
    for (int idir = 0; idir < 3; idir++)
    {
        ll[idir] = 1-cp->mbc;
    }
    ur[0] = cp->mx + cp->mbc;
    ur[1] = cp->my + cp->mbc;
    ur[2] = cp->mz + cp->mbc;

    Box box(ll,ur,cp->patch_dim);    

    // This will destroy any existing memory n griddata.
    cp->griddata.define(box, cp->meqn);
    if (fclaw_opt->subcycle)
        cp->griddata_time_interpolated.define(box, cp->meqn);

    if (fclaw_opt->compute_error)
    {
        cp->griderror.define(box,cp->meqn);
        cp->exactsolution.define(box,cp->meqn);
    }

    if (clawpatch_opt->maux > 0)
    {        
        cp->aux.define(box,cp->maux);
        if (clawpatch_opt->save_aux)
            cp->aux_save.define(box,cp->maux);
    }

    if (clawpatch_opt->rhs_fields > 0)
    {
        cp->rhs.define(box,cp->mfields);
        if (fclaw_opt->compute_error)
        {
            cp->elliptic_error.define(box,cp->mfields);
            cp->elliptic_soln.define(box,cp->mfields);
        }
    }

    if (fclaw_opt->manifold)
    {
        /*  
            Allocate space for all metric terms. 

            Note: We pass in detailed info so that the metric patch 
            doesn't have to know about a clawpatch 
        */

        if(cp->patch_dim == 2)
        {
            fclaw2d_metric_patch_define(glob,patch, 
                                        cp->mx, cp->my, cp->mbc, 
                                        cp->dx, cp->dy,
                                        cp->xlower,cp->ylower, cp->xupper, cp->yupper, 
                                        blockno, patchno, build_mode);
        }
        else
        {
            fclaw3d_metric_patch_define(glob,patch, 
                                        cp->mx, cp->my, cp->mz, cp->mbc, 
                                        cp->dx, cp->dy, cp->dz,
                                        cp->xlower,cp->ylower, cp->zlower, 
                                        cp->xupper, cp->yupper, cp->zupper,
                                        blockno, patchno, build_mode);
        }
    }
    
    /* Build interface registers needed for conservation */
    if(cp->patch_dim == 2)
    {
        fclaw2d_clawpatch_time_sync_new(glob,patch,
                                      blockno,patchno,&cp->d2->registers);
    }

    if (build_mode != FCLAW_BUILD_FOR_UPDATE)
        /* If we are building ghost patches, we don't need all the patch memory */
        return;

    cp->griddata_last.define(box, cp->meqn);
    cp->griddata_save.define(box, cp->meqn);

}

static
void clawpatch_build(fclaw_global_t *glob,
                     fclaw_patch_t *patch,
                     int blockno,
                     int patchno,
                     void *user)
{
    fclaw_build_mode_t build_mode =  *((fclaw_build_mode_t*) user);
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    const fclaw_clawpatch_options_t* clawpatch_opt = 
                         fclaw_clawpatch_get_options(glob);

    /* Allocate space for all required arrays */
    clawpatch_define(glob,patch,blockno,patchno,build_mode);

    if (fclaw_opt->manifold)
    { 
        /* Computes all metric terms : (xp,yp,zp,xd,yd,zd), 
        areas/volumes, basis vectors and (in 2d) curvatures
        and surface normals. 

        Note : Area/volume computations are expensive, since 
        they are all computed on finest levels 

        Note : When building from fine, call routine below 
        where areas/volumes are averaged from fine patches */
        if(clawpatch_opt->patch_dim == 2)
        {
            fclaw2d_metric_patch_build(glob,patch,blockno,patchno);
        }
        else 
        {
            fclaw3d_metric_patch_build(glob,patch,blockno,patchno);
        }
    }

    if(clawpatch_opt->patch_dim == 2)
    {
        /* Setup for conservation correction */
        fclaw2d_clawpatch_time_sync_setup(glob,patch,blockno,patchno);
    }
}

static
void clawpatch_build_from_fine(fclaw_global_t *glob,
                               fclaw_patch_t *fine_patches,
                               fclaw_patch_t *coarse_patch,
                               int blockno,
                               int coarse_patchno,
                               int fine0_patchno,
                               fclaw_build_mode_t build_mode)
{
    const fclaw_clawpatch_options_t* clawpatch_opt = 
                         fclaw_clawpatch_get_options(glob);
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    clawpatch_define(glob,coarse_patch,blockno,coarse_patchno,build_mode);

    if (fclaw_opt->manifold)
    {
        /* Average areas/volume from fine grid to coarse grid */
        if(clawpatch_opt->patch_dim == 2)
        {
            fclaw2d_metric_patch_build_from_fine(glob,fine_patches,coarse_patch,
                                                 blockno,coarse_patchno,
                                                 fine0_patchno);
        }
        else 
        {
            fclaw3d_metric_patch_build_from_fine(glob,fine_patches,coarse_patch,
                                                 blockno,coarse_patchno,
                                                 fine0_patchno);
        }
    }

    if(clawpatch_opt->patch_dim == 2)
    {
        /* Setup for conservation correction */
        fclaw2d_clawpatch_time_sync_setup(glob,coarse_patch,blockno,coarse_patchno);
    }
}


/* -------------------------------- time stepping ------------------------------------- */

static
void clawpatch_save_step(fclaw_global_t* glob,
                         fclaw_patch_t* patch)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    cp->griddata_save = cp->griddata;

    /* Some aux arrays are time dependent, or contain part of the solution.  In this case, 
       we should save the aux array in case we need to re-take a time step */
    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
    if (clawpatch_opt->save_aux)
        cp->aux_save = cp->aux;
}


static
void clawpatch_restore_step(fclaw_global_t* glob,
                            fclaw_patch_t* patch)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    cp->griddata = cp->griddata_save;

    /* Restore the aux array after before retaking a time step */
    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
    if (clawpatch_opt->save_aux)
        cp->aux = cp->aux_save;
}

static
void clawpatch_setup_timeinterp(fclaw_global_t *glob,
                                fclaw_patch_t *patch,
                                double alpha)
{
    /* We use the pack size here to make sure we are setting
       everything correctly;  it isn't needed for memory
       allocation */
    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);

    int mx = clawpatch_opt->mx;
    int my = clawpatch_opt->my;
    int mz = clawpatch_opt->mz;

    int meqn = clawpatch_opt->meqn;
    int mbc = clawpatch_opt->mbc;
    int mint = clawpatch_opt->interp_stencil_width/2+1;  

    int hole = (mx - 2*mint)*(my - 2*mint);  /* Hole in center */
    int wg = mx*my;  /* Whole grid but no ghost cells.  
                        Ghost cells will be averaged from finer level. */
    if(clawpatch_opt->patch_dim == 3)
    {
        hole *= mz;
        wg *= mz;
    }

    FCLAW_ASSERT(hole >= 0);

    int psize = (wg - hole)*meqn;
    FCLAW_ASSERT(psize > 0);

    /* Store time interpolated data that will be use in coarse grid
       exchanges */
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    double *qlast = cp->griddata_last.dataPtr();
    double *qcurr = cp->griddata.dataPtr();
    double *qinterp = cp->griddata_time_interpolated.dataPtr();

    int ierror;

    /* Do interpolation only on interior, since ghost cells in qcurr
       are invalid and will lead to floating point exceptions.
       We do a ghost cell update at the intermediate time level.  The
       neighboring fine grid will average to ghost cells of the interpolated
       level, then the interpolated level is used to interpolate to fine grid
       ghost cells. */
    fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);
    if(clawpatch_opt->patch_dim == 2)
    {
        clawpatch_vt->d2->fort_timeinterp(&mx,&my,&mbc,&meqn,&psize,
                                          qcurr,qlast,qinterp,&alpha,&ierror);
    }
    else
    {
        clawpatch_vt->d3->fort_timeinterp(&mx,&my,&mz, &mbc,&meqn,&psize,
                                          qcurr,qlast,qinterp,&alpha,&ierror);
    }
}


/* ------------------------------------- Ghost filling  ------------------------------- */

static
void clawpatch_copy_face(fclaw_global_t *glob,
                         fclaw_patch_t *patch,
                         fclaw_patch_t *neighbor_patch,
                         int iface,
                         int time_interp,
                         fclaw_patch_transform_data_t *transform_data)

{

    const fclaw_clawpatch_options_t *clawpatch_opt = 
                fclaw_clawpatch_get_options(glob);

    int mbc = clawpatch_opt->mbc;

    /* This routine might be called between two time-sync patches */
    int meqn;
    double *qthis;
    fclaw_clawpatch_timesync_data(glob,patch,time_interp,&qthis,&meqn);

    double *qneighbor;
    fclaw_clawpatch_timesync_data(glob,neighbor_patch,time_interp,&qneighbor,&meqn);

    if (fill_ghost(glob,time_interp))
    {
        fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);
        if(clawpatch_opt->patch_dim == 2)
        {
            int mx = clawpatch_opt->mx;
            int my = clawpatch_opt->my;
            clawpatch_vt->d2->fort_copy_face(&mx,&my,&mbc,&meqn,qthis,
                                             qneighbor,&iface,&transform_data);

        }
        else 
        {
            int mx = clawpatch_opt->mx;
            int my = clawpatch_opt->my;
            int mz = clawpatch_opt->mz;
            clawpatch_vt->d3->fort_copy_face(&mx,&my,&mz,&mbc,&meqn,qthis,
                                             qneighbor,&iface,&transform_data);
        }
    }
}

static
void clawpatch_average_face(fclaw_global_t *glob,
                            fclaw_patch_t *coarse_patch,
                            fclaw_patch_t *fine_patch,
                            int idir,
                            int iface_coarse,
                            int p4est_refineFactor,
                            int refratio,
                            int time_interp,
                            int igrid,
                            fclaw_patch_transform_data_t* transform_data)
{
    int meqn;
    double *qcoarse;
    fclaw_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
    double *qfine = fclaw_clawpatch_get_q(glob,fine_patch);


    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
    int mbc = clawpatch_opt->mbc;

    const fclaw_options_t* fclaw_opt = fclaw_get_options(glob);
    int manifold = fclaw_opt->manifold;
    fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);


    if(clawpatch_opt->patch_dim == 2)
    {
        double *areacoarse = clawpatch_get_area(glob, coarse_patch);
        double *areafine = clawpatch_get_area(glob, fine_patch);

        int mx = clawpatch_opt->mx;
        int my = clawpatch_opt->my;
        /* These will be empty for non-manifolds cases */
        clawpatch_vt->d2->fort_average_face(&mx,&my,&mbc,&meqn,qcoarse,qfine,
                                            areacoarse,areafine, &idir,&iface_coarse, 
                                            &p4est_refineFactor, &refratio,
                                            &igrid,&manifold,&transform_data);
    }
    else
    {
        double *volcoarse = clawpatch_get_volume(glob, coarse_patch);
        double *volfine = clawpatch_get_volume(glob, fine_patch);
        int mx = clawpatch_opt->mx;
        int my = clawpatch_opt->my;
        int mz = clawpatch_opt->mz;
        clawpatch_vt->d3->fort_average_face(&mx,&my,&mz,&mbc,&meqn,qcoarse,qfine,
                                            volcoarse,volfine, &idir,&iface_coarse, 
                                            &p4est_refineFactor, &refratio,
                                            &igrid,&manifold,&transform_data);
    }
}

static
void clawpatch_interpolate_face(fclaw_global_t *glob,
                                fclaw_patch_t *coarse_patch,
                                fclaw_patch_t *fine_patch,
                                int idir,
                                int iface_coarse,
                                int p4est_refineFactor,
                                int refratio,
                                int time_interp,
                                int igrid,
                                fclaw_patch_transform_data_t* transform_data)
{

    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);

    int meqn;
    double *qcoarse;
    fclaw_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
    double *qfine = fclaw_clawpatch_get_q(glob,fine_patch);

    int mbc = clawpatch_opt->mbc;

    if (fill_ghost(glob,time_interp))
    {
        fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);
        if(clawpatch_opt->patch_dim == 2)
        {    
            int mx = clawpatch_opt->mx;
            int my = clawpatch_opt->my;
            clawpatch_vt->d2->fort_interpolate_face(&mx,&my,&mbc,&meqn,qcoarse,qfine,&idir,
                                                    &iface_coarse, &p4est_refineFactor,
                                                    &refratio,&igrid,&transform_data);
        }
        else
        {
            int mx = clawpatch_opt->mx;
            int my = clawpatch_opt->my;
            int mz = clawpatch_opt->mz;
            clawpatch_vt->d3->fort_interpolate_face(&mx,&my,&mz, &mbc,&meqn,qcoarse,qfine,
                                                    &idir,&iface_coarse, &p4est_refineFactor,
                                                    &refratio,&igrid,&transform_data);
        }
    }
}

static
void clawpatch_copy_edge(fclaw_global_t *glob,
                         fclaw_patch_t *patch,
                         fclaw_patch_t *edge_patch,
                         int coarse_blockno,
                         int fine_blockno,
                         int iedge,
                         int time_interp,
                         fclaw_patch_transform_data_t *transform_data)
{

    const fclaw_clawpatch_options_t *clawpatch_opt = 
                 fclaw_clawpatch_get_options(glob);
    int mbc = clawpatch_opt->mbc;

    int meqn;
    double *qthis;
    fclaw_clawpatch_timesync_data(glob,patch,time_interp,&qthis,&meqn);

    double *qedge;
    fclaw_clawpatch_timesync_data(glob,edge_patch,time_interp,&qedge,&meqn);

    if (fill_ghost(glob,time_interp))
    {
        fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);
        FCLAW_ASSERT(clawpatch_opt->patch_dim == 3);
        int mx = clawpatch_opt->mx;
        int my = clawpatch_opt->my;
        int mz = clawpatch_opt->mz;
        clawpatch_vt->d3->fort_copy_edge(&mx,&my,&mz, &mbc,&meqn,qthis,qedge,
                                         &iedge, &transform_data);
    }
}

static
void clawpatch_average_edge(fclaw_global_t *glob,
                            fclaw_patch_t *coarse_patch,
                            fclaw_patch_t *fine_patch,
                            int coarse_edge,
                            int time_interp,
                            fclaw_patch_transform_data_t* transform_data)
{
    int meqn;
    double *qcoarse;
    fclaw_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);

    double *qfine = fclaw_clawpatch_get_q(glob,fine_patch);

    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
    int mbc = clawpatch_opt->mbc;

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    int manifold = fclaw_opt->manifold;
    if (fill_ghost(glob,time_interp))
    {
        int refratio = 2;
        fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);
        FCLAW_ASSERT(clawpatch_opt->patch_dim == 3);
        /* These will be empty for non-manifolds cases */
        double *volcoarse = clawpatch_get_volume(glob, coarse_patch);
        double *volfine = clawpatch_get_volume(glob, fine_patch);

        int mx = clawpatch_opt->mx;
        int my = clawpatch_opt->my;
        int mz = clawpatch_opt->mz;        
        clawpatch_vt->d3->fort_average_edge(&mx,&my,&mz,&mbc,&meqn,
                                            &refratio,qcoarse,qfine,
                                            volcoarse,volfine,
                                            &manifold,&coarse_edge,&transform_data);
    }
}

static
void clawpatch_interpolate_edge(fclaw_global_t *glob,
                                fclaw_patch_t *coarse_patch,
                                fclaw_patch_t *fine_patch,
                                int coarse_edge,
                                int time_interp,
                                fclaw_patch_transform_data_t* transform_data)
{
    int meqn;
    double *qcoarse;
    fclaw_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);

    double *qfine = fclaw_clawpatch_get_q(glob,fine_patch);

    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
    int mbc = clawpatch_opt->mbc;

    if (fill_ghost(glob,time_interp))
    {
        int refratio = 2;
        fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);
        FCLAW_ASSERT(clawpatch_opt->patch_dim == 3);
        /* These will be empty for non-manifolds cases */

        int mx = clawpatch_opt->mx;
        int my = clawpatch_opt->my;
        int mz = clawpatch_opt->mz;        
        clawpatch_vt->d3->fort_interpolate_edge(&mx,&my,&mz,&mbc,&meqn,
                                                &refratio,qcoarse,qfine,
                                                &coarse_edge,&transform_data);    
    }
}

static
void clawpatch_copy_corner(fclaw_global_t *glob,
                           fclaw_patch_t *patch,
                           fclaw_patch_t *corner_patch,
                           int coarse_blockno,
                           int fine_blockno,
                           int icorner,
                           int time_interp,
                           fclaw_patch_transform_data_t *transform_data)
{

    const fclaw_clawpatch_options_t *clawpatch_opt = 
                 fclaw_clawpatch_get_options(glob);
    int mbc = clawpatch_opt->mbc;

    int meqn;
    double *qthis;
    fclaw_clawpatch_timesync_data(glob,patch,time_interp,&qthis,&meqn);

    double *qcorner;
    fclaw_clawpatch_timesync_data(glob,corner_patch,time_interp,&qcorner,&meqn);

    if (fill_ghost(glob,time_interp))
    {
        fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);
        if(clawpatch_opt->patch_dim == 2)
        {
            int mx = clawpatch_opt->mx;
            int my = clawpatch_opt->my;
            clawpatch_vt->d2->fort_copy_corner(&mx,&my,&mbc,&meqn,qthis,qcorner,
                                               &icorner, &transform_data);
        }
        else
        {
            int mx = clawpatch_opt->mx;
            int my = clawpatch_opt->my;
            int mz = clawpatch_opt->mz;
            clawpatch_vt->d3->fort_copy_corner(&mx,&my,&mz, &mbc,&meqn,qthis,qcorner,
                                               &icorner, &transform_data);
        }
    }
}


static
void clawpatch_average_corner(fclaw_global_t *glob,
                              fclaw_patch_t *coarse_patch,
                              fclaw_patch_t *fine_patch,
                              int coarse_blockno,
                              int fine_blockno,
                              int coarse_corner,
                              int time_interp,
                              fclaw_patch_transform_data_t* transform_data)
{
    int meqn;
    double *qcoarse;
    fclaw_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);

    double *qfine = fclaw_clawpatch_get_q(glob,fine_patch);

    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
    int mbc = clawpatch_opt->mbc;

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    int manifold = fclaw_opt->manifold;
    if (fill_ghost(glob,time_interp))
    {
        int refratio = 2;
        fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);

        if(clawpatch_opt->patch_dim == 2)
        {
            double *areacoarse = clawpatch_get_area(glob, coarse_patch);
            double *areafine = clawpatch_get_area(glob, fine_patch);

            int mx = clawpatch_opt->mx;
            int my = clawpatch_opt->my;
               /* These will be empty for non-manifolds cases */
            clawpatch_vt->d2->fort_average_corner(&mx,&my,&mbc,&meqn,
                                                  &refratio,qcoarse,qfine,
                                                  areacoarse,areafine,
                                                  &manifold,&coarse_corner,&transform_data);
        }
        else
        {
             /* These will be empty for non-manifolds cases */
            double *volcoarse = clawpatch_get_volume(glob, coarse_patch);
            double *volfine = clawpatch_get_volume(glob, fine_patch);

            int mx = clawpatch_opt->mx;
            int my = clawpatch_opt->my;
            int mz = clawpatch_opt->mz;        
            clawpatch_vt->d3->fort_average_corner(&mx,&my,&mz,&mbc,&meqn,
                                                  &refratio,qcoarse,qfine,
                                                  volcoarse,volfine,
                                                  &manifold,&coarse_corner,&transform_data);
        }
    }
}

static
void clawpatch_interpolate_corner(fclaw_global_t* glob,
                                  fclaw_patch_t* coarse_patch,
                                  fclaw_patch_t* fine_patch,
                                  int coarse_blockno,
                                  int fine_blockno,
                                  int coarse_corner,
                                  int time_interp,
                                  fclaw_patch_transform_data_t* transform_data)

{
    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
    int mbc = clawpatch_opt->mbc;

    int meqn;
    double *qcoarse;
    fclaw_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);

    double *qfine = fclaw_clawpatch_get_q(glob,fine_patch);

    if (fill_ghost(glob,time_interp))
    {
        int refratio = 2;
        fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);
        if(clawpatch_opt->patch_dim == 2)
        {
            int mx = clawpatch_opt->mx;
            int my = clawpatch_opt->my;
            clawpatch_vt->d2->fort_interpolate_corner(&mx,&my,&mbc,&meqn,
                                                      &refratio,qcoarse,qfine,
                                                      &coarse_corner,&transform_data);    
        }
        else
        {
            int mx = clawpatch_opt->mx;
            int my = clawpatch_opt->my;
            int mz = clawpatch_opt->mz;
            clawpatch_vt->d3->fort_interpolate_corner(&mx,&my,&mz,&mbc,&meqn,
                                                      &refratio,qcoarse,qfine,
                                                      &coarse_corner,&transform_data);    
        }
    }

}



/* -------------------------------- Regridding functions ------------------------------ */

static
int clawpatch_tag4refinement(fclaw_global_t *glob,
                             fclaw_patch_t *patch,
                             int blockno, int patchno,
                             int initflag)
{
    int meqn;
    double *q;
    fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    double refine_threshold = fclaw_opt->refine_threshold;

    int tag_patch;
    if (refine_threshold < 0) 
    {
        /* Always refine */
        tag_patch = 1;
    }
    else
    {
        tag_patch = 0;    

        /* This allows the user to specify a "exceeds_th" */
        fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);

        fclaw_global_set_static(glob);
        if(clawpatch_vt->patch_dim == 2)
        {

            int mx,my,mbc;
            double xlower,ylower,dx,dy;
            fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                        &xlower,&ylower,&dx,&dy);
            clawpatch_vt->d2->fort_tag4refinement(&mx,&my,&mbc,&meqn,&xlower,&ylower,
                                                  &dx,&dy, &blockno, q,
                                                  &refine_threshold,
                                                  &initflag,&tag_patch);
        }
        else
        {
            int mx,my,mz,mbc;
            double xlower,ylower,zlower,dx,dy,dz;
            fclaw_clawpatch_3d_grid_data(glob,patch,&mx,&my,&mz, &mbc,
                                        &xlower,&ylower,&zlower, &dx,&dy,&dz);

            clawpatch_vt->d3->fort_tag4refinement(&mx,&my,&mz, &mbc,&meqn,
                                                  &xlower,&ylower,&zlower,
                                                  &dx,&dy, &dz, &blockno, q,
                                                  &refine_threshold,
                                                  &initflag,&tag_patch);
        }
        fclaw_global_clear_static();
    }
    return tag_patch;
}

static
int clawpatch_tag4coarsening(fclaw_global_t *glob,
                             fclaw_patch_t *fine_patches,
                             int blockno,
                             int patchno,
                             int initflag)
{
    const int num_children = fclaw_domain_num_children(glob->domain);
    fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);
    int mx,my,mz,mbc,meqn;
    double xlower[num_children],ylower[num_children],zlower[num_children];
    double dx,dy,dz;

    double *q[num_children];
    for (int igrid = 0; igrid < num_children; igrid++)
    {
        fclaw_clawpatch_soln_data(glob,&fine_patches[igrid],&q[igrid],&meqn);
        if(clawpatch_vt->patch_dim == 2)
        {
            fclaw_clawpatch_2d_grid_data(glob,&fine_patches[igrid],&mx,&my,&mbc,
                                        &xlower[igrid],&ylower[igrid],&dx,&dy);
        }
        else 
        {
            /* For extruded meshes, zlower is the same for all patches. */
            fclaw_clawpatch_3d_grid_data(glob,&fine_patches[igrid],&mx,&my,&mz,&mbc,
                                        &xlower[igrid],&ylower[igrid],&zlower[igrid],&dx,&dy,&dz);
        }
    }

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    double coarsen_threshold = fclaw_opt->coarsen_threshold;

    int tag_patch = 0;
    if (coarsen_threshold > 0) 
    {        
        fclaw_global_set_static(glob);
        if(clawpatch_vt->patch_dim == 2)
        {
            /* Get xlower,ylower for each grid. */
            clawpatch_vt->d2->fort_tag4coarsening(&mx,&my,&mbc,&meqn,
                                                  xlower,ylower,&dx,&dy,
                                                  &blockno, q[0],q[1],q[2],q[3],
                                                  &coarsen_threshold,&initflag,&tag_patch);
        }
        else if(glob->domain->refine_dim == 2)
        {
            clawpatch_vt->d3->fort_tag4coarsening_3dx(&mx,&my,&mz,&mbc,&meqn,
                                                  xlower,ylower,zlower,&dx,&dy,&dz,
                                                  &blockno, q[0],q[1],q[2],q[3],
                                                  &coarsen_threshold,&initflag,&tag_patch);
        }
        else
        {
            clawpatch_vt->d3->fort_tag4coarsening(&mx,&my,&mz,&mbc,&meqn,
                                                  xlower,ylower,zlower,&dx,&dy,&dz,
                                                  &blockno, q[0],q[1],q[2],q[3],
                                                  q[4],q[5],q[6],q[7],
                                                  &coarsen_threshold,&initflag,&tag_patch);
        }
        fclaw_global_clear_static();
    }
    else
    {
        /* Never coarsen */
    }
    return tag_patch == 1;
}

static
void clawpatch_interpolate2fine(fclaw_global_t* glob,
                                fclaw_patch_t *coarse_patch,
                                fclaw_patch_t* fine_patches,
                                int this_blockno, int coarse_patchno,
                                int fine0_patchno)
{

    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
    int mbc = clawpatch_opt->mbc;
    int meqn = clawpatch_opt->meqn;

    /* Loop over four siblings (z-ordering) */
    for (int igrid = 0; igrid < fclaw_domain_num_children(glob->domain); igrid++)
    {
        fclaw_patch_t *fine_patch = &fine_patches[igrid];
        double *qfine = fclaw_clawpatch_get_q(glob,fine_patch);

        const fclaw_options_t* fclaw_opt = fclaw_get_options(glob);

        fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);


        double *qcoarse = fclaw_clawpatch_get_q(glob,coarse_patch);

        if(clawpatch_opt->patch_dim == 2)
        {
            double *areacoarse = clawpatch_get_area(glob, coarse_patch);
            double *areafine = NULL;
            if (fclaw_opt->manifold)
                areafine = clawpatch_get_area(glob, fine_patch);

            int mx = clawpatch_opt->mx;
            int my = clawpatch_opt->my;
            clawpatch_vt->d2->fort_interpolate2fine(&mx,&my,&mbc,&meqn,qcoarse,qfine,
                                                    areacoarse, areafine, &igrid,
                                                    &fclaw_opt->manifold);
        }
        else 
        {
            double *volcoarse = clawpatch_get_volume(glob, coarse_patch);
            double *volfine = NULL;
            if (fclaw_opt->manifold)
                 volfine = clawpatch_get_volume(glob, fine_patch);

            int mx = clawpatch_opt->mx;
            int my = clawpatch_opt->my;
            int mz = clawpatch_opt->mz;
            clawpatch_vt->d3->fort_interpolate2fine(&mx,&my,&mz, &mbc,&meqn,qcoarse,qfine,
                                                    volcoarse, volfine, &igrid,
                                                    &fclaw_opt->manifold);
        }
    }
}

static
void clawpatch_average2coarse(fclaw_global_t *glob,
                              fclaw_patch_t *fine_patches,
                              fclaw_patch_t *coarse_patch,
                              int blockno, int fine0_patchno,
                              int coarse_patchno)

{
    const fclaw_clawpatch_options_t *clawpatch_opt = 
                    fclaw_clawpatch_get_options(glob);    
    int mbc = clawpatch_opt->mbc;
    int meqn = clawpatch_opt->meqn;


    for(int igrid = 0; igrid < fclaw_domain_num_children(glob->domain); igrid++)
    {
        fclaw_patch_t *fine_patch = &fine_patches[igrid];
        double *qfine = fclaw_clawpatch_get_q(glob,fine_patch);
        double *qcoarse = fclaw_clawpatch_get_q(glob,coarse_patch);

        const fclaw_options_t* fclaw_opt = fclaw_get_options(glob);

        fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);


        if(clawpatch_opt->patch_dim == 2)
        {
            double *areacoarse = clawpatch_get_area(glob, coarse_patch);

            double *areafine =  NULL;
            if (fclaw_opt->manifold)
                areafine = clawpatch_get_area(glob, fine_patch);

            int mx = clawpatch_opt->mx;
            int my = clawpatch_opt->my;
            clawpatch_vt->d2->fort_average2coarse(&mx,&my,&mbc,&meqn,qcoarse,qfine,
                                                  areacoarse, areafine, &igrid,
                                                  &fclaw_opt->manifold);
        }
        else
        {
            double *volcoarse = clawpatch_get_volume(glob, coarse_patch);

            double *volfine = NULL;
            if (fclaw_opt->manifold)
                volfine = clawpatch_get_volume(glob, fine_patch);

            int mx = clawpatch_opt->mx;
            int my = clawpatch_opt->my;
            int mz = clawpatch_opt->mz;
            clawpatch_vt->d3->fort_average2coarse(&mx,&my,&mz, &mbc,&meqn,qcoarse,qfine,
                                                  volcoarse, volfine, &igrid,
                                                  &fclaw_opt->manifold);
        }

    }
}

/* ------------------------------ Parallel ghost patches ------------------------------ */

/* This is called just to get a count of how much to pack */
static
size_t clawpatch_ghost_pack_elems(fclaw_global_t* glob)
{
    const fclaw_clawpatch_options_t *clawpatch_opt = 
                             fclaw_clawpatch_get_options(glob);

    int mx,my,mz;
    if(clawpatch_opt->patch_dim == 2)
    {
        mx = clawpatch_opt->mx;
        my = clawpatch_opt->my;
        mz = 1;
    }
    else
    {
        mx = clawpatch_opt->mx;
        my = clawpatch_opt->my;
        mz = clawpatch_opt->mz;
    }
    int mbc = clawpatch_opt->mbc;
    int meqn = clawpatch_opt->meqn;

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    int refratio = fclaw_opt->refratio;
    int packextra = fclaw_opt->ghost_patch_pack_numextrafields;
    int packarea = fclaw_opt->ghost_patch_pack_area && fclaw_opt->manifold;
    int packregisters = fclaw_opt->time_sync;

    int mint = refratio*mbc;
    int nghost = mbc;

    /* Include size of conservation registers.  Save fluxes on each size, 
       even though only one or two sides may be used. */
    // int frsize = 12*meqn*(mx + my); 

    int frsize = packregisters ? 2*(4*meqn+2)*(mx + my) : 0;

    if (clawpatch_opt->patch_dim == 3 && packregisters)
    {        
        fclaw_global_essentialf("clawpatch_ghost_comm: Conservation fix not yet " \
                                "implemented in 3d\n");
        exit(0);
    }

    int wg = (2*nghost + mx)*(2*nghost + my);  /* Whole grid     */
    int hole = (mx - 2*mint)*(my - 2*mint);    /* Hole in center */

    if(glob->domain->refine_dim == 2 && clawpatch_opt->patch_dim == 3)
    {
        wg *= (mz + 2*nghost);
        hole *= (mz + 2*nghost);  
    }
    else if (clawpatch_opt->patch_dim == 3)
    {
        wg *= (mz + 2*nghost);
        hole *= (mz - 2*mint);  
    }

    FCLAW_ASSERT(hole >= 0);

    size_t psize = (wg - hole)*(meqn + packarea + packextra) + frsize;
    FCLAW_ASSERT(psize >= 0);

    return psize;
}    


static
void clawpatch_ghost_comm(fclaw_global_t* glob,
                          fclaw_patch_t* patch,
                          void *unpack_from_here, int time_interp,
                          int packmode)
{
    const fclaw_clawpatch_options_t *clawpatch_opt = 
                            fclaw_clawpatch_get_options(glob);

    int mx = clawpatch_opt->mx;
    int my = clawpatch_opt->my;
    int mz = clawpatch_opt->mz;

    int mbc = clawpatch_opt->mbc;
    int meqn = clawpatch_opt->meqn;

    int packarea = packmode/2;   // (0,1)/2 = 0;  (2,3)/2 = 1;

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    int packextra = fclaw_opt->ghost_patch_pack_numextrafields;
    int packregisters = fclaw_opt->time_sync;  // For conservation
    int refratio = fclaw_opt->refratio;

    int mint = refratio*mbc;   /* # interior cells needed for averaging */
    int nghost = mbc;          /* # ghost values needed for interpolation */

    /* Include size of conservation registers.  Save fluxes on each size, 
       even though only one or two sides may be used.

       y-face : gp,gm,f(qfront),f(qback) + edge-length + area : 4*meqn + 2
       x-face : fp,fm,f(ql),f(qr) + edge-length + area        : 4*meqn + 2
    
       Multiply everything by mz to get fluxes on each face of cube. 
       (Time sync not yet implemented in 3d, though). 

       */
    int frsize = packregisters ? 2*(4*meqn+2)*(mx + my) : 0;

    if (clawpatch_opt->patch_dim == 3 && packregisters)
    {        
        fclaw_global_essentialf("clawpatch_ghost_comm: Conservation fix not yet " \
                                "implemented in 3d\n");
        exit(0);
    }

    /* wg   : whole grid
       hole : part we don't store */ 
    /* This is computed twice - here, and in fclaw2d_clawpatch_ghost_packsize */
    int wg = (2*nghost + mx)*(2*nghost + my);
    int hole = (mx - 2*mint)*(my - 2*mint);  /* Hole in center */

    if(patch->refine_dim == 2 && clawpatch_opt->patch_dim == 3)
    {
        wg *= (mz + 2*nghost);
        hole *= (mz + 2*nghost);  
    }
    else if(clawpatch_opt->patch_dim == 3)
    {
        wg *= (mz + 2*nghost);
        hole *= (mz - 2*mint);  
    }

    FCLAW_ASSERT(hole >= 0);

    size_t psize = (wg - hole)*(meqn + packarea + packextra) + frsize;
    FCLAW_ASSERT(psize > 0);

    size_t psize_check = clawpatch_ghost_pack_elems(glob);

    /* Check with routine that is used to allocate space */
    FCLAW_ASSERT(psize == psize_check);

    int ierror;
    fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);

    double *qthis;
    fclaw_clawpatch_timesync_data(glob,patch,time_interp,&qthis,&meqn);
    double *qpack = (double*) unpack_from_here;
    int size;
    if(clawpatch_opt->patch_dim == 2)
    {
        size = (wg - hole)*(meqn + packarea);
        double *area = clawpatch_get_area(glob, patch);    
        clawpatch_vt->d2->fort_local_ghost_pack(&mx,&my,&mbc,&meqn,&mint,qthis,area,
                                                qpack,&size,&packmode,&ierror);
    }
    else
    {
        size = (wg - hole)*(meqn + packarea);
        double *vol = clawpatch_get_volume(glob, patch);    
        clawpatch_vt->d3->fort_local_ghost_pack(&mx,&my,&mz,&mbc,&meqn,&mint,qthis,vol,
                                                qpack,&size,&packmode,&ierror);
    }
    FCLAW_ASSERT(ierror == 0);

    /* Advance pointer */
    qpack += size;

    int extrasize = (wg - hole)*(packextra);
    if (packextra)
    {
        FCLAW_ASSERT(extrasize > 0);
        FCLAW_ASSERT(clawpatch_vt->local_ghost_pack_aux != NULL);
        clawpatch_vt->local_ghost_pack_aux(glob,patch,mint,
                                           qpack,extrasize,
                                           packmode,&ierror);
        FCLAW_ASSERT(ierror == 0);
        qpack += extrasize; /* Advance pointer */
    }

    FCLAW_ASSERT(ierror == 0);
    if (packregisters)
    {
        FCLAW_ASSERT(frsize > 0);
        FCLAW_ASSERT(clawpatch_vt->time_sync_pack_registers != NULL);
        fclaw_clawpatch_packmode_t frpackmode = packmode % 2 == 0 ?  
                                                    CLAWPATCH_REGISTER_PACK : 
                                                    CLAWPATCH_REGISTER_UNPACK;
        clawpatch_vt->time_sync_pack_registers(glob, patch,
                                               qpack,frsize,frpackmode,
                                               &ierror);
        FCLAW_ASSERT(ierror == 0);
    }


    if (ierror > 0)
    {
        fclaw_global_essentialf("clawpatch_ghost_comm  : ierror = %d\n",ierror);
        exit(0);
    }
}


static size_t clawpatch_ghost_packsize(fclaw_global_t* glob)
{
    size_t esize = clawpatch_ghost_pack_elems(glob);
    return esize*sizeof(double);
}

static
void clawpatch_local_ghost_pack(fclaw_global_t *glob,
                                fclaw_patch_t *patch,
                                void *patch_data,
                                int time_interp)
{
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    int packarea = fclaw_opt->ghost_patch_pack_area && fclaw_opt->manifold;
    int packmode = 2*packarea;  // 0 or 2  (for pack)

    clawpatch_ghost_comm(glob,patch,patch_data, time_interp,packmode);
}

static
void clawpatch_remote_ghost_unpack(fclaw_global_t* glob,
                                   fclaw_patch_t* patch,
                                   int blockno,
                                   int patchno,
                                   void *qdata, int time_interp)
{
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    int packarea = fclaw_opt->ghost_patch_pack_area && fclaw_opt->manifold;
    int packmode = 2*packarea + 1;  // 1 or 3  (for unpack)

    clawpatch_ghost_comm(glob,patch,qdata,time_interp,packmode);
}

static
void clawpatch_remote_ghost_build(fclaw_global_t *glob,
                                  fclaw_patch_t *patch,
                                  int blockno,
                                  int patchno,
                                  fclaw_build_mode_t build_mode)
{
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    const fclaw_clawpatch_options_t* clawpatch_opt = 
                             fclaw_clawpatch_get_options(glob);

    clawpatch_define(glob,patch,blockno,patchno,build_mode);

    if (fclaw_opt->manifold)
    {
        if (build_mode != FCLAW_BUILD_FOR_GHOST_AREA_PACKED)
        {
            /* Cell areas/volumes are not sent as MPI messages and so must
               be recomputed
            */
            //fclaw2d_metric_patch_build(glob,patch,blockno,patchno);
            if(clawpatch_opt->patch_dim == 2)
            {
                fclaw2d_metric_patch_compute_area(glob,patch,blockno,patchno);
            }
            else 
            {
                fclaw3d_metric_patch_compute_volume(glob,patch,blockno,patchno);
            }
        }
    }
    /* Any metric terms we might need for the registers are packed */
#if 0    
    /* Setup for conservation correction */
    fclaw2d_clawpatch_time_sync_setup(glob,this_patch,blockno,patchno);
#endif    
}

static
void clawpatch_remote_ghost_delete(void *patchcp)
{
    clawpatch_delete(patchcp);
}

/* ---------------------------- Parallel partitioning --------------------------------- */

static
size_t clawpatch_partition_packsize(fclaw_global_t* glob)
{
    const fclaw_clawpatch_options_t *clawpatch_opt 
                              = fclaw_clawpatch_get_options(glob);
    int mbc = clawpatch_opt->mbc;
    int meqn = clawpatch_opt->meqn;

    size_t psize;
    if(clawpatch_opt->patch_dim == 2)
    {
        int mx = clawpatch_opt->mx;
        int my = clawpatch_opt->my;
        psize = meqn*(2*mbc + mx)*(2*mbc + my);  /* Store area */
    }
    else
    {
        int mx = clawpatch_opt->mx;
        int my = clawpatch_opt->my;
        int mz = clawpatch_opt->mz;
        psize = meqn*(2*mbc + mx)*(2*mbc + my)*(2*mbc + mz);
    }

    return psize*sizeof(double);
}

static
void clawpatch_partition_pack(fclaw_global_t *glob,
                              fclaw_patch_t *patch,
                              int blockno,
                              int patchno,
                              void *pack_data_here)
    {
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    FCLAW_ASSERT(cp != NULL);

    cp->griddata.copyToMemory((double*) pack_data_here);
}

static
void clawpatch_partition_unpack(fclaw_global_t *glob,  
                                fclaw_domain_t *new_domain,
                                fclaw_patch_t *patch,
                                int blockno,
                                int patchno,
                                void *unpack_data_from_here)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);

    /* Time interp is false, since we only partition when all grids
       are time synchronized and all flux registers are set to 
       zero.  After copying data, we re-build patch with any 
       data needed.  */
    cp->griddata.copyFromMemory((double*)unpack_data_from_here);
}

/* ------------------------------------ Virtual table  -------------------------------- */

static
fclaw_clawpatch_vtable_t* clawpatch_vt_new(int dim)
{
    fclaw_clawpatch_vtable_t * vt = FCLAW_ALLOC_ZERO (fclaw_clawpatch_vtable_t, 1);
    if(dim == 2)
    {
        vt->patch_dim = 2;
        vt->d2 = FCLAW_ALLOC_ZERO (struct fclaw_clawpatch_vtable_2d, 1);
    }
    else
    {
        vt->patch_dim = 3;
        vt->d3 = FCLAW_ALLOC_ZERO (struct fclaw_clawpatch_vtable_3d, 1);
    }
    return vt;
}

static
void clawpatch_vt_destroy(void* vt)
{
    fclaw_clawpatch_vtable_t* clawpatch_vt = (fclaw_clawpatch_vtable_t*) vt;
    if(clawpatch_vt->patch_dim == 2)
    {
        FCLAW_FREE (clawpatch_vt->d2);
    }else{
        FCLAW_FREE (clawpatch_vt->d3);
    }
    FCLAW_FREE (clawpatch_vt);
}

static 
void initialize_2d_claw46_fort_vt(fclaw_clawpatch_vtable_t* clawpatch_vt)
{
    /* Clawpatch settings functions */
    clawpatch_vt->d2->fort_average2coarse        = FCLAW2D_CLAWPATCH46_FORT_AVERAGE2COARSE;
    clawpatch_vt->d2->fort_interpolate2fine      = FCLAW2D_CLAWPATCH46_FORT_INTERPOLATE2FINE;

    clawpatch_vt->d2->fort_tag4refinement        = FCLAW2D_CLAWPATCH46_FORT_TAG4REFINEMENT;
    clawpatch_vt->d2->fort_tag4coarsening        = FCLAW2D_CLAWPATCH46_FORT_TAG4COARSENING;

    /* output functions */
    clawpatch_vt->fort_header_ascii              = FCLAW2D_CLAWPATCH46_FORT_HEADER_ASCII;
    clawpatch_vt->d2->fort_output_ascii          = FCLAW2D_CLAWPATCH46_FORT_OUTPUT_ASCII;

    /* Diagnostic functions */
    clawpatch_vt->conservation_check             = fclaw_clawpatch_diagnostics_cons_default;
    clawpatch_vt->compute_error                  = fclaw_clawpatch_diagnostics_error_default;
    clawpatch_vt->d2->fort_compute_patch_error   = NULL;   /* User defined */
    clawpatch_vt->d2->fort_compute_error_norm    = FCLAW2D_CLAWPATCH46_FORT_COMPUTE_ERROR_NORM;
    clawpatch_vt->d2->fort_compute_patch_area    = FCLAW2D_CLAWPATCH46_FORT_COMPUTE_PATCH_AREA;
    clawpatch_vt->d2->fort_conservation_check    = FCLAW2D_CLAWPATCH46_FORT_CONSERVATION_CHECK;

    /* Ghost cell exchange functions */
    clawpatch_vt->d2->fort_copy_face             = FCLAW2D_CLAWPATCH46_FORT_COPY_FACE;
    clawpatch_vt->d2->fort_average_face          = FCLAW2D_CLAWPATCH46_FORT_AVERAGE_FACE;
    clawpatch_vt->d2->fort_interpolate_face      = FCLAW2D_CLAWPATCH46_FORT_INTERPOLATE_FACE;

    clawpatch_vt->d2->fort_copy_corner           = FCLAW2D_CLAWPATCH46_FORT_COPY_CORNER;
    clawpatch_vt->d2->fort_average_corner        = FCLAW2D_CLAWPATCH46_FORT_AVERAGE_CORNER;
    clawpatch_vt->d2->fort_interpolate_corner    = FCLAW2D_CLAWPATCH46_FORT_INTERPOLATE_CORNER;

    clawpatch_vt->local_ghost_pack_aux           = NULL;
    clawpatch_vt->d2->fort_local_ghost_pack      = FCLAW2D_CLAWPATCH46_FORT_LOCAL_GHOST_PACK;

    clawpatch_vt->d2->fort_timeinterp            = FCLAW2D_CLAWPATCH46_FORT_TIMEINTERP;
}

static 
void initialize_2d_claw5_fort_vt(fclaw_clawpatch_vtable_t* clawpatch_vt)
{
    /* Clawpatch settings functions */
    clawpatch_vt->d2->fort_average2coarse     = FCLAW2D_CLAWPATCH5_FORT_AVERAGE2COARSE;
    clawpatch_vt->d2->fort_interpolate2fine   = FCLAW2D_CLAWPATCH5_FORT_INTERPOLATE2FINE;

    clawpatch_vt->d2->fort_tag4refinement     = FCLAW2D_CLAWPATCH5_FORT_TAG4REFINEMENT;
    clawpatch_vt->d2->fort_tag4coarsening     = FCLAW2D_CLAWPATCH5_FORT_TAG4COARSENING;

    /* output functions */
    clawpatch_vt->fort_header_ascii           = FCLAW2D_CLAWPATCH5_FORT_HEADER_ASCII;
    clawpatch_vt->d2->fort_output_ascii       = FCLAW2D_CLAWPATCH5_FORT_OUTPUT_ASCII;

    /* Diagnostic functions */
    clawpatch_vt->conservation_check          = fclaw_clawpatch_diagnostics_cons_default;
    clawpatch_vt->compute_error               = fclaw_clawpatch_diagnostics_error_default;
    clawpatch_vt->d2->fort_compute_patch_error  = NULL;   /* User defined */
    clawpatch_vt->d2->fort_compute_error_norm = FCLAW2D_CLAWPATCH5_FORT_COMPUTE_ERROR_NORM;
    clawpatch_vt->d2->fort_compute_patch_area = FCLAW2D_CLAWPATCH5_FORT_COMPUTE_PATCH_AREA;
    clawpatch_vt->d2->fort_conservation_check = FCLAW2D_CLAWPATCH5_FORT_CONSERVATION_CHECK;

    /* Ghost cell exchange functions */
    clawpatch_vt->d2->fort_copy_face          = FCLAW2D_CLAWPATCH5_FORT_COPY_FACE;
    clawpatch_vt->d2->fort_average_face       = FCLAW2D_CLAWPATCH5_FORT_AVERAGE_FACE;
    clawpatch_vt->d2->fort_interpolate_face   = FCLAW2D_CLAWPATCH5_FORT_INTERPOLATE_FACE;

    clawpatch_vt->d2->fort_copy_corner        = FCLAW2D_CLAWPATCH5_FORT_COPY_CORNER;
    clawpatch_vt->d2->fort_average_corner     = FCLAW2D_CLAWPATCH5_FORT_AVERAGE_CORNER;
    clawpatch_vt->d2->fort_interpolate_corner = FCLAW2D_CLAWPATCH5_FORT_INTERPOLATE_CORNER;

    clawpatch_vt->local_ghost_pack_aux        = NULL;
    clawpatch_vt->d2->fort_local_ghost_pack   = FCLAW2D_CLAWPATCH5_FORT_LOCAL_GHOST_PACK;

    clawpatch_vt->d2->fort_timeinterp         = FCLAW2D_CLAWPATCH5_FORT_TIMEINTERP;
}

static 
void initialize_3dx_claw46_fort_vt(fclaw_clawpatch_vtable_t* clawpatch_vt)
{
    /* Clawpatch settings functions */
    clawpatch_vt->d3->fort_average2coarse        = FCLAW3DX_CLAWPATCH46_FORT_AVERAGE2COARSE;
    clawpatch_vt->d3->fort_interpolate2fine      = FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE2FINE;

    clawpatch_vt->d3->fort_tag4refinement        = FCLAW3D_CLAWPATCH46_FORT_TAG4REFINEMENT;
    clawpatch_vt->d3->fort_tag4coarsening_3dx        = FCLAW3DX_CLAWPATCH46_FORT_TAG4COARSENING;

    /* output functions */
    clawpatch_vt->fort_header_ascii              = FCLAW3D_CLAWPATCH46_FORT_HEADER_ASCII;
    clawpatch_vt->d3->fort_output_ascii          = FCLAW3D_CLAWPATCH46_FORT_OUTPUT_ASCII;

    /* Ghost cell exchange functions */
    clawpatch_vt->d3->fort_copy_face             = FCLAW3DX_CLAWPATCH46_FORT_COPY_FACE;
    clawpatch_vt->d3->fort_average_face          = FCLAW3DX_CLAWPATCH46_FORT_AVERAGE_FACE;
    clawpatch_vt->d3->fort_interpolate_face      = FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE_FACE;

    clawpatch_vt->d3->fort_copy_corner           = FCLAW3DX_CLAWPATCH46_FORT_COPY_CORNER;
    clawpatch_vt->d3->fort_average_corner        = FCLAW3DX_CLAWPATCH46_FORT_AVERAGE_CORNER;
    clawpatch_vt->d3->fort_interpolate_corner    = FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE_CORNER;

    clawpatch_vt->local_ghost_pack_aux           = NULL;
    clawpatch_vt->d3->fort_local_ghost_pack      = FCLAW3DX_CLAWPATCH46_FORT_LOCAL_GHOST_PACK;

    clawpatch_vt->d3->fort_timeinterp            = FCLAW3D_CLAWPATCH46_FORT_TIMEINTERP;
}

static 
void initialize_3d_claw46_fort_vt(fclaw_clawpatch_vtable_t* clawpatch_vt)
{
    /* Clawpatch settings functions */
    clawpatch_vt->d3->fort_average2coarse        = FCLAW3D_CLAWPATCH46_FORT_AVERAGE2COARSE;
    clawpatch_vt->d3->fort_interpolate2fine      = FCLAW3D_CLAWPATCH46_FORT_INTERPOLATE2FINE;

    clawpatch_vt->d3->fort_tag4refinement        = FCLAW3D_CLAWPATCH46_FORT_TAG4REFINEMENT;
    clawpatch_vt->d3->fort_tag4coarsening        = FCLAW3D_CLAWPATCH46_FORT_TAG4COARSENING;

    /* output functions */
    clawpatch_vt->fort_header_ascii              = FCLAW3D_CLAWPATCH46_FORT_HEADER_ASCII;
    clawpatch_vt->d3->fort_output_ascii          = FCLAW3D_CLAWPATCH46_FORT_OUTPUT_ASCII;

    /* Ghost cell exchange functions */
    clawpatch_vt->d3->fort_copy_face             = FCLAW3D_CLAWPATCH46_FORT_COPY_FACE;
    clawpatch_vt->d3->fort_average_face          = FCLAW3D_CLAWPATCH46_FORT_AVERAGE_FACE;
    clawpatch_vt->d3->fort_interpolate_face      = FCLAW3D_CLAWPATCH46_FORT_INTERPOLATE_FACE;

    clawpatch_vt->d3->fort_copy_edge             = FCLAW3D_CLAWPATCH46_FORT_COPY_EDGE;
    clawpatch_vt->d3->fort_average_edge          = FCLAW3D_CLAWPATCH46_FORT_AVERAGE_EDGE;
    clawpatch_vt->d3->fort_interpolate_edge      = FCLAW3D_CLAWPATCH46_FORT_INTERPOLATE_EDGE;

    clawpatch_vt->d3->fort_copy_corner           = FCLAW3D_CLAWPATCH46_FORT_COPY_CORNER;
    clawpatch_vt->d3->fort_average_corner        = FCLAW3D_CLAWPATCH46_FORT_AVERAGE_CORNER;
    clawpatch_vt->d3->fort_interpolate_corner    = FCLAW3D_CLAWPATCH46_FORT_INTERPOLATE_CORNER;

    clawpatch_vt->local_ghost_pack_aux           = NULL;
    clawpatch_vt->d3->fort_local_ghost_pack      = FCLAW3D_CLAWPATCH46_FORT_LOCAL_GHOST_PACK;

    clawpatch_vt->d3->fort_timeinterp            = FCLAW3D_CLAWPATCH46_FORT_TIMEINTERP;
}

static
void initialize_2d_patch_vt(fclaw_patch_vtable_t* patch_vt, int patch_dim)
{
    /* Ghost filling */
    patch_vt->copy_face            = clawpatch_copy_face;
    patch_vt->average_face         = clawpatch_average_face;
    patch_vt->interpolate_face     = clawpatch_interpolate_face;

    patch_vt->copy_corner          = clawpatch_copy_corner;
    patch_vt->average_corner       = clawpatch_average_corner;
    patch_vt->interpolate_corner   = clawpatch_interpolate_corner;

    /* Assume regular block corners;  Change by calling 'fclaw2d_clawpatch_use_pillowsphere' */    
    patch_vt->copy_block_corner          = clawpatch_copy_corner;
    patch_vt->average_block_corner       = clawpatch_average_corner;
    patch_vt->interpolate_block_corner   = clawpatch_interpolate_corner;

    if(patch_dim == 2)
    {
        /* Timing syncing module for conservation */
        patch_vt->time_sync_f2c       = fclaw2d_clawpatch_time_sync_f2c;
        patch_vt->time_sync_samesize  = fclaw2d_clawpatch_time_sync_samesize;
        patch_vt->time_sync_reset     = fclaw2d_clawpatch_time_sync_reset;
    }

    patch_vt->transform_init_data  = fclaw2d_clawpatch_transform_init_data;
    patch_vt->transform_face       = fclaw2d_clawpatch_face_transformation;
    patch_vt->transform_face_intra = fclaw2d_clawpatch_face_transformation_intra;
    
}

static
void initialize_3d_patch_vt(fclaw_patch_vtable_t* patch_vt)
{
    /* Ghost filling */
    patch_vt->copy_face            = clawpatch_copy_face;
    patch_vt->average_face         = clawpatch_average_face;
    patch_vt->interpolate_face     = clawpatch_interpolate_face;

    patch_vt->copy_edge            = clawpatch_copy_edge;
    patch_vt->average_edge         = clawpatch_average_edge;
    patch_vt->interpolate_edge     = clawpatch_interpolate_edge;

    patch_vt->copy_corner          = clawpatch_copy_corner;
    patch_vt->average_corner       = clawpatch_average_corner;
    patch_vt->interpolate_corner   = clawpatch_interpolate_corner;

    /* Assume regular block corners;  Change by calling 'fclaw2d_clawpatch_use_pillowsphere' */    
    patch_vt->copy_block_corner          = clawpatch_copy_corner;
    patch_vt->average_block_corner       = clawpatch_average_corner;
    patch_vt->interpolate_block_corner   = clawpatch_interpolate_corner;

    patch_vt->transform_init_data  = fclaw3d_clawpatch_transform_init_data;
    patch_vt->transform_face       = fclaw3d_clawpatch_face_transformation;
    patch_vt->transform_face_intra = fclaw3d_clawpatch_face_transformation_intra;

}

void fclaw_clawpatch_vtable_initialize(fclaw_global_t* glob, 
                                       int claw_version)
{
    fclaw_domain_t* domain = glob->domain;
    fclaw_clawpatch_options_t * clawpatch_opt = 
                             fclaw_clawpatch_get_options(glob);
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);

    if(clawpatch_opt->patch_dim == 2 && domain->refine_dim == 3)
    {
        fclaw_abortf("fclaw_clawpatch_vtable_initialize : " \
                    "2d patches not supported with 3d oct-tree refinement\n");
    }

    if(clawpatch_opt->patch_dim == 2)
    {
        fclaw2d_metric_vtable_initialize(glob);
    }
    else 
    {
        fclaw3d_metric_vtable_initialize(glob);

    }

    fclaw_clawpatch_vtable_t *clawpatch_vt = clawpatch_vt_new(clawpatch_opt->patch_dim);

    /* Patch setup */
    patch_vt->patch_new             = (clawpatch_vt->patch_dim == 2) ? clawpatch_new_2d : clawpatch_new_3d;
    patch_vt->patch_delete          = clawpatch_delete;
    patch_vt->build                 = clawpatch_build;
    patch_vt->build_from_fine       = clawpatch_build_from_fine;    

    /* Time stepping */
    patch_vt->restore_step          = clawpatch_restore_step;
    patch_vt->save_step             = clawpatch_save_step;
    patch_vt->setup_timeinterp      = clawpatch_setup_timeinterp;

    if(domain->refine_dim == 2)
    {
        initialize_2d_patch_vt(patch_vt, clawpatch_opt->patch_dim);
    }
    else 
    {
        initialize_3d_patch_vt(patch_vt);
    }
    if(clawpatch_opt->patch_dim == 2)
    {
        patch_vt->metric_patch = clawpatch_get_metric_patch_2d;
    }
    else 
    {
        patch_vt->metric_patch = clawpatch_get_metric_patch_3d;
    }

    /* Regridding  functions */
    patch_vt->tag4refinement       = clawpatch_tag4refinement;
    patch_vt->tag4coarsening       = clawpatch_tag4coarsening;

    patch_vt->average2coarse       = clawpatch_average2coarse;
    patch_vt->interpolate2fine     = clawpatch_interpolate2fine;

    /* ghost patch */
    patch_vt->ghost_packsize       = clawpatch_ghost_packsize;
    patch_vt->local_ghost_pack     = clawpatch_local_ghost_pack;
    patch_vt->remote_ghost_build   = clawpatch_remote_ghost_build;
    patch_vt->remote_ghost_unpack  = clawpatch_remote_ghost_unpack;
    patch_vt->remote_ghost_delete  = clawpatch_remote_ghost_delete;


    /* partitioning */
    patch_vt->partition_packsize   = clawpatch_partition_packsize;
    patch_vt->partition_pack       = clawpatch_partition_pack;
    patch_vt->partition_unpack     = clawpatch_partition_unpack;

    /* output functions */
    clawpatch_vt->time_header_ascii  = fclaw_clawpatch_time_header_ascii;
    clawpatch_vt->cb_output_ascii    = cb_clawpatch_output_ascii; 

    if(clawpatch_opt->patch_dim == 2)
    {
        /* Ghost pack for registers (doesn't depend on clawpack version) */
        clawpatch_vt->time_sync_pack_registers = fclaw2d_clawpatch_time_sync_pack_registers;
    }

    /* Fortran functions that depend on data layout (version 4.6 or 5.0) */

    if (clawpatch_vt->patch_dim == 2 && claw_version == 4)
    {
        initialize_2d_claw46_fort_vt(clawpatch_vt);
    }
    else if (clawpatch_vt->patch_dim == 2 && claw_version == 5)
    {
        initialize_2d_claw5_fort_vt(clawpatch_vt);
    }
    else if (domain->refine_dim == 2 && clawpatch_vt->patch_dim == 3 && claw_version == 4)
    {
        initialize_3dx_claw46_fort_vt(clawpatch_vt);
    }
    else if (domain->refine_dim == 2 && clawpatch_vt->patch_dim == 3 && claw_version == 5)
    {
        fclaw_global_essentialf("clawpatch_vtable_initialize : Version 5 not yet " \
                                "implemented\n");
        exit(0);
    }
    else if (domain->refine_dim == 3 && clawpatch_vt->patch_dim == 3 && claw_version == 4)
    {
        initialize_3d_claw46_fort_vt(clawpatch_vt);
    }
    else if (domain->refine_dim == 3 && clawpatch_vt->patch_dim == 3 && claw_version == 5)
    {
        fclaw_global_essentialf("clawpatch_vtable_initialize : Version 5 not yet " \
                                "implemented\n");
        exit(0);
    }


    fclaw_clawpatch_diagnostics_vtable_initialize(glob);

    /* Set the virtual table, even if it isn't used */
    fclaw_clawpatch_pillow_vtable_initialize(glob, claw_version);

    clawpatch_vt->is_set = 1;

    if(fclaw_pointer_map_get(glob->vtables, "fclaw_clawpatch") != NULL)
    {
        fclaw_abortf("fclaw_clawpatch_vtable_initialize : " \
                    "fclaw_clawpatch already initialized\n");
    }
    fclaw_pointer_map_insert(glob->vtables, "fclaw_clawpatch", clawpatch_vt, clawpatch_vt_destroy);
}

/* ------------------------------- Public access functions ---------------------------- */

/* These functions are not virtualized and are not defined by the 
   patch interface */

int fclaw_clawpatch_dim(fclaw_patch_t* patch)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    return cp->patch_dim;
}

fclaw_clawpatch_vtable_t* fclaw_clawpatch_vt(fclaw_global_t* glob)
{

    fclaw_clawpatch_vtable_t* clawpatch_vt = (fclaw_clawpatch_vtable_t*) 
                          fclaw_pointer_map_get(glob->vtables, "fclaw_clawpatch");
    FCLAW_ASSERT(clawpatch_vt != nullptr);
    FCLAW_ASSERT(clawpatch_vt->is_set != 0);
    return clawpatch_vt;
}

/* Called from clawpack 4.6 and 5.0 */
void fclaw_clawpatch_save_current_step(fclaw_global_t* glob,
                                       fclaw_patch_t* patch)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    cp->griddata_last = cp->griddata;
}


fclaw_clawpatch_t* 
fclaw_clawpatch_get_clawpatch(fclaw_patch_t* patch)
{
    return get_clawpatch(patch);
}


fclaw2d_metric_patch_t* 
fclaw_clawpatch_get_2d_metric_patch(fclaw_patch_t* patch)
{
    return get_metric_patch_2d(patch);
}

fclaw3d_metric_patch_t* 
fclaw3d_clawpatch_get_metric_patch(fclaw_patch_t* patch)
{
    return get_metric_patch_3d(patch);
}

void fclaw_clawpatch_2d_grid_data(fclaw_global_t* glob,
                                 fclaw_patch_t* patch,
                                 int* mx, int* my, int* mbc,
                                 double* xlower, double* ylower,
                                 double* dx, double* dy)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    *mx = cp->mx;
    *my = cp->my;
    *mbc = cp->mbc;
    *xlower = cp->xlower;
    *ylower = cp->ylower;
    *dx = cp->dx;
    *dy = cp->dy;
}

void fclaw_clawpatch_3d_grid_data(fclaw_global_t* glob,
                                  fclaw_patch_t* patch,
                                  int* mx, int* my, int* mz, int* mbc,
                                  double* xlower, double* ylower,
                                  double* zlower, 
                                  double* dx, double* dy, double* dz)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    *mx = cp->mx;
    *my = cp->my;
    *mz = cp->mz;
    *mbc = cp->mbc;
    *xlower = cp->xlower;
    *ylower = cp->ylower;
    *zlower = cp->zlower;
    *dx = cp->dx;
    *dy = cp->dy;
    *dz = cp->dz;
}

void fclaw_clawpatch_aux_data(fclaw_global_t *glob,
                                fclaw_patch_t *patch,
                                double **aux, int* maux)
{
    fclaw_clawpatch_t *cp = get_clawpatch (patch);

    *maux = cp->maux;
    *aux = cp->aux.dataPtr();
}

void fclaw_clawpatch_soln_data(fclaw_global_t* glob,
                                 fclaw_patch_t* patch,
                                 double **q, int* meqn)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    *q = cp->griddata.dataPtr();
    *meqn = cp->meqn;
}

void fclaw_clawpatch_rhs_data(fclaw_global_t* glob,
                                 fclaw_patch_t* patch,
                                 double **rhs, int *mfields)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    *rhs = cp->rhs.dataPtr();
    *mfields = cp->mfields;
}

void fclaw_clawpatch_elliptic_error_data(fclaw_global_t* glob,
                                         fclaw_patch_t* patch,
                                         double **err, int *mfields)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    *err = cp->elliptic_error.dataPtr();
    *mfields = cp->mfields;
}

void fclaw_clawpatch_elliptic_soln_data(fclaw_global_t* glob,
                                        fclaw_patch_t* patch,
                                        double **soln, int *mfields)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    *soln = cp->elliptic_soln.dataPtr();
    *mfields = cp->mfields;
}


double *fclaw_clawpatch_get_q(fclaw_global_t* glob,
                              fclaw_patch_t* patch)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    return cp->griddata.dataPtr();
}

fclaw2d_clawpatch_registers_t* 
fclaw_clawpatch_get_2d_registers(fclaw_global_t* glob,
                                fclaw_patch_t* patch)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    return cp->d2->registers;
}


double* fclaw_clawpatch_get_error(fclaw_global_t* glob,
                                  fclaw_patch_t* patch)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    return cp->griderror.dataPtr();
}

double* fclaw_clawpatch_get_exactsoln(fclaw_global_t* glob,
                                      fclaw_patch_t* patch)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    return cp->exactsolution.dataPtr();
}

void* fclaw_clawpatch_get_user_data(fclaw_global_t* glob,
                                    fclaw_patch_t* patch)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    return cp->user_data;
}

void fclaw_clawpatch_set_user_data(fclaw_global_t* glob,
                                   fclaw_patch_t* patch,
                                   void *udata)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    cp->user_data = udata;
}

void* fclaw_clawpatch_get_solver_data(fclaw_global_t* glob,
                                      fclaw_patch_t* patch)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    return cp->solver_data;
}

void fclaw_clawpatch_set_solver_data(fclaw_global_t* glob,
                                     fclaw_patch_t* patch,
                                     void *sdata)
{
    fclaw_clawpatch_t *cp = get_clawpatch(patch);
    cp->solver_data = sdata;
}

size_t fclaw_clawpatch_size(fclaw_global_t *glob)
{
    const fclaw_clawpatch_options_t *clawpatch_opt = 
                     fclaw_clawpatch_get_options(glob);
    int meqn = clawpatch_opt->meqn;
    int mbc = clawpatch_opt->mbc;

    if(clawpatch_opt->patch_dim == 2)
    {
        int mx = clawpatch_opt->mx;
        int my = clawpatch_opt->my;
        return (mx+2*mbc)*(my+2*mbc)*meqn;
    }
    else 
    {
        int mx = clawpatch_opt->mx;
        int my = clawpatch_opt->my;
        int mz = clawpatch_opt->mz;
        return (mx+2*mbc)*(my+2*mbc)*(mz+2*mbc)*meqn;
    }
}

void fclaw_clawpatch_timesync_data(fclaw_global_t* glob,
                                   fclaw_patch_t* patch,
                                   int time_interp,
                                   double **q, int* meqn)
{
    fclaw_clawpatch_t* cp = get_clawpatch(patch);
    *q = q_time_sync(patch, time_interp);
    *meqn = cp->meqn;
}

double* fclaw_clawpatch_get_q_timesync(fclaw_global_t* glob,
                                       fclaw_patch_t* patch,
                                       int time_interp)
{
    return q_time_sync(patch, time_interp);
}

double* fclaw_clawpatch_get_2d_area(fclaw_global_t* glob,
                                 fclaw_patch_t* patch)
{
    return clawpatch_get_area(glob, patch);
}

double* fclaw_clawpatch_get_3d_volume(fclaw_global_t* glob,
                                     fclaw_patch_t* patch)
{
    return clawpatch_get_volume(glob, patch);
}

void fclaw_clawpatch_2d_metric_scalar(fclaw_global_t* glob,
                                     fclaw_patch_t* patch,
                                     double **area, double** edgelengths,
                                     double **curvature)
{
    fclaw2d_metric_patch_scalar(glob,patch,area,edgelengths,
                                curvature);
}

void fclaw_clawpatch_2d_metric_vector(struct fclaw_global* glob,
                                     struct fclaw_patch* patch,
                                     double **xnormals, double **ynormals,
                                     double **xtangents, double **ytangents,
                                     double **surfnormals)
{
    fclaw2d_metric_patch_vector(glob,patch,xnormals,ynormals,
                                xtangents,ytangents,surfnormals);
}




void fclaw_clawpatch_2d_metric_data(fclaw_global_t* glob,
                                   fclaw_patch_t* patch,
                                   double **xp, double **yp, double **zp,
                                   double **xd, double **yd, double **zd,
                                   double **area)
{
    fclaw2d_metric_patch_mesh_data(glob,patch,xp,yp,zp,xd,yd,zd,area);
}

void fclaw_clawpatch_2d_metric_data2(fclaw_global_t* glob,
                                    fclaw_patch_t* patch,
                                    double **xnormals, double **ynormals,
                                    double **xtangents, double **ytangents,
                                    double **surfnormals,
                                    double **edgelengths, double **curvature)
{
    fclaw2d_metric_patch_mesh_data2(glob,patch,xnormals,ynormals,
                                    xtangents,ytangents,surfnormals,
                                    edgelengths,curvature);
}

void fclaw_clawpatch_3d_metric_scalar(fclaw_global_t* glob,
                                     fclaw_patch_t* patch,
                                     double **volume, double** faceareas)
{
    fclaw3d_metric_patch_scalar(glob,patch,volume,faceareas);
}

void fclaw2d_clawpatch_metric_basis(struct fclaw_global* glob,
                                     struct fclaw_patch* patch,
                                     double **xrot, double** yrot, double** zrot)
{
    fclaw3d_metric_patch_basis(glob,patch,xrot, yrot, zrot);
}




void fclaw_clawpatch_3d_mesh_data(fclaw_global_t* glob,
                                 fclaw_patch_t* patch,
                                 double **xp, double **yp, double **zp,
                                 double **xd, double **yd, double **zd,
                                 double **volume, double **faceareas)
{
    fclaw3d_metric_patch_mesh_data(glob,patch,xp,yp,zp,xd,yd,zd,
                                   volume, faceareas);
}

