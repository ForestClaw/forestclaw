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

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch.hpp>

#include <fclaw2d_clawpatch_diagnostics.h>
#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_clawpatch_output_ascii.h> 
#include <fclaw2d_clawpatch_output_vtk.h>
#include <fclaw2d_clawpatch_fort.h>


#include <fclaw2d_patch.h>  /* Needed to get enum for build modes */

#include <fclaw2d_defs.h>
#include <fclaw2d_global.h>

#include <fclaw2d_global.h>
#include <fclaw2d_vtable.h>
#include <fclaw2d_options.h>
#include <fclaw2d_transform.h>

#include <fclaw2d_timeinterp.h>
#include <fclaw2d_diagnostics.h>

#include <fclaw2d_metric.h>
#include <fclaw2d_map_query.h>

/* Needed only for MB_BLOCK_CORNER_GHOST */
#include <fclaw2d_neighbors_fort.h>


/* ------------------------------- Static function defs ------------------------------- */

static fclaw2d_clawpatch_vtable_t s_clawpatch_vt;

static
fclaw2d_clawpatch_vtable_t* clawpatch_vt()
{
    FCLAW_ASSERT(s_clawpatch_vt.is_set != 0);
    return &s_clawpatch_vt;
}

static
fclaw2d_clawpatch_t* clawpatch_data(fclaw2d_patch_t *this_patch)
{
    fclaw2d_clawpatch_t *cp = (fclaw2d_clawpatch_t*) 
                     fclaw2d_patch_get_user_patch(this_patch);
    return cp;
}

static
double* q_time_sync(fclaw2d_clawpatch_t* cp, int time_interp);

/* ----------------------------- Creating/deleting patches ---------------------------- */

static
void* clawpatch_new()
{
    fclaw2d_clawpatch_t *cp = new fclaw2d_clawpatch_t;
    return (void*) cp;
}

static
void clawpatch_delete(void *patchcp)
{
    FCLAW_ASSERT(patchcp != NULL);
    fclaw2d_clawpatch_t* cp = (fclaw2d_clawpatch_t*) patchcp;
    delete cp;
    patchcp = NULL;
}


static
void setup_area_storage(fclaw2d_clawpatch_t* cp);

static
void setup_metric_storage(fclaw2d_clawpatch_t* cp);


static
void metric_setup(fclaw2d_global_t* glob,
                  fclaw2d_patch_t* this_patch,
                  int blockno,
                  int patchno);

/* Maybe this should just be a 'build' function? */
static
void clawpatch_define(fclaw2d_global_t* glob,
                      fclaw2d_patch_t *this_patch,
                      int blockno, int patchno,
                      fclaw2d_build_mode_t build_mode)
{
    /* We are getting closer to getting rid the class fclaw2d_clawpatch_t */
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);

    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    cp->mx = clawpatch_opt->mx;
    cp->my = clawpatch_opt->my;
    cp->mbc = clawpatch_opt->mbc;
    cp->blockno = blockno;
    cp->meqn = clawpatch_opt->meqn;
    cp->maux = clawpatch_opt->maux;

    for (int icorner=0; icorner < 4; icorner++)
    {
        fclaw2d_patch_set_block_corner_count(glob,this_patch,icorner,0);
    }

    fclaw2d_map_context_t* cont = glob->cont;

    int is_brick = FCLAW2D_MAP_IS_BRICK(&cont);

    cp->manifold = fclaw_opt->manifold;

    if (cp->manifold)
    {
        cp->xlower = this_patch->xlower;
        cp->ylower = this_patch->ylower;
        cp->xupper = this_patch->xupper;
        cp->yupper = this_patch->yupper;
    }
    else
    {
        double ax = fclaw_opt->ax;
        double bx = fclaw_opt->bx;
        double ay = fclaw_opt->ay;
        double by = fclaw_opt->by;

        double xl = this_patch->xlower;
        double yl = this_patch->ylower;
        double xu = this_patch->xupper;
        double yu = this_patch->yupper;

        double xlower, ylower, xupper, yupper;

        if (is_brick)
        {
            double z;
            /* Scale to [0,1]x[0,1], based on blockno */
            fclaw2d_map_c2m_nomap_brick(cont,cp->blockno,xl,yl,&xlower,&ylower,&z);
            fclaw2d_map_c2m_nomap_brick(cont,cp->blockno,xu,yu,&xupper,&yupper,&z);
        }
        else
        {
            xlower = xl;
            ylower = yl;
            xupper = xu;
            yupper = yu;
        }

        cp->xlower = ax + (bx - ax)*xlower;
        cp->xupper = ax + (bx - ax)*xupper;
        cp->ylower = ay + (by - ay)*ylower;
        cp->yupper = ay + (by - ay)*yupper;
    }

    cp->dx = (cp->xupper - cp->xlower)/cp->mx;
    cp->dy = (cp->yupper - cp->ylower)/cp->my;

    int ll[2];
    int ur[2];
    for (int idir = 0; idir < 2; idir++)
    {
        ll[idir] = 1-cp->mbc;
    }
    ur[0] = cp->mx + cp->mbc;
    ur[1] = cp->my + cp->mbc;
    Box box(ll,ur);

    // This will destroy any existing memory n griddata.
    cp->griddata.define(box, cp->meqn);
    if (fclaw_opt->subcycle)
    {
        cp->griddata_time_interpolated.define(box, cp->meqn);
    }
    if (fclaw_opt->compute_error)
    {
        cp->griderror.define(box,cp->meqn);
    }

    if (clawpatch_opt->maux > 0)
    {
      cp->aux.define(box,cp->maux);
    }
    // Set up storage for metric terms, if needed.
    if (fclaw_opt->manifold)
    {
        setup_area_storage(cp);
        if (build_mode != FCLAW2D_BUILD_FOR_GHOST_AREA_PACKED)
        {
            /* Don't need any more manifold info for ghost patches */
            if (build_mode == FCLAW2D_BUILD_FOR_UPDATE)
            {
                setup_metric_storage(cp);
            }
        }
    }
    
    if (build_mode != FCLAW2D_BUILD_FOR_UPDATE)
    {
        return;
    }

    cp->griddata_last.define(box, cp->meqn);
    cp->griddata_save.define(box, cp->meqn);

}

static
void clawpatch_build(fclaw2d_global_t *glob,
                     fclaw2d_patch_t *this_patch,
                     int blockno,
                     int patchno,
                     void *user)
{
    fclaw2d_build_mode_t build_mode =  *((fclaw2d_build_mode_t*) user);
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

    clawpatch_define(glob,this_patch,blockno,patchno,build_mode);

    if (fclaw_opt->manifold)
    {
        fclaw2d_vt()->metric_compute_area(glob,this_patch,blockno,patchno);
        metric_setup(glob,this_patch,blockno,patchno);
    }
}

static
void clawpatch_build_from_fine(fclaw2d_global_t *glob,
                               fclaw2d_patch_t *fine_patches,
                               fclaw2d_patch_t *coarse_patch,
                               int blockno,
                               int coarse_patchno,
                               int fine0_patchno,
                               fclaw2d_build_mode_t build_mode)
{
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

    clawpatch_define(glob,coarse_patch,blockno,coarse_patchno,build_mode);

    if (fclaw_opt->manifold)
    {
        /* Don't recompute the area, but rather average from finer areas */
        fclaw2d_metric_average_area(glob,fine_patches,coarse_patch,
                                    blockno, coarse_patchno, fine0_patchno);

        metric_setup(glob,coarse_patch,blockno,coarse_patchno);
    }
}

static
void metric_setup(fclaw2d_global_t* glob,
                  fclaw2d_patch_t* this_patch,
                  int blockno,
                  int patchno)
{
    /* vt.patch_manifold_setup_mesh(...) */
    fclaw2d_metric_setup_mesh(glob,this_patch,blockno,patchno);

    /* vt.patch_manifold_compute_normals(...) */
    fclaw2d_metric_compute_normals(glob,this_patch,blockno,patchno);
}

static void setup_area_storage(fclaw2d_clawpatch_t* cp)
{
    int mx = cp->mx;
    int my = cp->my;
    int mbc = cp->mbc;

    int ll[2];
    int ur[2];
    for (int idir = 0; idir < 2; idir++)
    {
        ll[idir] = -mbc;
    }
    ur[0] = mx + mbc + 1;
    ur[1] = my + mbc + 1;

    Box box_p(ll,ur);
    cp->area.define(box_p,1);
}

static 
void setup_metric_storage(fclaw2d_clawpatch_t* cp)
{
    /* 24 additional field variables needed for all metric terms
       xp,yp,zp           : 3
       xd,yd,zd           : 3
       surf_normals       : 3
       curvature          : 1
       <xy>face normals   : 6
       <xy>face tangents  : 6
       edge lengths       : 2
       -----------------------
       Total              : 24

       We should come up with a way to store only what is needed */

    int mx = cp->mx;
    int my = cp->my;
    int mbc = cp->mbc;

    int ll[2];
    int ur[2];
    for (int idir = 0; idir < 2; idir++)
    {
        ll[idir] = -cp->mbc;
    }
    ur[0] = mx + mbc + 1;
    ur[1] = my + mbc + 1;

    Box box_p(ll,ur);   /* Store cell centered values here */

    /* Mesh cell centers of physical mesh */
    cp->xp.define(box_p,1);
    cp->yp.define(box_p,1);
    cp->zp.define(box_p,1);
    cp->surf_normals.define(box_p,3);
    cp->curvature.define(box_p,1);

    /* Node centered values */
    for (int idir = 0; idir < 2; idir++)
    {
        ll[idir] = -mbc;
    }
    ur[0] = mx + mbc + 2;
    ur[1] = my + mbc + 2;
    Box box_d(ll,ur);

    cp->xd.define(box_d,1);
    cp->yd.define(box_d,1);
    cp->zd.define(box_d,1);

    /* Face centered values */
    cp->xface_normals.define(box_d,3);
    cp->yface_normals.define(box_d,3);
    cp->xface_tangents.define(box_d,3);
    cp->yface_tangents.define(box_d,3);
    cp->edge_lengths.define(box_d,2);
}

/* -------------------------------- time stepping ------------------------------------- */

static
void clawpatch_save_step(fclaw2d_global_t* glob,
                         fclaw2d_patch_t* this_patch)
{
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    cp->griddata_save = cp->griddata;
}


static
void clawpatch_restore_step(fclaw2d_global_t* glob,
                            fclaw2d_patch_t* this_patch)
{
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    cp->griddata = cp->griddata_save;
}

static
void clawpatch_setup_timeinterp(fclaw2d_global_t *glob,
                                fclaw2d_patch_t *this_patch,
                                double alpha)
{
    /* We use the pack size here to make sure we are setting
       everything correctly;  it isn't needed for memory
       allocation */
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    int mx = clawpatch_opt->mx;
    int my = clawpatch_opt->my;
    int meqn = clawpatch_opt->meqn;
    int mbc = clawpatch_opt->mbc;
    int mint = clawpatch_opt->interp_stencil_width/2+1;  

    int hole = (mx - 2*mint)*(my - 2*mint);  /* Hole in center */
    FCLAW_ASSERT(hole >= 0);

    int wg = mx*my;  /* Whole grid but no ghost cells.  
                        Ghost cells will be averaged from finer level. */
    int psize = (wg - hole)*meqn;
    FCLAW_ASSERT(psize > 0);

    /* Store time interpolated data that will be use in coarse grid
       exchanges */
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
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

    clawpatch_vt->fort_timeinterp(&mx,&my,&mbc,&meqn,&psize,
                                    qcurr,qlast,qinterp,&alpha,&ierror);

}


/* ------------------------------------- Ghost filling  ------------------------------- */

static
void clawpatch_copy_face(fclaw2d_global_t *glob,
                         fclaw2d_patch_t *this_patch,
                         fclaw2d_patch_t *neighbor_patch,
                         int iface,
                         int time_interp,
                         fclaw2d_transform_data_t *transform_data)

{
    int meqn,mx,my,mbc;
    double *qthis, *qneighbor;
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    mx = clawpatch_opt->mx;
    my = clawpatch_opt->my;
    mbc = clawpatch_opt->mbc;

    fclaw2d_clawpatch_timesync_data(glob,this_patch,time_interp,&qthis,&meqn);
    fclaw2d_clawpatch_timesync_data(glob,neighbor_patch,time_interp,&qneighbor,&meqn);

    clawpatch_vt()->fort_copy_face(&mx,&my,&mbc,&meqn,qthis,qneighbor,&iface,&transform_data);
}

static
void clawpatch_average_face(fclaw2d_global_t *glob,
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
    double *qcoarse, *qfine;
    double *areacoarse, *areafine;

    const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    fclaw2d_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
    qfine = fclaw2d_clawpatch_get_q(glob,fine_patch);

    /* These will be empty for non-manifolds cases */
    areacoarse = fclaw2d_clawpatch_get_area(glob,coarse_patch);
    areafine = fclaw2d_clawpatch_get_area(glob,fine_patch);

    mx = clawpatch_opt->mx;
    my = clawpatch_opt->my;
    mbc = clawpatch_opt->mbc;

    int manifold = fclaw_opt->manifold;
    clawpatch_vt()->fort_average_face(&mx,&my,&mbc,&meqn,qcoarse,qfine,areacoarse,areafine,
                                      &idir,&iface_coarse, &p4est_refineFactor, &refratio,
                                      &igrid,&manifold,&transform_data);


}

static
void clawpatch_interpolate_face(fclaw2d_global_t *glob,
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
    int meqn,mx,my,mbc;
    double *qcoarse, *qfine;
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    fclaw2d_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
    qfine = fclaw2d_clawpatch_get_q(glob,fine_patch);

    mx = clawpatch_opt->mx;
    my = clawpatch_opt->my;
    mbc = clawpatch_opt->mbc;

    clawpatch_vt()->fort_interpolate_face(&mx,&my,&mbc,&meqn,qcoarse,qfine,&idir,&iside,
                                          &p4est_refineFactor,&refratio,&igrid,&transform_data);
}

static
void clawpatch_copy_corner(fclaw2d_global_t *glob,
                           fclaw2d_patch_t *this_patch,
                           fclaw2d_patch_t *corner_patch,
                           int icorner,
                           int time_interp,
                           fclaw2d_transform_data_t *transform_data)
{
    int meqn,mx,my,mbc;
    double *qthis, *qcorner;
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    mx = clawpatch_opt->mx;
    my = clawpatch_opt->my;
    mbc = clawpatch_opt->mbc;

    fclaw2d_clawpatch_timesync_data(glob,this_patch,time_interp,&qthis,&meqn);
    fclaw2d_clawpatch_timesync_data(glob,corner_patch,time_interp,&qcorner,&meqn);

    clawpatch_vt()->fort_copy_corner(&mx,&my,&mbc,&meqn,qthis,qcorner,&icorner,&transform_data);

}

static
void clawpatch_average_corner(fclaw2d_global_t *glob,
                              fclaw2d_patch_t *coarse_patch,
                              fclaw2d_patch_t *fine_patch,
                              int coarse_corner,
                              int refratio,
                              int time_interp,
                              fclaw2d_transform_data_t* transform_data)
{
    int meqn,mx,my,mbc;
    double *qcoarse, *qfine;
    double *areacoarse, *areafine;

    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    fclaw2d_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
    qfine = fclaw2d_clawpatch_get_q(glob,fine_patch);

    areacoarse = fclaw2d_clawpatch_get_area(glob,coarse_patch);
    areafine = fclaw2d_clawpatch_get_area(glob,fine_patch);

    mx = clawpatch_opt->mx;
    my = clawpatch_opt->my;
    mbc = clawpatch_opt->mbc;

    int manifold = fclaw_opt->manifold;
    clawpatch_vt()->fort_average_corner(&mx,&my,&mbc,&meqn,&refratio,
                                        qcoarse,qfine,areacoarse,areafine,
                                        &manifold,&coarse_corner,&transform_data);
}

static
void clawpatch_interpolate_corner(fclaw2d_global_t* glob,
                                  fclaw2d_patch_t* coarse_patch,
                                  fclaw2d_patch_t* fine_patch,
                                  int coarse_corner,
                                  int refratio,
                                  int time_interp,
                                  fclaw2d_transform_data_t* transform_data)

{
    int meqn,mx,my,mbc;
    double *qcoarse, *qfine;
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    mx = clawpatch_opt->mx;
    my = clawpatch_opt->my;
    mbc = clawpatch_opt->mbc;
    meqn = clawpatch_opt->meqn;


    fclaw2d_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
    qfine = fclaw2d_clawpatch_get_q(glob,fine_patch);

    clawpatch_vt()->fort_interpolate_corner(&mx,&my,&mbc,&meqn,
                                            &refratio,qcoarse,qfine,
                                            &coarse_corner,&transform_data);
}

/* -------------------------------- Regridding functions ------------------------------ */

static
int clawpatch_tag4refinement(fclaw2d_global_t *glob,
                             fclaw2d_patch_t *this_patch,
                             int blockno, int patchno,
                             int initflag)
    {
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *q;
    int tag_patch;
    double refine_threshold;

    refine_threshold = fclaw_opt->refine_threshold;

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);

    tag_patch = 0;
    fclaw2d_clawpatch_vt()->fort_tag4refinement(&mx,&my,&mbc,&meqn,&xlower,&ylower,&dx,&dy,
                                                &blockno, q,&refine_threshold,
                                                &initflag,&tag_patch);
    return tag_patch;
}

static
int clawpatch_tag4coarsening(fclaw2d_global_t *glob,
                             fclaw2d_patch_t *fine_patches,
                             int blockno,
                             int patchno)
{
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *q[4];
    int tag_patch,igrid;
    double coarsen_threshold;
    fclaw2d_patch_t *patch0;

    patch0 = &fine_patches[0];

    coarsen_threshold = fclaw_opt->coarsen_threshold;

    fclaw2d_clawpatch_grid_data(glob,patch0,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    for (igrid = 0; igrid < 4; igrid++)
    {
        fclaw2d_clawpatch_soln_data(glob,&fine_patches[igrid],&q[igrid],&meqn);
    }

    tag_patch = 0;
    fclaw2d_clawpatch_vt()->fort_tag4coarsening(&mx,&my,&mbc,&meqn,&xlower,&ylower,&dx,&dy,
                                                &blockno, q[0],q[1],q[2],q[3],
                                                &coarsen_threshold,&tag_patch);
    return tag_patch == 1;
}

static
void clawpatch_interpolate2fine(fclaw2d_global_t* glob,
                                fclaw2d_patch_t *coarse_patch,
                                fclaw2d_patch_t* fine_patches,
                                int this_blockno, int coarse_patchno,
                                int fine0_patchno)
{
    int mx,my,mbc,meqn;
    double *qcoarse,*qfine;
    double *areacoarse,*areafine;
    double *xp,*yp,*zp,*xd,*yd,*zd;
    int igrid;

    const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    fclaw2d_patch_t* fine_patch;

    mx = clawpatch_opt->mx;
    my = clawpatch_opt->my;
    mbc = clawpatch_opt->mbc;

    fclaw2d_clawpatch_metric_data(glob,coarse_patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&areacoarse);
    fclaw2d_clawpatch_soln_data(glob,coarse_patch,&qcoarse,&meqn);

    /* Loop over four siblings (z-ordering) */
    for (igrid = 0; igrid < 4; igrid++)
    {
        fine_patch = &fine_patches[igrid];

        fclaw2d_clawpatch_soln_data(glob,fine_patch,&qfine,&meqn);

        if (fclaw_opt->manifold)
        {
            fclaw2d_clawpatch_metric_data(glob,fine_patch,&xp,&yp,&zp,
                                          &xd,&yd,&zd,&areafine);
        }
        else
        {
            areafine = NULL;
        }

        fclaw2d_clawpatch_vt()->fort_interpolate2fine(&mx,&my,&mbc,&meqn,qcoarse,qfine,
                                                      areacoarse, areafine, &igrid,
                                                      &fclaw_opt->manifold);

    }
}

static
void clawpatch_average2coarse(fclaw2d_global_t *glob,
                              fclaw2d_patch_t *fine_patches,
                              fclaw2d_patch_t *coarse_patch,
                              int blockno, int fine0_patchno,
                              int coarse_patchno)

{
    const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
    
    int mx,my, mbc, meqn;
    double *qcoarse, *qfine;
    double *areacoarse, *areafine;
    double *xp,*yp,*zp,*xd,*yd,*zd;
    int igrid;
    fclaw2d_patch_t *fine_patch;

    mx = clawpatch_opt->mx;
    my = clawpatch_opt->my;
    mbc = clawpatch_opt->mbc;

    fclaw2d_clawpatch_metric_data(glob,coarse_patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&areacoarse);
    fclaw2d_clawpatch_soln_data(glob,coarse_patch,&qcoarse,&meqn);

    for(igrid = 0; igrid < 4; igrid++)
    {
        fine_patch = &fine_patches[igrid];

        fclaw2d_clawpatch_soln_data(glob,fine_patch,&qfine,&meqn);

        if (fclaw_opt->manifold)
        {
            fclaw2d_clawpatch_metric_data(glob,fine_patch,&xp,&yp,&zp,
                                          &xd,&yd,&zd,&areafine);
        }
        else
        {
            areafine = NULL;
        }

        fclaw2d_clawpatch_vt()->fort_average2coarse(&mx,&my,&mbc,&meqn,qcoarse,qfine,
                                                    areacoarse, areafine, &igrid,
                                                    &fclaw_opt->manifold);

    }
}

/* ------------------------------ Parallel ghost patches ------------------------------ */

static
void clawpatch_ghost_comm(fclaw2d_global_t* glob,
                          fclaw2d_patch_t* this_patch,
                          void *unpack_from_here, int time_interp,
                          int packmode)
{
    double *qpack = (double*) unpack_from_here;
    int meqn;
    double *qthis;
    double *area;
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    int ierror;

    int packextra = fclaw_opt->ghost_patch_pack_numextrafields;
    int packarea = packmode/2;   // (0,1)/2 = 0;  (2,3)/2 = 1;

    fclaw2d_clawpatch_timesync_data(glob,this_patch,time_interp,&qthis,&meqn);
    area = fclaw2d_clawpatch_get_area(glob,this_patch);

    int mx = clawpatch_opt->mx;
    int my = clawpatch_opt->my;
    int mbc = clawpatch_opt->mbc;
    int refratio = fclaw_opt->refratio;

    int mint = mbc*refratio;   /* # interior cells needed for averaging */
    int nghost = mbc;          /* # ghost values needed for interpolation */

    /* This is computed twice - here, and in fclaw2d_clawpatch_ghost_packsize */
    int wg = (2*nghost + mx)*(2*nghost + my);
    int hole = (mx - 2*mint)*(my - 2*mint);  /* Hole in center */
    FCLAW_ASSERT(hole >= 0);

    int psize = (wg - hole)*(meqn + packarea + packextra);
    FCLAW_ASSERT(psize > 0);

    int qareasize = (wg - hole)*(meqn + packarea);
    clawpatch_vt()->fort_local_ghost_pack(&mx,&my,&mbc,&meqn,&mint,qthis,area,
                                         qpack,&qareasize,&packmode,&ierror);
    FCLAW_ASSERT(ierror == 0);
    if (packextra)
    {
        qpack += qareasize;
        int extrasize = psize - qareasize;
        FCLAW_ASSERT(extrasize > 0);
        FCLAW_ASSERT(clawpatch_vt()->fort_local_ghost_pack_aux != NULL);
        clawpatch_vt()->fort_local_ghost_pack_aux(glob,this_patch,mint,
                                                  qpack,extrasize,
                                                  packmode,&ierror);
        FCLAW_ASSERT(ierror == 0);
    }

    if (ierror > 0)
    {
        fclaw_global_essentialf("clawpatch_ghost_comm  : ierror = %d\n",ierror);
        exit(0);
    }
}

static
size_t clawpatch_ghost_pack_elems(fclaw2d_global_t* glob)
{
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    const fclaw2d_clawpatch_options_t *clawpatch_opt = 
                         fclaw2d_clawpatch_get_options(glob);

    int mx = clawpatch_opt->mx;
    int my = clawpatch_opt->my;
    int mbc = clawpatch_opt->mbc;
    int meqn = clawpatch_opt->meqn;
    int refratio = fclaw_opt->refratio;

    int packarea = fclaw_opt->ghost_patch_pack_area && fclaw_opt->manifold;
    int packextra = fclaw_opt->ghost_patch_pack_numextrafields;

    int mint = refratio*mbc;
    int nghost = mbc;

    int wg = (2*nghost + mx)*(2*nghost + my);  /* Whole grid     */
    int hole = (mx - 2*mint)*(my - 2*mint);    /* Hole in center */
    FCLAW_ASSERT(hole >= 0);

    size_t psize = (wg - hole)*(meqn + packarea + packextra);
    FCLAW_ASSERT(psize > 0);

    return psize;
}    

size_t clawpatch_ghost_packsize(fclaw2d_global_t* glob)
{
    size_t esize = clawpatch_ghost_pack_elems(glob);
    return esize*sizeof(double);
}

static
void clawpatch_local_ghost_pack(fclaw2d_global_t *glob,
                                fclaw2d_patch_t *this_patch,
                                void *patch_data,
                                int time_interp)
{
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    int packarea = fclaw_opt->ghost_patch_pack_area && fclaw_opt->manifold;
    int packmode = 2*packarea;  // 0 or 2  (for pack)

    clawpatch_ghost_comm(glob,this_patch,patch_data, time_interp,packmode);
}

static
void clawpatch_remote_ghost_unpack(fclaw2d_global_t* glob,
                                   fclaw2d_patch_t* this_patch,
                                   int this_block_idx,
                                   int this_patch_idx,
                                   void *qdata, int time_interp)
{
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    int packarea = fclaw_opt->ghost_patch_pack_area && fclaw_opt->manifold;
    int packmode = 2*packarea + 1;  // 1 or 3  (for unpack)

    clawpatch_ghost_comm(glob,this_patch,qdata,time_interp,packmode);
}

static
void clawpatch_remote_ghost_build(fclaw2d_global_t *glob,
                                  fclaw2d_patch_t *this_patch,
                                  int blockno,
                                  int patchno,
                                  void *user)
{
    fclaw2d_build_mode_t build_mode =  *((fclaw2d_build_mode_t*) user);
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

    clawpatch_define(glob,this_patch,blockno,patchno,build_mode);

    if (fclaw_opt->manifold)
    {
        if (build_mode != FCLAW2D_BUILD_FOR_GHOST_AREA_PACKED)
        {
            fclaw2d_vt()->metric_compute_area(glob,this_patch,blockno,patchno);
        }
    }
}

static
void clawpatch_remote_ghost_delete(void *patchcp)
{
    FCLAW_ASSERT(patchcp != NULL);
    fclaw2d_clawpatch_t* cp = (fclaw2d_clawpatch_t*) patchcp;
    delete cp;
    patchcp = NULL;
}

/* ---------------------------- Parallel partitioning --------------------------------- */

static
size_t clawpatch_partition_packsize(fclaw2d_global_t* glob)
{
    const fclaw2d_clawpatch_options_t *clawpatch_opt 
                              = fclaw2d_clawpatch_get_options(glob);
    int mx = clawpatch_opt->mx;
    int my = clawpatch_opt->my;
    int mbc = clawpatch_opt->mbc;
    int meqn = clawpatch_opt->meqn;
    size_t psize = (2*mbc + mx)*(2*mbc + my)*meqn;  /* Store area */
    return psize*sizeof(double);
}

static
void clawpatch_partition_pack(fclaw2d_global_t *glob,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx,
                              void *pack_data_here)
    {
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    FCLAW_ASSERT(cp != NULL);

    cp->griddata.copyToMemory((double*) pack_data_here);
}

static
void clawpatch_partition_unpack(fclaw2d_global_t *glob,  
                                fclaw2d_domain_t *new_domain,
                                fclaw2d_patch_t *this_patch,
                                int this_block_idx,
                                int this_patch_idx,
                                void *unpack_data_from_here)
{
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);

    /* Time interp is false, since we only partition when all grids
       are time synchronized */
    cp->griddata.copyFromMemory((double*)unpack_data_from_here);
}

/* ------------------------------------ Virtual table  -------------------------------- */

static
fclaw2d_clawpatch_vtable_t* clawpatch_vt_init()
{
    FCLAW_ASSERT(s_clawpatch_vt.is_set == 0);
    return &s_clawpatch_vt;
}

void fclaw2d_clawpatch_vtable_initialize()
{
    fclaw2d_patch_vtable_initialize();
    fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();

    fclaw2d_clawpatch_vtable_t *clawpatch_vt = clawpatch_vt_init();

    /* These may be redefined by the user */
    /* Problem setup */
    patch_vt->patch_new             = clawpatch_new;
    patch_vt->patch_delete          = clawpatch_delete;
    patch_vt->build                 = clawpatch_build;
    patch_vt->build_from_fine       = clawpatch_build_from_fine;

    /* Time stepping */
    patch_vt->restore_step          = clawpatch_restore_step;
    patch_vt->save_step             = clawpatch_save_step;
    patch_vt->setup_timeinterp      = clawpatch_setup_timeinterp;

    /* Ghost filling - solver specific */
    patch_vt->copy_face            = clawpatch_copy_face;
    patch_vt->average_face         = clawpatch_average_face;
    patch_vt->interpolate_face     = clawpatch_interpolate_face;

    patch_vt->copy_corner          = clawpatch_copy_corner;
    patch_vt->average_corner       = clawpatch_average_corner;
    patch_vt->interpolate_corner   = clawpatch_interpolate_corner;

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
    patch_vt->partition_packsize    = clawpatch_partition_packsize;
    patch_vt->partition_pack        = clawpatch_partition_pack;
    patch_vt->partition_unpack      = clawpatch_partition_unpack;

    /* output functions */
    clawpatch_vt->cb_output_ascii   = cb_clawpatch_output_ascii; 
    clawpatch_vt->fort_header_ascii = NULL;    /* Defined in clawpack solvers */ 
    clawpatch_vt->fort_output_ascii = NULL;    /* Defined in clawpack solvers */

    fclaw2d_clawpatch_diagnostics_vtable_initialize();

    clawpatch_vt->is_set = 1;
}


/* --------------------------------- Misc access functions ---------------------------- */

/* These functions are not virtualized and are not defined by the 
   patch interface */

fclaw2d_clawpatch_vtable_t* fclaw2d_clawpatch_vt()
{
    FCLAW_ASSERT(s_clawpatch_vt.is_set != 0);
    return &s_clawpatch_vt;
}

/* Called from clawpack 4.6 and 5.0 */
void fclaw2d_clawpatch_save_current_step(fclaw2d_global_t* glob,
                                         fclaw2d_patch_t* this_patch)
{
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    cp->griddata_last = cp->griddata;
}


fclaw2d_clawpatch_t* fclaw2d_clawpatch_get_cp(fclaw2d_patch_t* this_patch)

{
    return clawpatch_data(this_patch);
}

void fclaw2d_clawpatch_grid_data(fclaw2d_global_t* glob,
                                 fclaw2d_patch_t* this_patch,
                                 int* mx, int* my, int* mbc,
                                 double* xlower, double* ylower,
                                 double* dx, double* dy)
{
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    *mx = cp->mx;
    *my = cp->my;
    *mbc = cp->mbc;
    *xlower = cp->xlower;
    *ylower = cp->ylower;
    *dx = cp->dx;
    *dy = cp->dy;
}

void fclaw2d_clawpatch_aux_data(fclaw2d_global_t *glob,
                                fclaw2d_patch_t *this_patch,
                                double **aux, int* maux)
{
    fclaw2d_clawpatch_t *cp = fclaw2d_clawpatch_get_cp (this_patch);

    *maux = cp->maux;
    *aux = cp->aux.dataPtr();
}

void fclaw2d_clawpatch_soln_data(fclaw2d_global_t* glob,
                                 fclaw2d_patch_t* this_patch,
                                 double **q, int* meqn)
{
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    *q = cp->griddata.dataPtr();
    *meqn = cp->meqn;
}

double *fclaw2d_clawpatch_get_q(fclaw2d_global_t* glob,
                                fclaw2d_patch_t* this_patch)
{
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    return cp->griddata.dataPtr();
}

double* fclaw2d_clawpatch_get_error(fclaw2d_global_t* glob,
                                    fclaw2d_patch_t* this_patch)
{
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    return cp->griderror.dataPtr();
}

size_t fclaw2d_clawpatch_size(fclaw2d_global_t *glob)
{
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
    int mx = clawpatch_opt->mx;
    int my = clawpatch_opt->my;
    int meqn = clawpatch_opt->meqn;
    int mbc = clawpatch_opt->mbc;
    size_t size = (mx+2*mbc)*(my+2*mbc)*meqn;

    return size;
}

void fclaw2d_clawpatch_timesync_data(fclaw2d_global_t* glob,
                                     fclaw2d_patch_t* this_patch,
                                     int time_interp,
                                     double **q, int* meqn)
{
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    *q = q_time_sync(cp, time_interp);
    *meqn = cp->meqn;
}


/* Return a pointer to either time interpolated data or regular grid data */
static double* q_time_sync(fclaw2d_clawpatch_t* cp, int time_interp)
{
    if (time_interp)
        return cp->griddata_time_interpolated.dataPtr();
    else
        return cp->griddata.dataPtr();
}


double* fclaw2d_clawpatch_get_q_timesync(fclaw2d_global_t* glob,
                                         fclaw2d_patch_t* this_patch,
                                         int time_interp)
{
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    return q_time_sync(cp, time_interp);
}

void fclaw2d_clawpatch_metric_data(fclaw2d_global_t* glob,
                                   fclaw2d_patch_t* this_patch,
                                   double **xp, double **yp, double **zp,
                                   double **xd, double **yd, double **zd,
                                   double **area)
{
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    *xp = cp->xp.dataPtr();
    *yp = cp->yp.dataPtr();
    *zp = cp->zp.dataPtr();
    *xd = cp->xd.dataPtr();
    *yd = cp->yd.dataPtr();
    *zd = cp->zd.dataPtr();
    *area = cp->area.dataPtr();
}

void fclaw2d_clawpatch_metric_data2(fclaw2d_global_t* glob,
                                    fclaw2d_patch_t* this_patch,
                                    double **xnormals, double **ynormals,
                                    double **xtangents, double **ytangents,
                                    double **surfnormals,
                                    double **edgelengths, double **curvature)
{
    /* or just call the member functions? */
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    *xnormals    = cp->xface_normals.dataPtr();
    *ynormals    = cp->yface_normals.dataPtr();
    *xtangents   = cp->xface_tangents.dataPtr();
    *ytangents   = cp->yface_tangents.dataPtr();
    *surfnormals = cp->surf_normals.dataPtr();
    *curvature   = cp->curvature.dataPtr();
    *edgelengths = cp->edge_lengths.dataPtr();
}

double* fclaw2d_clawpatch_get_area(fclaw2d_global_t* glob,
                                   fclaw2d_patch_t* this_patch)
{
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    return cp->area.dataPtr();
}


/* --------------------------------- Pillow grid  ------------------------------------- */

void fclaw2d_clawpatch_t::mb_exchange_block_corner_ghost(const int& a_corner,
                                                         fclaw2d_clawpatch_t *cp_corner,
                                                         int time_interp)
{
    double *qthis = q_time_sync(this,time_interp);
    double *qcorner = cp_corner->griddata.dataPtr();


    FCLAW2D_FORT_MB_EXCHANGE_BLOCK_CORNER_GHOST(mx, my, mbc, meqn, qthis, qcorner,
                                    a_corner, blockno);

}

// internal corners only a block boundaries.
void fclaw2d_clawpatch_t::mb_average_block_corner_ghost(const int& a_coarse_corner,
                                                        const int& a_refratio,
                                                        fclaw2d_clawpatch_t *cp_corner,
                                                        int a_time_interp)
{
    // 'this' is the finer grid; 'cp_corner' is the coarser grid.
    double *qcoarse = q_time_sync(this,a_time_interp);


    double *areacoarse = this->area.dataPtr();
    double *areafine = cp_corner->area.dataPtr();
    double *qfine = cp_corner->griddata.dataPtr();

    FCLAW2D_FORT_MB_AVERAGE_BLOCK_CORNER_GHOST(mx,my,mbc,meqn,
                                               a_refratio,qcoarse,qfine,
                                               areacoarse,areafine,
                                               a_coarse_corner,blockno);
}


void fclaw2d_clawpatch_t::mb_interpolate_block_corner_ghost(const int& a_coarse_corner,
                                                            const int& a_refratio,
                                                            fclaw2d_clawpatch_t *cp_corner,
                                                            int a_time_interp)

{
    double *qcoarse = q_time_sync(this, a_time_interp);

    /* qcorner is the finer level. */
    double *qfine = cp_corner->griddata.dataPtr();

    FCLAW2D_FORT_MB_INTERPOLATE_BLOCK_CORNER_GHOST(mx, my, mbc, meqn,
                                                   a_refratio, qcoarse, qfine,
                                                   a_coarse_corner, blockno);
}

