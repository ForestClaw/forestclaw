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

#include <fclaw2d_forestclaw.h>
#include <fclaw_clawpatch3.hpp>

#include <fclaw2d_timeinterp.h>
#include <fclaw2d_metric.h>
#include <fclaw2d_neighbors_fort.h>


static
fclaw_clawpatch3_vtable_t* clawpatch3_vt()
{
    static fclaw_clawpatch3_vtable_t s_clawpatch3_vt;
    return &s_clawpatch3_vt;
}

fclaw_clawpatch3_vtable_t* fclaw_clawpatch3_vt()
{
    return clawpatch3_vt();
}


#if 0
static
void setup_area_storage(fclaw_clawpatch3_t* cp);

static
void setup_metric_storage(fclaw_clawpatch3_t* cp);
#endif

static
void metric_setup(fclaw2d_global_t* glob,
                  fclaw2d_patch_t* this_patch,
                  int blockno,
                  int patchno);

static
void ghost_comm(fclaw2d_global_t* glob,
                fclaw2d_patch_t* this_patch,
                double *qpack, int time_interp,
                int packmode);

static
fclaw_clawpatch3_t* clawpatch3_data(fclaw2d_patch_t *this_patch);

static
double* q_time_sync(fclaw_clawpatch3_t* cp, int time_interp);


/* ------------------------------------------------------------
   Solution access functions
   ---------------------------------------------------------- */
void* fclaw_clawpatch3_new_patch()
{
    fclaw_clawpatch3_t *cp = new fclaw_clawpatch3_t;
    FCLAW_ASSERT(cp != NULL);
    return (void*) cp;
}

void fclaw_clawpatch3_delete_patch(void *patchcp)
{
    FCLAW_ASSERT(patchcp != NULL);
    fclaw_clawpatch3_t* cp = (fclaw_clawpatch3_t*) patchcp;
    delete cp;
    patchcp = NULL;
}

/* Access clawpatch stored in patch->user_data */
static
fclaw_clawpatch3_t* clawpatch3_data(fclaw2d_patch_t *this_patch)
{
    fclaw_clawpatch3_t *cp = (fclaw_clawpatch3_t*) fclaw2d_patch_get_user_patch(this_patch);
    return cp;
}


fclaw_clawpatch3_t* fclaw_clawpatch3_get_cp(fclaw2d_patch_t* this_patch)

{
    return clawpatch3_data(this_patch);
}

void fclaw_clawpatch3_grid_data(fclaw2d_global_t* glob,
                                  fclaw2d_patch_t* this_patch,
                                  int* mx, int* my, int* mz, int* mbc,
                                  double* xlower, double* ylower, double* zlower,
                                  double* dx, double* dy, double* dz)
{
    fclaw_clawpatch3_t *cp = clawpatch3_data(this_patch);
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

void fclaw_clawpatch3_aux_data(fclaw2d_global_t *glob,
                                fclaw2d_patch_t *this_patch,
                                double **aux, int* maux)
{
    fclaw_clawpatch3_t *cp = fclaw_clawpatch3_get_cp (this_patch);

    *maux = cp->maux;
    *aux = cp->aux.dataPtr();
}

void fclaw_clawpatch3_soln_data(fclaw2d_global_t* glob,
                                 fclaw2d_patch_t* this_patch,
                                 double **q, int* meqn)
{
    fclaw_clawpatch3_t *cp = clawpatch3_data(this_patch);
    *q = cp->griddata.dataPtr();
    *meqn = cp->meqn;
}

double *fclaw_clawpatch3_get_q(fclaw2d_global_t* glob,
                                 fclaw2d_patch_t* this_patch)
{
    fclaw_clawpatch3_t *cp = clawpatch3_data(this_patch);
    return cp->griddata.dataPtr();
}

double* fclaw_clawpatch3_get_error(fclaw2d_global_t* glob,
                                     fclaw2d_patch_t* this_patch)
{
    fclaw_clawpatch3_t *cp = clawpatch3_data(this_patch);
    return cp->griderror.dataPtr();
}

size_t fclaw_clawpatch3_size(fclaw2d_global_t *glob)
{
    const fclaw_clawpatch3_options_t *clawpatch3_opt;
    clawpatch3_opt = fclaw_clawpatch3_get_options(glob);
    
    int mx = clawpatch3_opt->mx;
    int my = clawpatch3_opt->my;
    int mz = clawpatch3_opt->mz;
    int meqn = clawpatch3_opt->meqn;
    int mbc = clawpatch3_opt->mbc;
    size_t size = (mx+2*mbc)*(my+2*mbc)*(mz+2*mbc)*meqn;

    return size;
}

/* ------------------------------------------------
   Functions for handling time-interpolated data
   ------------------------------------------------ */
void fclaw_clawpatch3_setup_timeinterp(fclaw2d_global_t *glob,
                                         fclaw2d_patch_t *this_patch,
                                         double alpha)
{
    /* We use the pack size here to make sure we are setting
       everything correctly;  it isn't needed for memory
       allocation */
    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(glob);
    const fclaw_clawpatch3_options_t *clawpatch3_opt = fclaw_clawpatch3_get_options(glob);

    int mx = clawpatch3_opt->mx;
    int my = clawpatch3_opt->my;
    int mz = clawpatch3_opt->mz;

    int meqn = clawpatch3_opt->meqn;
    int mbc = clawpatch3_opt->mbc;
    int mint = gparms->interp_stencil_width/2+1;  /* Assume interp stencils have odd width */

    /* 2.5 modification */
    int hole = (mx - 2*mint)*(my - 2*mint)*mz;  /* Hole in center */
    FCLAW_ASSERT(hole >= 0);

    int wg = mx*my*mz;  /* Whole grid but no ghost cells.  Ghost cells will be averaged from finer
                      * level. */
    int psize = (wg - hole)*meqn;
    FCLAW_ASSERT(psize > 0);

    /* Store time interpolated data that will be use in coarse grid
       exchanges */
    fclaw_clawpatch3_t *cp = clawpatch3_data(this_patch);
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

    clawpatch3_vt()->fort_timeinterp(&mx,&my,&mz,&mbc,&meqn,&psize,
                                     qcurr,qlast,qinterp,&alpha,&ierror); 

}

void fclaw_clawpatch3_timesync_data(fclaw2d_global_t* glob,
                                      fclaw2d_patch_t* this_patch,
                                      fclaw_bool time_interp,
                                      double **q, int* meqn)
{
    fclaw_clawpatch3_t *cp = clawpatch3_data(this_patch);
    *q = q_time_sync(cp, time_interp);
    *meqn = cp->meqn;
}

double* fclaw_clawpatch3_get_q_timesync(fclaw2d_global_t* glob,
                                          fclaw2d_patch_t* this_patch,
                                          int time_interp)
{
    fclaw_clawpatch3_t *cp = clawpatch3_data(this_patch);
    return q_time_sync(cp, time_interp);
}

/* This is called from libraries routines (clawpack4.6, clawpack5, etc) */
void fclaw_clawpatch3_save_current_step(fclaw2d_global_t* glob,
                                          fclaw2d_patch_t* this_patch)
{
    fclaw_clawpatch3_t *cp = clawpatch3_data(this_patch);
    cp->griddata_last = cp->griddata;
}

/* This is called from libraries routines (clawpack4.6, clawpack5, etc) */
void fclaw_clawpatch3_restore_step(fclaw2d_global_t* glob,
                                     fclaw2d_patch_t* this_patch)
{
    fclaw_clawpatch3_t *cp = clawpatch3_data(this_patch);
    cp->griddata = cp->griddata_save;
}

/* This is called from libraries routines (clawpack4.6, clawpack5, etc) */
void fclaw_clawpatch3_save_step(fclaw2d_global_t* glob,
                                  fclaw2d_patch_t* this_patch)
{
    fclaw_clawpatch3_t *cp = clawpatch3_data(this_patch);
    cp->griddata_save = cp->griddata;
}


/* ------------------------------------------------------------------
   Manifold setup and access
 ------------------------------------------------------------------ */

static
void metric_setup(fclaw2d_global_t* glob,
                  fclaw2d_patch_t* this_patch,
                  int blockno,
                  int patchno)
{
    fclaw_global_essentialf("metric_setup (fclaw_clawpatch3.cpp): This function should not be called\n");
    exit(0);
#if 0
    /* vt.patch_manifold_setup_mesh(...) */
    fclaw2d_metric_setup_mesh(glob,this_patch,blockno,patchno);

    /* vt.patch_manifold_compute_normals(...) */
    fclaw2d_metric_compute_normals(glob,this_patch,blockno,patchno);
#endif
}

#if 0
void fclaw_clawpatch3_metric_data(fclaw2d_global_t* glob,
                                    fclaw2d_patch_t* this_patch,
                                    double **xp, double **yp, double **zp,
                                    double **xd, double **yd, double **zd,
                                    double **area)
{
    fclaw_clawpatch3_t *cp = clawpatch3_data(this_patch);
    *xp = cp->xp.dataPtr();
    *yp = cp->yp.dataPtr();
    *zp = cp->zp.dataPtr();
    *xd = cp->xd.dataPtr();
    *yd = cp->yd.dataPtr();
    *zd = cp->zd.dataPtr();
    *area = cp->area.dataPtr();
}

void fclaw_clawpatch3_metric_data2(fclaw2d_global_t* glob,
                                     fclaw2d_patch_t* this_patch,
                                     double **xnormals, double **ynormals,
                                     double **xtangents, double **ytangents,
                                     double **surfnormals,
                                     double **edgelengths, double **curvature)
{
    /* or just call the member functions? */
    fclaw_clawpatch3_t *cp = clawpatch3_data(this_patch);
    *xnormals    = cp->xface_normals.dataPtr();
    *ynormals    = cp->yface_normals.dataPtr();
    *xtangents   = cp->xface_tangents.dataPtr();
    *ytangents   = cp->yface_tangents.dataPtr();
    *surfnormals = cp->surf_normals.dataPtr();
    *curvature   = cp->curvature.dataPtr();
    *edgelengths = cp->edge_lengths.dataPtr();
}
#endif
double* fclaw_clawpatch3_get_area(fclaw2d_global_t* glob,
                                    fclaw2d_patch_t* this_patch)
{
    fclaw_clawpatch3_t *cp = clawpatch3_data(this_patch);
    return cp->area.dataPtr();
}

/* ------------------------------------------------------------------
   Re-partition and rebuild new domains, or construct initial domain
 -------------------------------------------------------------------- */

void fclaw_clawpatch3_define(fclaw2d_global_t* glob,
                               fclaw2d_patch_t *this_patch,
                               int blockno, int patchno,
                               fclaw2d_build_mode_t build_mode)
{
    /* We are getting closer to getting rid the class fclaw_clawpatch3_t */
    fclaw_clawpatch3_t *cp = clawpatch3_data(this_patch);

    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(glob);
    const fclaw_clawpatch3_options_t *clawpatch3_opt = fclaw_clawpatch3_get_options(glob);

    cp->mx = clawpatch3_opt->mx;
    cp->my = clawpatch3_opt->my;
    cp->mz = clawpatch3_opt->mz;
    cp->mbc = clawpatch3_opt->mbc;
    cp->blockno = blockno;
    cp->meqn = clawpatch3_opt->meqn;
    cp->maux = clawpatch3_opt->maux;

    for (int icorner=0; icorner < 4; icorner++)
    {
        fclaw2d_patch_set_block_corner_count(glob,this_patch,icorner,0);
    }

    fclaw2d_map_context_t* cont = glob->cont;

    int is_brick = FCLAW2D_MAP_IS_BRICK(&cont);

    cp->manifold = gparms->manifold;

    if (cp->manifold)
    {
        cp->xlower = this_patch->xlower; 
        cp->ylower = this_patch->ylower;
        cp->xupper = this_patch->xupper;
        cp->yupper = this_patch->yupper;
    }
    else
    {
        double ax = gparms->ax;
        double bx = gparms->bx;
        double ay = gparms->ay;
        double by = gparms->by;
        double az = gparms->az;
        double bz = gparms->bz;

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
        cp->zlower = az;
        cp->zupper = bz;
    }

    cp->dx = (cp->xupper - cp->xlower)/cp->mx;
    cp->dy = (cp->yupper - cp->ylower)/cp->my;
    cp->dz = (cp->zupper - cp->zlower)/cp->mz;

    int ll[3];
    int ur[3];
    for (int idir = 0; idir < 3; idir++)
    {
        ll[idir] = 1-cp->mbc;
    }
    ur[0] = cp->mx + cp->mbc;
    ur[1] = cp->my + cp->mbc;
    ur[2] = cp->mz + cp->mbc;
    Box box(ll,ur);

    // This will destroy any existing memory n griddata.
    cp->griddata.define(box, cp->meqn);
    if (gparms->subcycle)
    {
        cp->griddata_time_interpolated.define(box, cp->meqn);
    }
    if (gparms->compute_error)
    {
        cp->griderror.define(box,cp->meqn);
    }

    if (clawpatch3_opt->maux > 0)
    {
      cp->aux.define(box,cp->maux);
    }
    // Set up storage for metric terms, if needed.
#if 0
    if (gparms->manifold)
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
#endif  
    if (build_mode != FCLAW2D_BUILD_FOR_UPDATE)
    {
        return;
    }

    cp->griddata_last.define(box, cp->meqn);
    cp->griddata_save.define(box, cp->meqn);

}


/* This is called from fclaw2d_regrid_new_domain_setup (in fclaw2d_regrid.cpp)
   every time a new domain is created after regridding.  This is followed by
   an "iterate_adapted" step, in which data in these patches is copied, averaged or
   interpolated. */

void fclaw_clawpatch3_build(fclaw2d_global_t *glob,
                             fclaw2d_patch_t *this_patch,
                             int blockno,
                             int patchno,
                             void *user)
{
    fclaw2d_build_mode_t build_mode =  *((fclaw2d_build_mode_t*) user);
    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(glob);

    fclaw_clawpatch3_define(glob,this_patch,blockno,patchno,build_mode);

    if (gparms->manifold)
    {
        fclaw2d_vt()->metric_compute_area(glob,this_patch,blockno,patchno);
        metric_setup(glob,this_patch,blockno,patchno);
    }
}

void fclaw_clawpatch3_build_from_fine(fclaw2d_global_t *glob,
                                       fclaw2d_patch_t *fine_patches,
                                       fclaw2d_patch_t *coarse_patch,
                                       int blockno,
                                       int coarse_patchno,
                                       int fine0_patchno,
                                       fclaw2d_build_mode_t build_mode)
{
    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(glob);

    fclaw_clawpatch3_define(glob,coarse_patch,blockno,coarse_patchno,build_mode);

    if (gparms->manifold)
    {
        /* Don't recompute the area, but rather average from finer areas */
        fclaw2d_metric_average_area(glob,fine_patches,coarse_patch,
                                    blockno, coarse_patchno, fine0_patchno);

        metric_setup(glob,coarse_patch,blockno,coarse_patchno);
    }
}

void fclaw_clawpatch3_build_ghost(fclaw2d_global_t *glob,
                                   fclaw2d_patch_t *this_patch,
                                   int blockno,
                                   int patchno,
                                   void *user)
{
    fclaw2d_build_mode_t build_mode =  *((fclaw2d_build_mode_t*) user);
    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(glob);

    fclaw_clawpatch3_define(glob,this_patch,blockno,patchno,build_mode);

    if (gparms->manifold)
    {
        if (build_mode != FCLAW2D_BUILD_FOR_GHOST_AREA_PACKED)
        {
            fclaw2d_vt()->metric_compute_area(glob,this_patch,blockno,patchno);
        }
    }
}

/* --------------------------------------------------------------
   Parallel ghost exchanges.

   Note this is different from the packing/unpacking that happens
   when partitioning the domain
   --------------------------------------------------------------*/

size_t fclaw_clawpatch3_ghost_packsize(fclaw2d_global_t* glob)
{
    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(glob);
    const fclaw_clawpatch3_options_t *clawpatch3_opt = fclaw_clawpatch3_get_options(glob);

    int mx = clawpatch3_opt->mx;
    int my = clawpatch3_opt->my;
    int mz = clawpatch3_opt->mz;
    int mbc = clawpatch3_opt->mbc;
    int meqn = clawpatch3_opt->meqn;
    int refratio = gparms->refratio;

    int mint = refratio*mbc;

    int wg = (2*mbc + mx)*(2*mbc + my)*(2*mbc + mz);  /* Whole grid     */
    /* if patch_dim = 3, refine_dim = 2, then hole = ... * mz
       otherwise, hole = ... * (mz - 2*mint) */
    int hole = (mx - 2*mint)*(my - 2*mint)*mz;    /* Hole in center */
    FCLAW_ASSERT(hole >= 0);

    int packarea = gparms->ghost_patch_pack_area && gparms->manifold;
    int packextra = gparms->ghost_patch_pack_numextrafields;
    int nfields = meqn + packarea + packextra;
    size_t psize = (wg - hole)*nfields;
    FCLAW_ASSERT(psize > 0);

    return psize*sizeof(double);
}

void fclaw_clawpatch3_local_ghost_alloc(fclaw2d_global_t* glob,
                                         fclaw2d_patch_t* this_patch,
                                         void** q)
{
    /* Create contiguous block for data and area */
    int msize = fclaw_clawpatch3_ghost_packsize(glob);
    *q = (void*) FCLAW_ALLOC(double,msize);
    FCLAW_ASSERT(*q != NULL);

}

void fclaw_clawpatch3_local_ghost_free(fclaw2d_global_t* glob,
                                        void **q)
{
    FCLAW_FREE(*q);
    *q = NULL;
}


static
void ghost_comm(fclaw2d_global_t* glob,
                fclaw2d_patch_t* this_patch,
                double *qpack, int time_interp,
                int packmode)
{
    int meqn;
    double *qthis;
    double *area;
    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(glob);
    const fclaw_clawpatch3_options_t *clawpatch3_opt = fclaw_clawpatch3_get_options(glob);

    int ierror;

    int packarea = packmode/2;   // (0,1)/2 = 0;  (2,3)/2 = 1;

    fclaw_clawpatch3_timesync_data(glob,this_patch,time_interp,&qthis,&meqn);
    area = fclaw_clawpatch3_get_area(glob,this_patch);

    int mx = clawpatch3_opt->mx;
    int my = clawpatch3_opt->my;
    int mz = clawpatch3_opt->mz;
    int mbc = clawpatch3_opt->mbc;
    int refratio = gparms->refratio;

    int mint = mbc*refratio;   /* # interior cells needed for averaging */
    int nghost = mbc;          /* # ghost values needed for interpolation */

    /* This is computed twice - here, and in fclaw_clawpatch3_ghost_packsize */
    int wg = (2*nghost + mx)*(2*nghost + my)*(2*nghost + mz);
    /* 2.5 check needed */
    int hole = (mx - 2*mint)*(my - 2*mint)*mz;  /* Hole in center */
    FCLAW_ASSERT(hole >= 0);

    int psize = (wg - hole)*(meqn + packarea + gparms->ghost_patch_pack_numextrafields);
    FCLAW_ASSERT(psize > 0);

    int qareasize = (wg - hole)*(meqn + packarea);
    clawpatch3_vt()->fort_ghostpack_qarea(&mx,&my,&mz,&mbc,&meqn,&mint,qthis,area,
                                          qpack,&qareasize,&packmode,&ierror);
    FCLAW_ASSERT(ierror == 0);
    if (gparms->ghost_patch_pack_extra)
    {
      qpack += qareasize;
      int extrasize = psize - qareasize;
      FCLAW_ASSERT(clawpatch3_vt()->ghostpack_extra != NULL);
      clawpatch3_vt()->ghostpack_extra(glob,this_patch,mint,qpack,extrasize,packmode,&ierror);
      FCLAW_ASSERT(ierror == 0);
    }

    if (ierror > 0)
    {
        fclaw_global_essentialf("ghost_pack (fclaw2d_clawpatch.cpp) : ierror = %d\n",ierror);
        exit(0);
    }
}

void fclaw_clawpatch3_ghost_pack(fclaw2d_global_t *glob,
                                  fclaw2d_patch_t *this_patch,
                                  double *patch_data,
                                  int time_interp)
{
    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(glob);
    int packarea = gparms->ghost_patch_pack_area && gparms->manifold;
    int packmode = 2*packarea;  // 0 or 2  (for pack)

    ghost_comm(glob,this_patch,patch_data, time_interp,packmode);
}


void fclaw_clawpatch3_ghost_unpack(fclaw2d_global_t* glob,
                                    fclaw2d_patch_t* this_patch,
                                    int this_block_idx,
                                    int this_patch_idx,
                                    double *qdata, fclaw_bool time_interp)
{
    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(glob);
    int packarea = gparms->ghost_patch_pack_area && gparms->manifold;
    int packmode = 2*packarea + 1;  // 1 or 3  (for unpack)

    ghost_comm(glob,this_patch,qdata,time_interp,packmode);
}

/* --------------------------------------------------------
   Domain partitioning.

   Pack local patches on this processor before re-partitioning;
   retrieve patches that migrated across parallel borders to
   load balance the computation.
   -------------------------------------------------------*/

size_t fclaw_clawpatch3_partition_packsize(fclaw2d_global_t* glob)
{
    const fclaw_clawpatch3_options_t *clawpatch3_opt = fclaw_clawpatch3_get_options(glob);
    int mx = clawpatch3_opt->mx;
    int my = clawpatch3_opt->my;
    int mz = clawpatch3_opt->mz;
    int mbc = clawpatch3_opt->mbc;
    int meqn = clawpatch3_opt->meqn;
    size_t size = (2*mbc + mx)*(2*mbc + my)*(2*mbc + mz)*meqn;  /* Store area */
    return size*sizeof(double);
}

void fclaw_clawpatch3_partition_pack(fclaw2d_global_t *glob,
                                      fclaw2d_patch_t *this_patch,
                                      int this_block_idx,
                                      int this_patch_idx,
                                      void *user)
{
    fclaw2d_domain_t *domain = glob->domain;

    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    int patch_num = this_block->num_patches_before + this_patch_idx;
    double* patch_data = (double*) ((void**)user)[patch_num];

    fclaw_clawpatch3_t *cp = clawpatch3_data(this_patch);
    FCLAW_ASSERT(cp != NULL);

    cp->griddata.copyToMemory(patch_data);
}

void fclaw_clawpatch3_partition_unpack(fclaw2d_global_t *glob,
                                        fclaw2d_domain_t *new_domain,
                                        fclaw2d_patch_t *this_patch,
                                        int this_block_idx,
                                        int this_patch_idx,
                                        void *user)
{
    fclaw2d_block_t *this_block = &new_domain->blocks[this_block_idx];
    int patch_num = this_block->num_patches_before + this_patch_idx;
    double* patch_data = (double*) ((void**)user)[patch_num];

    fclaw_clawpatch3_t *cp = clawpatch3_data(this_patch);

    /* Time interp is false, since we only partition when all grids
       are time synchronized */
    cp->griddata.copyFromMemory(patch_data);
}

void fclaw_clawpatch3_copy_face(fclaw2d_global_t *glob,
                                  fclaw2d_patch_t *this_patch,
                                  fclaw2d_patch_t *neighbor_patch,
                                  int iface,
                                  int time_interp,
                                  fclaw2d_transform_data_t *transform_data)

{
    int meqn,mx,my,mz,mbc;
    double *qthis, *qneighbor;
    const fclaw_clawpatch3_options_t *clawpatch3_opt = fclaw_clawpatch3_get_options(glob);

    mx = clawpatch3_opt->mx;
    my = clawpatch3_opt->my;
    mz = clawpatch3_opt->mz;
    mbc = clawpatch3_opt->mbc;

    fclaw_clawpatch3_timesync_data(glob,this_patch,time_interp,&qthis,&meqn);
    fclaw_clawpatch3_timesync_data(glob,neighbor_patch,time_interp,&qneighbor,&meqn);

    clawpatch3_vt()->fort_copy_face(&mx,&my,&mz,&mbc,&meqn,qthis,qneighbor,&iface,&transform_data);
}

void fclaw_clawpatch3_average_face(fclaw2d_global_t *glob,
                                    fclaw2d_patch_t *coarse_patch,
                                    fclaw2d_patch_t *fine_patch,
                                    int idir,
                                    int iface_coarse,
                                    int p4est_refineFactor,
                                    int refratio,
                                    fclaw_bool time_interp,
                                    int igrid,
                                    fclaw2d_transform_data_t* transform_data)
{
    int meqn,mx,my,mz,mbc;
    double *qcoarse, *qfine;
    double *areacoarse, *areafine;

    const amr_options_t* gparms = fclaw2d_forestclaw_get_options(glob);
    const fclaw_clawpatch3_options_t *clawpatch3_opt = fclaw_clawpatch3_get_options(glob);

    fclaw_clawpatch3_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
    qfine = fclaw_clawpatch3_get_q(glob,fine_patch);

    /* These will be empty for non-manifolds cases */
    areacoarse = fclaw_clawpatch3_get_area(glob,coarse_patch);
    areafine = fclaw_clawpatch3_get_area(glob,fine_patch);

    mx = clawpatch3_opt->mx;
    my = clawpatch3_opt->my;
    mz = clawpatch3_opt->mz;
    mbc = clawpatch3_opt->mbc;

    int manifold = gparms->manifold;
    clawpatch3_vt()->fort_average_face(&mx,&my,&mz,&mbc,&meqn,qcoarse,qfine,areacoarse,areafine,
                                      &idir,&iface_coarse, &p4est_refineFactor, &refratio,
                                      &igrid,&manifold,&transform_data);


}

void fclaw_clawpatch3_interpolate_face(fclaw2d_global_t *glob,
                                        fclaw2d_patch_t *coarse_patch,
                                        fclaw2d_patch_t *fine_patch,
                                        int idir,
                                        int iside,
                                        int p4est_refineFactor,
                                        int refratio,
                                        fclaw_bool time_interp,
                                        int igrid,
                                        fclaw2d_transform_data_t* transform_data)
{
    int meqn,mx,my,mz,mbc;
    double *qcoarse, *qfine;
    const fclaw_clawpatch3_options_t *clawpatch3_opt = fclaw_clawpatch3_get_options(glob);

    fclaw_clawpatch3_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
    qfine = fclaw_clawpatch3_get_q(glob,fine_patch);

    mx = clawpatch3_opt->mx;
    my = clawpatch3_opt->my;
    mz = clawpatch3_opt->mz;
    mbc = clawpatch3_opt->mbc;

    clawpatch3_vt()->fort_interpolate_face(&mx,&my,&mz,&mbc,&meqn,qcoarse,qfine,&idir,&iside,
                                           &p4est_refineFactor,&refratio,&igrid,&transform_data);
}


void fclaw_clawpatch3_copy_corner(fclaw2d_global_t *glob,
                                   fclaw2d_patch_t *this_patch,
                                   fclaw2d_patch_t *corner_patch,
                                   int icorner,
                                   int time_interp,
                                   fclaw2d_transform_data_t *transform_data)
{
    int meqn,mx,my,mz,mbc;
    double *qthis, *qcorner;
    const fclaw_clawpatch3_options_t *clawpatch3_opt = fclaw_clawpatch3_get_options(glob);

    mx = clawpatch3_opt->mx;
    my = clawpatch3_opt->my;
    mz = clawpatch3_opt->mz;
    mbc = clawpatch3_opt->mbc;

    fclaw_clawpatch3_timesync_data(glob,this_patch,time_interp,&qthis,&meqn);
    fclaw_clawpatch3_timesync_data(glob,corner_patch,time_interp,&qcorner,&meqn);

    clawpatch3_vt()->fort_copy_corner(&mx,&my,&mz,&mbc,&meqn,qthis,qcorner,&icorner,&transform_data);

}

void fclaw_clawpatch3_average_corner(fclaw2d_global_t *glob,
                                      fclaw2d_patch_t *coarse_patch,
                                      fclaw2d_patch_t *fine_patch,
                                      int coarse_corner,
                                      int refratio,
                                      fclaw_bool time_interp,
                                      fclaw2d_transform_data_t* transform_data)
{
    int meqn,mx,my,mz,mbc;
    double *qcoarse, *qfine;
    double *areacoarse, *areafine;

    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(glob);
    const fclaw_clawpatch3_options_t *clawpatch3_opt = fclaw_clawpatch3_get_options(glob);

    fclaw_clawpatch3_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
    qfine = fclaw_clawpatch3_get_q(glob,fine_patch);

    areacoarse = fclaw_clawpatch3_get_area(glob,coarse_patch);
    areafine = fclaw_clawpatch3_get_area(glob,fine_patch);

    mx = clawpatch3_opt->mx;
    my = clawpatch3_opt->my;
    mz = clawpatch3_opt->mz;
    mbc = clawpatch3_opt->mbc;

    int manifold = gparms->manifold;
    clawpatch3_vt()->fort_average_corner(&mx,&my,&mz,&mbc,&meqn,&refratio,
                                        qcoarse,qfine,areacoarse,areafine,
                                        &manifold,&coarse_corner,&transform_data);
}


void fclaw_clawpatch3_interpolate_corner(fclaw2d_global_t* glob,
                                          fclaw2d_patch_t* coarse_patch,
                                          fclaw2d_patch_t* fine_patch,
                                          int coarse_corner,
                                          int refratio,
                                          fclaw_bool time_interp,
                                          fclaw2d_transform_data_t* transform_data)

{
    int meqn,mx,my,mz,mbc;
    double *qcoarse, *qfine;
    const fclaw_clawpatch3_options_t *clawpatch3_opt = fclaw_clawpatch3_get_options(glob);

    mx = clawpatch3_opt->mx;
    my = clawpatch3_opt->my;
    mz = clawpatch3_opt->mz;
    mbc = clawpatch3_opt->mbc;
    meqn = clawpatch3_opt->meqn;


    fclaw_clawpatch3_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
    qfine = fclaw_clawpatch3_get_q(glob,fine_patch);

    clawpatch3_vt()->fort_interpolate_corner(&mx,&my,&mz,&mbc,&meqn,
                                            &refratio,qcoarse,qfine,
                                            &coarse_corner,&transform_data);
}




/* ----------------------------------------------------------------
   fclaw_clawpatch3_t member functions (previous in fclaw_clawpatch3_t.cpp)
   ---------------------------------------------------------------- */

Box fclaw_clawpatch3_t::dataBox()
{
    return griddata.box();
}

Box fclaw_clawpatch3_t::areaBox()
{
    return area.box();
}
#if 0
Box fclaw_clawpatch3_t::edgeBox()
{
    return edge_lengths.box();
}

Box fclaw_clawpatch3_t::nodeBox()
{
    return xp.box();
}
#endif

/* Return a pointer to either time interpolated data or regular grid data */
static double* q_time_sync(fclaw_clawpatch3_t* cp, int time_interp)
{
    if (time_interp)
        return cp->griddata_time_interpolated.dataPtr();
    else
        return cp->griddata.dataPtr();
}

/* ----------------------------------------------------------------
   Special case : Pillow grid ghost exchanges/average/interpolate
   ---------------------------------------------------------------- */
#if 0
void fclaw_clawpatch3_t::mb_exchange_block_corner_ghost(const int& a_corner,
                                                         fclaw_clawpatch3_t *cp_corner,
                                                         int time_interp)
{
    double *qthis = q_time_sync(this,time_interp);
    double *qcorner = cp_corner->griddata.dataPtr();


    FCLAW2D_FORT_MB_EXCHANGE_BLOCK_CORNER_GHOST(mx, my, mbc, meqn, qthis, qcorner,
                                    a_corner, blockno);

}

// internal corners only a block boundaries.
void fclaw_clawpatch3_t::mb_average_block_corner_ghost(const int& a_coarse_corner,
                                                        const int& a_refratio,
                                                        fclaw_clawpatch3_t *cp_corner,
                                                        fclaw_bool a_time_interp)
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


void fclaw_clawpatch3_t::mb_interpolate_block_corner_ghost(const int& a_coarse_corner,
                                                            const int& a_refratio,
                                                            fclaw_clawpatch3_t *cp_corner,
                                                            fclaw_bool a_time_interp)

{
    double *qcoarse = q_time_sync(this, a_time_interp);

    /* qcorner is the finer level. */
    double *qfine = cp_corner->griddata.dataPtr();

    FCLAW2D_FORT_MB_INTERPOLATE_BLOCK_CORNER_GHOST(mx, my, mbc, meqn,
                                                   a_refratio, qcoarse, qfine,
                                                   a_coarse_corner, blockno);
}
/* ----------------------------------------------------------------
   Mapped grids
   ---------------------------------------------------------------- */

static void setup_area_storage(fclaw_clawpatch3_t* cp)
{
    int mx = cp->mx;
    int my = cp->my;
    int mbc = cp->mbc;

    int ll[SpaceDim];
    int ur[SpaceDim];
    for (int idir = 0; idir < SpaceDim; idir++)
    {
        ll[idir] = -mbc;
    }
    ur[0] = mx + mbc + 1;
    ur[1] = my + mbc + 1;

    Box box_p(ll,ur);
    cp->area.define(box_p,1);
}

static void setup_metric_storage(fclaw_clawpatch3_t* cp)
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

    int ll[SpaceDim];
    int ur[SpaceDim];
    for (int idir = 0; idir < SpaceDim; idir++)
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
    for (int idir = 0; idir < SpaceDim; idir++)
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
#endif

/* ----------------------------------------------------------------
   Set defaults for clawpatch virtual table
   ---------------------------------------------------------------- */


void fclaw_clawpatch3_init_vtable_defaults()
{
    fclaw2d_diagnostics_vtable_t *diag_vt = fclaw2d_diagnostics_vt();
    fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();

    /* These must be redefined by the solver and user */
    patch_vt->initialize         = NULL;
    patch_vt->physical_bc        = NULL;
    patch_vt->single_step_update = NULL;

    /* These may be redefined by the user */
    /* Problem setup */
    patch_vt->patch_new             = fclaw_clawpatch3_new_patch;
    patch_vt->patch_delete          = fclaw_clawpatch3_delete_patch;
    patch_vt->setup                 = NULL;
    patch_vt->setup_ghost           = NULL;
    patch_vt->build                 = &fclaw_clawpatch3_build;
    patch_vt->build_from_fine       = &fclaw_clawpatch3_build_from_fine;
    patch_vt->restore_step          = &fclaw_clawpatch3_restore_step;
    patch_vt->save_step             = &fclaw_clawpatch3_save_step;

    /* Ghost filling - solver specific */
    /* 2.5 check needed, may need fclaw3d_clawpatch3_xxxx */
    patch_vt->copy_face            = fclaw_clawpatch3_copy_face;
    patch_vt->average_face         = fclaw_clawpatch3_average_face;
    patch_vt->interpolate_face     = fclaw_clawpatch3_interpolate_face;

    patch_vt->copy_corner          = fclaw_clawpatch3_copy_corner;
    patch_vt->average_corner       = fclaw_clawpatch3_average_corner;
    patch_vt->interpolate_corner   = fclaw_clawpatch3_interpolate_corner;

    /* Regridding  functions */
    patch_vt->tag4refinement    = &fclaw_clawpatch3_tag4refinement;
    patch_vt->tag4coarsening    = &fclaw_clawpatch3_tag4coarsening;

    patch_vt->average2coarse    = &fclaw_clawpatch3_average2coarse;
    patch_vt->interpolate2fine  = &fclaw_clawpatch3_interpolate2fine;

    /* Defaults for writing output */
    patch_vt->write_header             = &fclaw_clawpatch3_output_ascii_header;
    patch_vt->write_file               = &fclaw_clawpatch3_output_ascii;

    /* Time interpolation functions */
    patch_vt->setup_timeinterp         = &fclaw_clawpatch3_setup_timeinterp;

    /* ghost patch */
    patch_vt->ghost_pack        = &fclaw_clawpatch3_ghost_pack;
    patch_vt->ghost_unpack      = &fclaw_clawpatch3_ghost_unpack;
    patch_vt->build_ghost       = &fclaw_clawpatch3_build_ghost;
    patch_vt->ghost_packsize    = &fclaw_clawpatch3_ghost_packsize;
    patch_vt->local_ghost_alloc = &fclaw_clawpatch3_local_ghost_alloc;
    patch_vt->local_ghost_free  = &fclaw_clawpatch3_local_ghost_free;

    /* partitioning */
    patch_vt->partition_pack    = &fclaw_clawpatch3_partition_pack;
    patch_vt->partition_unpack  = &fclaw_clawpatch3_partition_unpack;
    patch_vt->partition_packsize= &fclaw_clawpatch3_partition_packsize;

    /* diagnostic functions that apply to patches (error, conservation) */
    diag_vt->patch_init_diagnostics      = &fclaw_clawpatch3_diagnostics_initialize;
    diag_vt->patch_compute_diagnostics   = &fclaw_clawpatch3_diagnostics_compute;
    diag_vt->patch_gather_diagnostics    = &fclaw_clawpatch3_diagnostics_gather;
    diag_vt->patch_reset_diagnostics     = &fclaw_clawpatch3_diagnostics_reset;
    diag_vt->patch_finalize_diagnostics  = &fclaw_clawpatch3_diagnostics_finalize;

    patch_vt->defaults_set = 1;
}
