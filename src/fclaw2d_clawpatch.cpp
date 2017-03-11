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
#include <fclaw2d_clawpatch.hpp>

#include <fclaw2d_timeinterp.h>
#include <fclaw2d_metric.h>
#include <fclaw2d_neighbors_fort.h>


/* Store virtual table for all clawpatch routines */
static fclaw2d_clawpatch_vtable_t s_clawpatch_vt;

static
void setup_area_storage(fclaw2d_clawpatch_t* cp);

static
void setup_metric_storage(fclaw2d_clawpatch_t* cp);

static
void ghost_comm(fclaw2d_domain_t* domain,
                fclaw2d_patch_t* this_patch,
                double *qpack, int time_interp,
                int packmode);

static
void metric_setup(fclaw2d_domain_t* domain,
                  fclaw2d_patch_t* this_patch,
                  int blockno,
                  int patchno);

static
fclaw2d_clawpatch_t* clawpatch_data(fclaw2d_patch_t *this_patch);

static
double* q_time_sync(fclaw2d_clawpatch_t* cp, int time_interp);

static
fclaw2d_clawpatch_vtable_t* clawpatch_vt()
{
    // FCLAW_ASSERT(s_clawpatch_vt.defaults_set);
    return &s_clawpatch_vt;
}

fclaw2d_clawpatch_vtable_t* fclaw2d_clawpatch_vt()
{
    /* This isn't a great check, since we never the 'set_defaults' to anything initially.  It is
       probably just getting called every time in this function */
    // FCLAW_ASSERT(s_clawpatch_vt.defaults_set);
    return &s_clawpatch_vt;
}


void fclaw2d_clawpatch_link_app(fclaw_app_t* app)
{
    fclaw2d_clawpatch_t::app = app;
}

void fclaw2d_clawpatch_link_global (fclaw2d_global_t * global)
{
    fclaw2d_clawpatch_t::global = global;
}

/* ------------------------------------------------------------
   Solution access functions
   ---------------------------------------------------------- */
void* fclaw2d_clawpatch_new_patch()
{
    fclaw2d_clawpatch_t *cp = new fclaw2d_clawpatch_t;
    cp->package_data_ptr = fclaw_package_data_new();
    // cp->clawp = new fclaw2d_clawpatch_t();
    return (void*) cp;
}

void fclaw2d_clawpatch_delete_patch(void *user_patch)
{
    FCLAW_ASSERT(user_patch != NULL);
    fclaw2d_clawpatch_t* cp = (fclaw2d_clawpatch_t*) user_patch;
    fclaw_package_patch_data_destroy(fclaw2d_clawpatch_t::app,
                                     cp->package_data_ptr);

    fclaw_package_data_destroy(cp->package_data_ptr);
    user_patch = NULL;
}

/* Access clawpatch stored in patch->user_data */
static
fclaw2d_clawpatch_t* clawpatch_data(fclaw2d_patch_t *this_patch)
{
    fclaw2d_clawpatch_t *cp = (fclaw2d_clawpatch_t*) fclaw2d_patch_get_user_patch(this_patch);
    return cp;
}


fclaw2d_clawpatch_t* fclaw2d_clawpatch_get_cp(fclaw2d_patch_t* this_patch)

{
    return clawpatch_data(this_patch);
}

void fclaw2d_clawpatch_grid_data(fclaw2d_domain_t* domain,
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

void fclaw2d_clawpatch_soln_data(fclaw2d_domain_t* domain,
                                 fclaw2d_patch_t* this_patch,
                                 double **q, int* meqn)
{
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    *q = cp->griddata.dataPtr();
    *meqn = cp->meqn;
}

double *fclaw2d_clawpatch_get_q(fclaw2d_domain_t* domain,
                                fclaw2d_patch_t* this_patch)
{
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    return cp->griddata.dataPtr();
}

double* fclaw2d_clawpatch_get_error(fclaw2d_domain_t* domain,
                                    fclaw2d_patch_t* this_patch)
{
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    return cp->griderror.dataPtr();
}




size_t fclaw2d_clawpatch_size(fclaw2d_domain_t *domain)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int meqn = gparms->meqn;
    int mbc = gparms->mbc;
    size_t size = (mx+2*mbc)*(my+2*mbc)*meqn;

    return size;
}

/* ------------------------------------------------
   Functions for handling time-interpolated data
   ------------------------------------------------ */

void fclaw2d_clawpatch_setup_timeinterp(fclaw2d_domain_t* domain,
                                        fclaw2d_patch_t *this_patch,
                                        double alpha)
{
    /* We use the pack size here to make sure we are setting
       everything correctly;  it isn't needed for memory
       allocation */
    const amr_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int meqn = gparms->meqn;
    int mbc = gparms->mbc;
    int mint = gparms->interp_stencil_width/2+1;  /* Assume interp stencils have odd width */

    int hole = (mx - 2*mint)*(my - 2*mint);  /* Hole in center */
    FCLAW_ASSERT(hole >= 0);

    int wg = mx*my;  /* Whole grid but no ghost cells.  Ghost cells will be averaged from finer
                      * level. */
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

    clawpatch_vt()->fort_timeinterp(&mx,&my,&mbc,&meqn,&psize,
                                    qcurr,qlast,qinterp,&alpha,&ierror);

}


void fclaw2d_clawpatch_timesync_data(fclaw2d_domain_t* domain,
                                     fclaw2d_patch_t* this_patch,
                                     fclaw_bool time_interp,
                                     double **q, int* meqn)
{
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    *q = q_time_sync(cp, time_interp);
    *meqn = cp->meqn;
}

double *fclaw2d_clawpatch_get_q_timesync(fclaw2d_domain_t* domain,
                                         fclaw2d_patch_t* this_patch,
                                         int time_interp)
{
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    return q_time_sync(cp, time_interp);
}



/* This is called from libraries routines (clawpack4.6, clawpack5, etc) */
void fclaw2d_clawpatch_save_current_step(fclaw2d_domain_t* domain,
                                         fclaw2d_patch_t* this_patch)
{
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    cp->griddata_last = cp->griddata;
}

/* This is called from libraries routines (clawpack4.6, clawpack5, etc) */
void fclaw2d_clawpatch_restore_step(fclaw2d_domain_t* domain,
                                    fclaw2d_patch_t* this_patch)
{
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    cp->griddata = cp->griddata_save;
}

/* This is called from libraries routines (clawpack4.6, clawpack5, etc) */
void fclaw2d_clawpatch_save_step(fclaw2d_domain_t* domain,
                                 fclaw2d_patch_t* this_patch)
{
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    cp->griddata_save = cp->griddata;
}


/* ------------------------------------------------------------------
   Manifold setup and access
 ------------------------------------------------------------------ */

static
void metric_setup(fclaw2d_domain_t* domain,
                  fclaw2d_patch_t* this_patch,
                  int blockno,
                  int patchno)
{
    /* vt.patch_manifold_setup_mesh(...) */
    fclaw2d_metric_setup_mesh(domain,this_patch,blockno,patchno);

    /* vt.patch_manifold_compute_normals(...) */
    fclaw2d_metric_compute_normals(domain,this_patch,blockno,patchno);
}


void fclaw2d_clawpatch_metric_data(fclaw2d_domain_t* domain,
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

void fclaw2d_clawpatch_metric_data2(fclaw2d_domain_t* domain,
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

double* fclaw2d_clawpatch_get_area(fclaw2d_domain_t* domain,
                                   fclaw2d_patch_t* this_patch)
{
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    return cp->area.dataPtr();
}

/* ------------------------------------------------------------------
   Re-partition and rebuild new domains, or construct initial domain
 -------------------------------------------------------------------- */

void fclaw2d_clawpatch_define(fclaw2d_domain_t* domain,
                              fclaw2d_patch_t *this_patch,
                              int blockno, int patchno,
                              fclaw2d_build_mode_t build_mode)
{
    /* We are getting closer to getting rid the class fclaw2d_clawpatch_t */
    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);

    const amr_options_t *gparms = get_domain_parms(domain);

    cp->mx = gparms->mx;
    cp->my = gparms->my;
    cp->mbc = gparms->mbc;
    cp->blockno = blockno;
    cp->meqn = gparms->meqn;
    for (int icorner=0; icorner < 4; icorner++)
    {
        fclaw2d_patch_set_block_corner_count(domain,this_patch,icorner,0);
    }

    fclaw2d_map_context_t* cont =
        fclaw2d_domain_get_map_context(domain);

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

    int ll[SpaceDim];
    int ur[SpaceDim];
    for (int idir = 0; idir < SpaceDim; idir++)
    {
        ll[idir] = 1-cp->mbc;
    }
    ur[0] = cp->mx + cp->mbc;
    ur[1] = cp->my + cp->mbc;
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

    // Set up storage for metric terms, if needed.
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


    fclaw_package_patch_data_new(fclaw2d_clawpatch_t::app,cp->package_data_ptr);

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

void fclaw2d_clawpatch_build(fclaw2d_domain_t *domain,
                             fclaw2d_patch_t *this_patch,
                             int blockno,
                             int patchno,
                             void *user)
{
    // fclaw2d_vtable_t vt = fclaw2d_get_vtable(domain);

    fclaw2d_build_mode_t build_mode =  *((fclaw2d_build_mode_t*) user);
    const amr_options_t *gparms = get_domain_parms(domain);

    fclaw2d_clawpatch_define(domain,this_patch,blockno,patchno,build_mode);

    if (gparms->manifold)
    {
        fclaw2d_vt()->metric_compute_area(domain,this_patch,blockno,patchno);
        metric_setup(domain,this_patch,blockno,patchno);
    }
}

void fclaw2d_clawpatch_build_from_fine(fclaw2d_domain_t *domain,
                                       fclaw2d_patch_t *fine_patches,
                                       fclaw2d_patch_t *coarse_patch,
                                       int blockno,
                                       int coarse_patchno,
                                       int fine0_patchno,
                                       fclaw2d_build_mode_t build_mode)
{

    const amr_options_t *gparms = get_domain_parms(domain);

    fclaw2d_clawpatch_define(domain,coarse_patch,blockno,coarse_patchno,build_mode);

    if (gparms->manifold)
    {
        /* Don't recompute the area, but rather average from finer areas */
        fclaw2d_metric_average_area(domain,fine_patches,coarse_patch,
                                    blockno, coarse_patchno, fine0_patchno);

        metric_setup(domain,coarse_patch,blockno,
                                       coarse_patchno);
    }
}

void fclaw2d_clawpatch_build_ghost(fclaw2d_domain_t *domain,
                                   fclaw2d_patch_t *this_patch,
                                   int blockno,
                                   int patchno,
                                   void *user)
{
    // fclaw2d_vtable_t vt = fclaw2d_get_vtable(domain);

    fclaw2d_build_mode_t build_mode =  *((fclaw2d_build_mode_t*) user);
    const amr_options_t *gparms = get_domain_parms(domain);

    fclaw2d_clawpatch_define(domain,this_patch,blockno,patchno,build_mode);

    if (gparms->manifold)
    {
        if (build_mode != FCLAW2D_BUILD_FOR_GHOST_AREA_PACKED)
        {
            fclaw2d_vt()->metric_compute_area(domain,this_patch,blockno,patchno);
        }
    }
}

/* --------------------------------------------------------------
   Parallel ghost exchanges.

   Note this is different from the packing/unpacking that happens
   when partitioning the domain
   --------------------------------------------------------------*/

size_t fclaw2d_clawpatch_ghost_packsize(fclaw2d_domain_t* domain)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;
    int refratio = gparms->refratio;

    int mint = refratio*mbc;
    int nghost = mbc;

    int wg = (2*nghost + mx)*(2*nghost + my);  /* Whole grid     */
    int hole = (mx - 2*mint)*(my - 2*mint);    /* Hole in center */
    FCLAW_ASSERT(hole >= 0);

    int packarea = gparms->ghost_patch_pack_area && gparms->manifold;
    int packextra = gparms->ghost_patch_pack_numextrafields;
    int nfields = meqn + packarea + packextra;
    size_t psize = (wg - hole)*nfields;
    FCLAW_ASSERT(psize > 0);

    return psize*sizeof(double);
}

void fclaw2d_clawpatch_local_ghost_alloc(fclaw2d_domain_t* domain,
                                           fclaw2d_patch_t* this_patch,
                                           void** q)
{
    /* Create contiguous block for data and area */
    int msize = fclaw2d_clawpatch_ghost_packsize(domain);
    *q = (void*) FCLAW_ALLOC(double,msize);
    FCLAW_ASSERT(*q != NULL);

}

void fclaw2d_clawpatch_local_ghost_free(fclaw2d_domain_t* domain,
                                        void **q)
{
    FCLAW_FREE(*q);
    *q = NULL;
}


static
void ghost_comm(fclaw2d_domain_t* domain,
                          fclaw2d_patch_t* this_patch,
                          double *qpack, int time_interp,
                          int packmode)
{
    int meqn;
    double *qthis;
    double *area;
    const amr_options_t *gparms = get_domain_parms(domain);

    int ierror;

    int packarea = packmode/2;   // (0,1)/2 = 0;  (2,3)/2 = 1;

    fclaw2d_clawpatch_timesync_data(domain,this_patch,time_interp,&qthis,&meqn);
    area = fclaw2d_clawpatch_get_area(domain,this_patch);

    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int refratio = gparms->refratio;

    int mint = mbc*refratio;   /* # interior cells needed for averaging */
    int nghost = mbc;          /* # ghost values needed for interpolation */

    /* This is computed twice - here, and in fclaw2d_clawpatch_ghost_packsize */
    int wg = (2*nghost + mx)*(2*nghost + my);
    int hole = (mx - 2*mint)*(my - 2*mint);  /* Hole in center */
    FCLAW_ASSERT(hole >= 0);

    int psize = (wg - hole)*(meqn + packarea + gparms->ghost_patch_pack_numextrafields);
    FCLAW_ASSERT(psize > 0);

    int qareasize = (wg - hole)*(meqn + packarea);
    clawpatch_vt()->fort_ghostpack_qarea(&mx,&my,&mbc,&meqn,&mint,qthis,area,
                                      qpack,&qareasize,&packmode,&ierror);
    FCLAW_ASSERT(ierror == 0);
    if (gparms->ghost_patch_pack_extra)
    {
      qpack += qareasize;
      int extrasize = psize - qareasize;
      FCLAW_ASSERT(clawpatch_vt()->ghostpack_extra != NULL);
      clawpatch_vt()->ghostpack_extra(domain,this_patch,mint,qpack,extrasize,packmode,&ierror);
      FCLAW_ASSERT(ierror == 0);
    }

    if (ierror > 0)
    {
        fclaw_global_essentialf("ghost_pack (fclaw2d_clawpatch.cpp) : ierror = %d\n",ierror);
        exit(0);
    }
}

void fclaw2d_clawpatch_ghost_pack(fclaw2d_domain_t *domain,
                                  fclaw2d_patch_t *this_patch,
                                  double *patch_data,
                                  int time_interp)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    int packarea = gparms->ghost_patch_pack_area && gparms->manifold;
    int packmode = 2*packarea;  // 0 or 2  (for pack)

    ghost_comm(domain,this_patch,patch_data, time_interp,packmode);
}


void fclaw2d_clawpatch_ghost_unpack(fclaw2d_domain_t* domain,
                                    fclaw2d_patch_t* this_patch,
                                    int this_block_idx,
                                    int this_patch_idx,
                                    double *qdata, fclaw_bool time_interp)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    int packarea = gparms->ghost_patch_pack_area && gparms->manifold;
    int packmode = 2*packarea + 1;  // 1 or 3  (for unpack)

    ghost_comm(domain,this_patch,qdata, time_interp,packmode);
}

/* --------------------------------------------------------
   Domain partitioning.

   Pack local patches on this processor before re-partitioning;
   retrieve patches that migrated across parallel borders to
   load balance the computation.
   -------------------------------------------------------*/

size_t fclaw2d_clawpatch_partition_packsize(fclaw2d_domain_t* domain)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;
    size_t size = (2*mbc + mx)*(2*mbc + my)*meqn;  /* Store area */
    return size*sizeof(double);
}

void fclaw2d_clawpatch_partition_pack(fclaw2d_domain_t *domain,
                                         fclaw2d_patch_t *this_patch,
                                         int this_block_idx,
                                         int this_patch_idx,
                                         void *user)
{
    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    int patch_num = this_block->num_patches_before + this_patch_idx;
    double* patch_data = (double*) ((void**)user)[patch_num];

    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);
    FCLAW_ASSERT(cp != NULL);

    cp->griddata.copyToMemory(patch_data);
}

void fclaw2d_clawpatch_partition_unpack(fclaw2d_domain_t *domain,
                                        fclaw2d_patch_t *this_patch,
                                        int this_block_idx,
                                        int this_patch_idx,
                                        void *user)
{
    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    int patch_num = this_block->num_patches_before + this_patch_idx;
    double* patch_data = (double*) ((void**)user)[patch_num];

    fclaw2d_clawpatch_t *cp = clawpatch_data(this_patch);

    /* Time interp is false, since we only partition when all grids
       are time synchronized */
    cp->griddata.copyFromMemory(patch_data);
}

/* ----------------------------------------------------------------
   Set defaults for clawpatch virtual table
   ---------------------------------------------------------------- */


void fclaw2d_clawpatch_init_vtable_defaults()
{

    fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();

    /* These must be redefined by the solver and user */
    patch_vt->patch_initialize         = NULL;
    patch_vt->patch_physical_bc        = NULL;
    patch_vt->patch_single_step_update = NULL;

    /* These may be redefined by the user */
    /* Problem setup */
    patch_vt->patch_new             = fclaw2d_clawpatch_new_patch;
    patch_vt->patch_delete          = fclaw2d_clawpatch_delete_patch;
    patch_vt->patch_setup           = NULL;
    patch_vt->patch_setup_ghost     = NULL;
    patch_vt->patch_build           = &fclaw2d_clawpatch_build;
    patch_vt->patch_build_from_fine = &fclaw2d_clawpatch_build_from_fine;
    patch_vt->patch_restore_step    = &fclaw2d_clawpatch_restore_step;
    patch_vt->patch_save_step       = &fclaw2d_clawpatch_save_step;

    /* Ghost filling - solver specific */
    patch_vt->copy_face            = fclaw2d_clawpatch_copy_face;
    patch_vt->average_face         = fclaw2d_clawpatch_average_face;
    patch_vt->interpolate_face     = fclaw2d_clawpatch_interpolate_face;

    patch_vt->copy_corner          = fclaw2d_clawpatch_copy_corner;
    patch_vt->average_corner       = fclaw2d_clawpatch_average_corner;
    patch_vt->interpolate_corner   = fclaw2d_clawpatch_interpolate_corner;

    /* Regridding  functions */
    patch_vt->patch_tag4refinement    = &fclaw2d_clawpatch_tag4refinement;
    patch_vt->patch_tag4coarsening    = &fclaw2d_clawpatch_tag4coarsening;

    patch_vt->patch_average2coarse    = &fclaw2d_clawpatch_average2coarse;
    patch_vt->patch_interpolate2fine  = &fclaw2d_clawpatch_interpolate2fine;

    /* Defaults for writing output */
    patch_vt->write_header             = &fclaw2d_clawpatch_output_ascii_header;
    patch_vt->patch_write_file         = &fclaw2d_clawpatch_output_ascii;

    /* Time interpolation functions */
    patch_vt->setup_timeinterp         = &fclaw2d_clawpatch_setup_timeinterp;

    /* ghost patch */
    patch_vt->ghost_pack        = &fclaw2d_clawpatch_ghost_pack;
    patch_vt->ghost_unpack      = &fclaw2d_clawpatch_ghost_unpack;
    patch_vt->build_ghost       = &fclaw2d_clawpatch_build_ghost;
    patch_vt->ghost_packsize    = &fclaw2d_clawpatch_ghost_packsize;
    patch_vt->local_ghost_alloc = &fclaw2d_clawpatch_local_ghost_alloc;
    patch_vt->local_ghost_free  = &fclaw2d_clawpatch_local_ghost_free;

    /* partitioning */
    patch_vt->partition_pack    = &fclaw2d_clawpatch_partition_pack;
    patch_vt->partition_unpack  = &fclaw2d_clawpatch_partition_unpack;
    patch_vt->partition_packsize= &fclaw2d_clawpatch_partition_packsize;

    patch_vt->compute_error     = &fclaw2d_clawpatch_compute_error;
    patch_vt->compute_norm      = &fclaw2d_clawpatch_compute_error_norm;

    patch_vt->defaults_set = 1;
    // s_clawpatch_vt.defaults_set = 1;
}


void fclaw2d_clawpatch_copy_face(fclaw2d_domain_t *domain,
                                 fclaw2d_patch_t *this_patch,
                                 fclaw2d_patch_t *neighbor_patch,
                                 int iface,
                                 int time_interp,
                                 fclaw2d_transform_data_t *transform_data)

{
    int meqn,mx,my,mbc;
    double *qthis, *qneighbor;
    const amr_options_t *gparms = get_domain_parms(domain);

    mx = gparms->mx;
    my = gparms->my;
    mbc = gparms->mbc;

    fclaw2d_clawpatch_timesync_data(domain,this_patch,time_interp,&qthis,&meqn);
    fclaw2d_clawpatch_timesync_data(domain,neighbor_patch,time_interp,&qneighbor,&meqn);

    clawpatch_vt()->fort_copy_face(&mx,&my,&mbc,&meqn,qthis,qneighbor,&iface,&transform_data);
}

void fclaw2d_clawpatch_average_face(fclaw2d_domain_t *domain,
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
    int meqn,mx,my,mbc;
    double *qcoarse, *qfine;
    double *areacoarse, *areafine;
    const amr_options_t *gparms = get_domain_parms(domain);

    fclaw2d_clawpatch_timesync_data(domain,coarse_patch,time_interp,&qcoarse,&meqn);
    qfine = fclaw2d_clawpatch_get_q(domain,fine_patch);

    /* These will be empty for non-manifolds cases */
    areacoarse = fclaw2d_clawpatch_get_area(domain,coarse_patch);
    areafine = fclaw2d_clawpatch_get_area(domain,fine_patch);

    mx = gparms->mx;
    my = gparms->my;
    mbc = gparms->mbc;

    int manifold = gparms->manifold;
    clawpatch_vt()->fort_average_face(&mx,&my,&mbc,&meqn,qcoarse,qfine,areacoarse,areafine,
                                      &idir,&iface_coarse, &p4est_refineFactor, &refratio,
                                      &igrid,&manifold,&transform_data);


}

void fclaw2d_clawpatch_interpolate_face(fclaw2d_domain_t *domain,
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
    int meqn,mx,my,mbc;
    double *qcoarse, *qfine;
    const amr_options_t *gparms = get_domain_parms(domain);

    fclaw2d_clawpatch_timesync_data(domain,coarse_patch,time_interp,&qcoarse,&meqn);
    qfine = fclaw2d_clawpatch_get_q(domain,fine_patch);

    mx = gparms->mx;
    my = gparms->my;
    mbc = gparms->mbc;

    clawpatch_vt()->fort_interpolate_face(&mx,&my,&mbc,&meqn,qcoarse,qfine,&idir,&iside,
                                          &p4est_refineFactor,&refratio,&igrid,&transform_data);
}


void fclaw2d_clawpatch_copy_corner(fclaw2d_domain_t *domain,
                                   fclaw2d_patch_t *this_patch,
                                   fclaw2d_patch_t *corner_patch,
                                   int icorner,
                                   int time_interp,
                                   fclaw2d_transform_data_t *transform_data)
{
    int meqn,mx,my,mbc;
    double *qthis, *qcorner;
    const amr_options_t *gparms = get_domain_parms(domain);

    mx = gparms->mx;
    my = gparms->my;
    mbc = gparms->mbc;

    fclaw2d_clawpatch_timesync_data(domain,this_patch,time_interp,&qthis,&meqn);
    fclaw2d_clawpatch_timesync_data(domain,corner_patch,time_interp,&qcorner,&meqn);

    clawpatch_vt()->fort_copy_corner(&mx,&my,&mbc,&meqn,qthis,qcorner,&icorner,&transform_data);

}

void fclaw2d_clawpatch_average_corner(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *coarse_patch,
                                      fclaw2d_patch_t *fine_patch,
                                      int coarse_corner,
                                      int refratio,
                                      fclaw_bool time_interp,
                                      fclaw2d_transform_data_t* transform_data)
{
    int meqn,mx,my,mbc;
    double *qcoarse, *qfine;
    double *areacoarse, *areafine;
    const amr_options_t *gparms = get_domain_parms(domain);

    fclaw2d_clawpatch_timesync_data(domain,coarse_patch,time_interp,&qcoarse,&meqn);
    qfine = fclaw2d_clawpatch_get_q(domain,fine_patch);

    areacoarse = fclaw2d_clawpatch_get_area(domain,coarse_patch);
    areafine = fclaw2d_clawpatch_get_area(domain,fine_patch);

    mx = gparms->mx;
    my = gparms->my;
    mbc = gparms->mbc;

    int manifold = gparms->manifold;
    clawpatch_vt()->fort_average_corner(&mx,&my,&mbc,&meqn,&refratio,
                           qcoarse,qfine,areacoarse,areafine,
                           &manifold,&coarse_corner,&transform_data);
}


void fclaw2d_clawpatch_interpolate_corner(fclaw2d_domain_t* domain,
                                          fclaw2d_patch_t* coarse_patch,
                                          fclaw2d_patch_t* fine_patch,
                                          int coarse_corner,
                                          int refratio,
                                          fclaw_bool time_interp,
                                          fclaw2d_transform_data_t* transform_data)

{
    int meqn,mx,my,mbc;
    double *qcoarse, *qfine;
    const amr_options_t *gparms = get_domain_parms(domain);

    mx = gparms->mx;
    my = gparms->my;
    mbc = gparms->mbc;
    meqn = gparms->meqn;


    fclaw2d_clawpatch_timesync_data(domain,coarse_patch,time_interp,&qcoarse,&meqn);
    qfine = fclaw2d_clawpatch_get_q(domain,fine_patch);

    clawpatch_vt()->fort_interpolate_corner(&mx,&my,&mbc,&meqn,
                                         &refratio,qcoarse,qfine,
                                         &coarse_corner,&transform_data);
}




/* ----------------------------------------------------------------
   fclaw2d_clawpatch_t member functions (previous in fclaw2d_clawpatch_t.cpp)
   ---------------------------------------------------------------- */

fclaw_app_t *fclaw2d_clawpatch_t::app;
fclaw2d_global_t *fclaw2d_clawpatch_t::global;

Box fclaw2d_clawpatch_t::dataBox()
{
    return griddata.box();
}

Box fclaw2d_clawpatch_t::areaBox()
{
    return area.box();
}

Box fclaw2d_clawpatch_t::edgeBox()
{
    return edge_lengths.box();
}

Box fclaw2d_clawpatch_t::nodeBox()
{
    return xp.box();
}


/* Return a pointer to either time interpolated data or regular grid data */
static double* q_time_sync(fclaw2d_clawpatch_t* cp, int time_interp)
{
    if (time_interp)
        return cp->griddata_time_interpolated.dataPtr();
    else
        return cp->griddata.dataPtr();
}

/* ----------------------------------------------------
   Solver data and functions
   ---------------------------------------------------*/
// Wave propagation algorithms

void* fclaw2d_clawpatch_t::clawpack_patch_data(int id)
{
    return fclaw_package_get_data(package_data_ptr,id);
}

/* ----------------------------------------------------------------
   Special case : Pillow grid ghost exchanges/average/interpolate
   ---------------------------------------------------------------- */

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


void fclaw2d_clawpatch_t::mb_interpolate_block_corner_ghost(const int& a_coarse_corner,
                                                            const int& a_refratio,
                                                            fclaw2d_clawpatch_t *cp_corner,
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

static void setup_area_storage(fclaw2d_clawpatch_t* cp)
{
    int mx = cp->mx;
    int my = cp->my;
    int mbc = cp->mbc;
xf
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

static void setup_metric_storage(fclaw2d_clawpatch_t* cp)
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
