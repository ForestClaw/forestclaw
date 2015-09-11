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
#include <fclaw2d_metric_default.h>
#include <fclaw2d_timeinterp.h>
#include <fclaw2d_metric.h>

#include <fclaw2d_timeinterp.h>
#include <fclaw2d_ghost_fill.h>
#include <fclaw2d_neighbors_fort.h>
#include <fclaw2d_metric_default_fort.h>
#include <fclaw2d_map_query.h>

#include <fclaw2d_clawpatch.hpp>



#if FCLAW_DEBUG
/* This is a duplicate of the function in fclaw2d_farraybox.cpp */
static
void set_snan(double& f)
{
    /* From :
      "NaNs, Uninitialized Variables, and C++"
      http://codingcastles.blogspot.fr/2008/12/nans-in-c.html
    */
    *((long long*)&f) = 0x7ff0000000000001LL;
}
#endif

void fclaw2d_clawpatch_link_app(fclaw_app_t* app)
{
    ClawPatch::app = app;
}

void fclaw2d_clawpatch_link_global (fclaw2d_global_t * global)
{
    ClawPatch::global = global;
}

/* ------------------------------------------------------------
   Solution access functions
   ---------------------------------------------------------- */

/* We also have a fclaw2d_patch_get_cp(this_patch) */
ClawPatch* fclaw2d_clawpatch_get_cp(fclaw2d_patch_t* this_patch)

{
    return fclaw2d_patch_get_cp(this_patch);
}

void fclaw2d_clawpatch_grid_data(fclaw2d_domain_t* domain,
                                 fclaw2d_patch_t* this_patch,
                                 int* mx, int* my, int* mbc,
                                 double* xlower, double* ylower,
                                 double* dx, double* dy)
{
    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    *mx = cp->mx();
    *my = cp->my();
    *mbc = cp->mbc();
    *xlower = cp->xlower();
    *ylower = cp->ylower();
    *dx = cp->dx();
    *dy = cp->dy();
}

void fclaw2d_clawpatch_soln_data(fclaw2d_domain_t* domain,
                                 fclaw2d_patch_t* this_patch,
                                 double **q, int* meqn)
{
    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    *q = cp->q();
    *meqn = cp->meqn();
}

double *fclaw2d_clawpatch_get_q(fclaw2d_domain_t* domain,
                                fclaw2d_patch_t* this_patch)
{
    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    return cp->q();
}

static
void cb_has_finegrid_neighbors(fclaw2d_domain_t* domain,
                               fclaw2d_patch_t *this_patch,
                               int blockno,
                               int patchno,
                               void* user)
{
    int y =
        fclaw2d_timeinterp_has_finegrid_neighbors(domain,blockno,patchno);
    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    cp->finegrid_neighbors(y);
}

void fclaw2d_clawpatch_finegrid_neighbors(fclaw2d_domain_t* domain)
{
    fclaw2d_domain_iterate_patches(domain,cb_has_finegrid_neighbors,
                                   (void*) NULL);
}

int fclaw2d_clawpatch_has_finegrid_neighbors(fclaw2d_domain_t* domain,
                                             fclaw2d_patch_t* this_patch)
{
    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    return cp->has_finegrid_neighbors();
}


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
    int mint = 2;  /* Number of internal layers needed */

    int ng = 0;  /* Number of ghost cells (mbc,1,2,..) */
    int wg = (2*ng + mx)*(2*ng + my);  /* Whole grid  (one layer of ghost cells)*/
    int hole = (mx - 2*mint)*(my - 2*mint);  /* Hole in center */
    FCLAW_ASSERT(hole >= 0);
    size_t psize = (wg - hole)*meqn;
    FCLAW_ASSERT(psize > 0);

    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    cp->setup_for_time_interpolation(alpha,psize);
}

void fclaw2d_clawpatch_timesync_data(fclaw2d_domain_t* domain,
                                     fclaw2d_patch_t* this_patch,
                                     fclaw_bool time_interp,
                                     double **q, int* meqn)
{
    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    *q = cp->q_time_sync(time_interp);
    *meqn = cp->meqn();
}

double *fclaw2d_clawpatch_get_q_timesync(fclaw2d_domain_t* domain,
                                         fclaw2d_patch_t* this_patch,
                                         int time_interp)
{
    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    return cp->q_time_sync(time_interp);
}



void fclaw2d_clawpatch_save_current_step(fclaw2d_domain_t* domain,
                                         fclaw2d_patch_t* this_patch)
{
    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    cp->save_current_step();
}

/* ------------------------------------------------------------------
   Miscellaneous
 ------------------------------------------------------------------ */
int* fclaw2d_clawpatch_corner_count(fclaw2d_domain_t* domain,
                                   fclaw2d_patch_t* this_patch)
{
    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    return cp->block_corner_count();
}

/* ------------------------------------------------------------------
   Manifold setup and access
 ------------------------------------------------------------------ */

static
void fclaw2d_clawpatch_metric_setup(fclaw2d_domain_t* domain,
                                      fclaw2d_patch_t* this_patch,
                                      int blockno,
                                      int patchno)
{
    fclaw2d_vtable_t vt;
    vt = fclaw2d_get_vtable(domain);

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
    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    *xp = cp->xp();
    *yp = cp->yp();
    *zp = cp->zp();
    *xd = cp->xd();
    *yd = cp->yd();
    *zd = cp->zd();
    *area = cp->area();
}

void fclaw2d_clawpatch_metric_data2(fclaw2d_domain_t* domain,
                                    fclaw2d_patch_t* this_patch,
                                    double **xnormals, double **ynormals,
                                    double **xtangents, double **ytangents,
                                    double **surfnormals,
                                    double **edgelengths, double **curvature)
{
    /* or just call the member functions? */
    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    *xnormals    = cp->xface_normals();
    *ynormals    = cp->yface_normals();
    *xtangents   = cp->xface_tangents();
    *ytangents   = cp->yface_tangents();
    *surfnormals = cp->surf_normals();
    *curvature   = cp->curvature();
    *edgelengths = cp->edge_lengths();
}

double* fclaw2d_clawpatch_get_area(fclaw2d_domain_t* domain,
                                   fclaw2d_patch_t* this_patch)
{
    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    return cp->area();
}

/* ------------------------------------------------------------------
   Re-partition and rebuild new domains, or construct initial domain
 -------------------------------------------------------------------- */

void fclaw2d_clawpatch_define(fclaw2d_domain_t* domain,
                              fclaw2d_patch_t *this_patch,
                              int blockno, int patchno,
                              fclaw2d_build_mode_t build_mode)
{
    /* We are getting closer to getting rid the class ClawPatch */
    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    cp->define(domain,this_patch,blockno,build_mode);

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
    fclaw2d_vtable_t vt;
    vt = fclaw2d_get_vtable(domain);

    fclaw2d_build_mode_t build_mode =  *((fclaw2d_build_mode_t*) user);
    const amr_options_t *gparms = get_domain_parms(domain);

    fclaw2d_clawpatch_define(domain,this_patch,blockno,patchno,build_mode);

    if (gparms->manifold)
    {
        if (build_mode != FCLAW2D_BUILD_FOR_GHOST_AREA_PACKED)
        {
            vt.metric_compute_area(domain,this_patch,blockno,patchno);

            /* Don't need any more manifold info for ghost patches */
            if (build_mode == FCLAW2D_BUILD_FOR_UPDATE)
            {
                fclaw2d_clawpatch_metric_setup(domain,this_patch,blockno,patchno);
            }
        }
    }

    vt = fclaw2d_get_vtable(domain);
    if (vt.patch_setup != NULL && build_mode == FCLAW2D_BUILD_FOR_UPDATE)
    {
        vt.patch_setup(domain,this_patch,blockno,patchno);
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
    fclaw2d_vtable_t vt;
    vt = fclaw2d_get_vtable(domain);

    const amr_options_t *gparms = get_domain_parms(domain);

    fclaw2d_clawpatch_define(domain,coarse_patch,blockno,coarse_patchno,build_mode);

    if (gparms->manifold)
    {
        fclaw2d_metric_average_area(domain,fine_patches,coarse_patch,
                                    blockno, coarse_patchno, fine0_patchno);

        fclaw2d_clawpatch_metric_setup(domain,coarse_patch,blockno,
                                       coarse_patchno);
    }

    vt = fclaw2d_get_vtable(domain);
    if (vt.patch_setup != NULL && build_mode == FCLAW2D_BUILD_FOR_UPDATE)
    {
        /* We might want to distinguish between new fine grid patches, and
           new coarse grid patches.  In the latter case, we might just average
           aux array info, for example, rather than recreate it from scratch.
           Something like a general "build from fine" routine might be needed */
        vt.patch_setup(domain,coarse_patch,blockno,coarse_patchno);
    }
}


/* -------------------------------------------------
   For debugging
   ----------------------------------------------- */
void cb_set_corners_to_value(fclaw2d_domain_t* domain,
                            fclaw2d_patch_t* this_patch,
                            int blockno,
                            int patchno,
                            void *user)
{
    struct set_bc { int time_interp; double value; };
    struct set_bc s = *((set_bc*) user);
    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    cp->set_corners_to_value(s.time_interp,s.value);
}


void fclaw2d_clawpatch_set_corners_to_value(fclaw2d_domain_t* domain,
                                            int minlevel,
                                            int maxlevel,
                                            int time_interp,
                                            double value)
{
    struct set_bc { int time_interp; double value; };
    struct set_bc s;
    s.time_interp = time_interp;
    s.value = value;

    for (int level = minlevel; level <= maxlevel; level++)
    {
        s.time_interp = 0;
        fclaw2d_domain_iterate_level(domain, level,cb_set_corners_to_value,
                                     (void *)  &s);
    }
    if (time_interp)
    {
        s.time_interp = time_interp;
        int time_interp_minlevel = minlevel-1;
        fclaw2d_domain_iterate_level(domain, time_interp_minlevel,
                                     cb_set_corners_to_value,
                                     (void *)  &s);
    }
}

void fclaw2d_clawpatch_set_corners_to_nan(fclaw2d_domain_t* domain,
                                           int minlevel,
                                           int maxlevel,
                                           int time_interp)
{
    double s;
    fclaw2d_farraybox_set_to_nan(s);
    fclaw2d_clawpatch_set_corners_to_value(domain,minlevel,maxlevel,
                                           time_interp,s);
}



void cb_set_boundary_to_value(fclaw2d_domain_t* domain,
                              fclaw2d_patch_t* this_patch,
                              int blockno,
                              int patchno,
                              void *user)
{
    struct set_bc { int time_interp; double value; };
    struct set_bc s = *((set_bc*) user);
    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    cp->set_boundary_to_value(s.time_interp,s.value);
}



void fclaw2d_clawpatch_set_boundary_to_value(fclaw2d_domain_t* domain,
                                             int minlevel,
                                             int maxlevel,
                                             int time_interp,
                                             double value)
{
    struct set_bc { int time_interp; double value; };
    struct set_bc s;
    s.time_interp = time_interp;
    s.value = value;

    for (int level = minlevel; level <= maxlevel; level++)
    {
        s.time_interp = 0;
        fclaw2d_domain_iterate_level(domain, level,cb_set_boundary_to_value,
                                     (void *)  &s);
    }
    if (time_interp)
    {
        s.time_interp = time_interp;
        int time_interp_minlevel = minlevel-1;
        fclaw2d_domain_iterate_level(domain, time_interp_minlevel,
                                     cb_set_boundary_to_value,
                                     (void *)  &s);
    }
}

void fclaw2d_clawpatch_set_boundary_to_nan(fclaw2d_domain_t* domain,
                                           int minlevel,
                                           int maxlevel,
                                           int time_interp)
{
    double s;
    fclaw2d_farraybox_set_to_nan(s);
    fclaw2d_clawpatch_set_boundary_to_value(domain,minlevel,maxlevel,
                                            time_interp,s);
}



void fclaw2d_clawpatch_initialize_after_regrid(fclaw2d_domain_t* domain,
                                               fclaw2d_patch_t* this_patch,
                                               int this_block_idx,
                                               int this_patch_idx)
{
    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    cp->initialize_after_regrid();
}


/* ----------------------------------------------------------
   Parallel ghost exchanges.

   Note this is different from the packing/unpacking that happens
   when partitioning the domain
   ----------------------------------------------------------*/

size_t fclaw2d_clawpatch_ghost_packsize(fclaw2d_domain_t* domain)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int meqn = gparms->meqn;
    int packarea = gparms->ghost_patch_pack_area && gparms->manifold;

    FCLAW_ASSERT(ClawPatch::pack_layers > 0);
    int mint = ClawPatch::pack_layers;
    int wg = (2 + mx)*(2 + my);  /* Whole grid  (one layer of ghost cells)*/
    int hole = (mx - 2*mint)*(my - 2*mint);  /* Hole in center */
    FCLAW_ASSERT(hole >= 0);
    size_t psize = (wg - hole)*(meqn + packarea);
    FCLAW_ASSERT(psize > 0);
    return psize*sizeof(double);
}

void fclaw2d_clawpatch_ghost_pack_location(fclaw2d_domain_t* domain,
                                           fclaw2d_patch_t* this_patch,
                                           void** q)
{
    /* Create contiguous block for data and area */
    int msize = fclaw2d_clawpatch_ghost_packsize(domain);
    *q = (void*) FCLAW_ALLOC(double,msize);
    FCLAW_ASSERT(*q != NULL);

#if FCLAW_DEBUG
    /* This should go in some kind of debug pre-processing block */
    double snan;
    set_snan(snan);
    double *qq = (double*) *q;
    for(int i = 0; i < msize; i++)
    {
        qq[i] = snan;
    }
#endif
}

void fclaw2d_clawpatch_ghost_free_pack_location(fclaw2d_domain_t* domain,
                                                void **q)
{
    FCLAW_FREE(*q);
    *q = NULL;
}


void fclaw2d_clawpatch_ghost_pack(fclaw2d_domain_t *domain,
                                  fclaw2d_patch_t *this_patch,
                                  double *patch_data,
                                  int time_interp)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    int packarea = gparms->ghost_patch_pack_area && gparms->manifold;
    int packmode = 2*packarea;  // 0 or 2  (for pack)

    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    cp->ghost_comm(patch_data,time_interp,packmode);
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

    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    cp->ghost_comm(qdata,time_interp,packmode);
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

void cb_fclaw2d_clawpatch_partition_pack(fclaw2d_domain_t *domain,
                                         fclaw2d_patch_t *this_patch,
                                         int this_block_idx,
                                         int this_patch_idx,
                                         void *user)
{
    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    int patch_num = this_block->num_patches_before + this_patch_idx;
    double* patch_data = (double*) ((void**)user)[patch_num];

    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    FCLAW_ASSERT(cp != NULL);
    cp->partition_pack(patch_data);
}

void cb_fclaw2d_clawpatch_partition_unpack(fclaw2d_domain_t *domain,
                                           fclaw2d_patch_t *this_patch,
                                           int this_block_idx,
                                           int this_patch_idx,
                                           void *user)
{
    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    int patch_num = this_block->num_patches_before + this_patch_idx;
    double* patch_data = (double*) ((void**)user)[patch_num];

    /* Create new data in 'user' pointer */
    fclaw2d_patch_data_new(domain,this_patch);

    fclaw2d_build_mode_t build_mode = FCLAW2D_BUILD_FOR_UPDATE;
    fclaw2d_clawpatch_build(domain,this_patch,this_block_idx,
                            this_patch_idx,(void*) &build_mode);

    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);

    /* Time interp is false, since we only partition when all grids
       are time synchronized */
    cp->partition_unpack(patch_data);
}

void fclaw2d_clawpatch_initialize_after_partition(fclaw2d_domain_t* domain,
                                                  fclaw2d_patch_t* this_patch,
                                                  int this_block_idx,
                                                  int this_patch_idx)
{
    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    cp->initialize_after_partition();
}


/* ----------------------------------------------------------------
   ClawPatch member functions (previous in ClawPatch.cpp)
   ---------------------------------------------------------------- */

fclaw_app_t *ClawPatch::app;
fclaw2d_global_t *ClawPatch::global;

int ClawPatch::pack_layers = 4;
int ClawPatch::ghost_patch_pack_area = -1;

ClawPatch::ClawPatch()
{
    m_package_data_ptr = fclaw_package_data_new();
}

ClawPatch::~ClawPatch()
{
    fclaw_package_patch_data_destroy(ClawPatch::app,
                                     m_package_data_ptr);

    fclaw_package_data_destroy(m_package_data_ptr);
}

void ClawPatch::define(fclaw2d_domain_t* domain,
                           fclaw2d_patch_t* this_patch,
                           int blockno,
                           fclaw2d_build_mode_t build_mode)

{
    const amr_options_t *gparms = get_domain_parms(domain);

    m_mx = gparms->mx;
    m_my = gparms->my;
    m_mbc = gparms->mbc;
    m_blockno = blockno;
    m_meqn = gparms->meqn;
    for (int m=0; m < 4; m++)
    {
        m_block_corner_count[m] = 0;
    }
    m_has_finegrid_neighbors = 0;

    fclaw2d_map_context_t* cont =
        fclaw2d_domain_get_map_context(domain);

    int is_brick = FCLAW2D_MAP_IS_BRICK(&cont);

    m_manifold = gparms->manifold;

    if (m_manifold)
    {
        m_xlower = this_patch->xlower;
        m_ylower = this_patch->ylower;
        m_xupper = this_patch->xupper;
        m_yupper = this_patch->yupper;
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
        double xlower,ylower;
        double xupper,yupper;

        if (is_brick)
        {
            double z;
            /* Scale to [0,1]x[0,1], based on blockno */
            fclaw2d_map_c2m_nomap_brick(cont,blockno,xl,yl,&xlower,&ylower,&z);
            fclaw2d_map_c2m_nomap_brick(cont,blockno,xu,yu,&xupper,&yupper,&z);
        }
        else
        {
            xlower = xl;
            ylower = yl;
            xupper = xu;
            yupper = yu;
        }

        m_xlower = ax + (bx - ax)*xlower;
        m_xupper = ax + (bx - ax)*xupper;
        m_ylower = ay + (by - ay)*ylower;
        m_yupper = ay + (by - ay)*yupper;
    }

    m_dx = (m_xupper - m_xlower)/m_mx;
    m_dy = (m_yupper - m_ylower)/m_my;

    int ll[SpaceDim];
    int ur[SpaceDim];
    for (int idir = 0; idir < SpaceDim; idir++)
    {
        ll[idir] = 1-m_mbc;
    }
    ur[0] = m_mx + m_mbc;
    ur[1] = m_my + m_mbc;
    Box box(ll,ur);

    // This will destroy any existing memory n m_griddata.
    m_griddata.define(box, m_meqn);
    m_griddata_time_interpolated.define(box, m_meqn);

    // Set up storage for metric terms, if needed.
    if (gparms->manifold)
    {
        setup_area_storage();
        if (build_mode != FCLAW2D_BUILD_FOR_GHOST_AREA_PACKED)
        {
            /* Don't need any more manifold info for ghost patches */
            if (build_mode == FCLAW2D_BUILD_FOR_UPDATE)
            {
                setup_metric_storage();
            }
        }
    }

    fclaw_package_patch_data_new(ClawPatch::app,m_package_data_ptr);
    ClawPatch::pack_layers = 4;

    if (build_mode != FCLAW2D_BUILD_FOR_UPDATE)
    {
        return;
    }

    m_griddata_last.define(box, m_meqn);
    m_griddata_save.define(box, m_meqn);

}

void ClawPatch::initialize_after_regrid()
{
    m_griddata_last.set_to_nan();
    m_griddata_save.set_to_nan();
    m_griddata_time_interpolated.set_to_nan();
}

void ClawPatch::initialize_after_partition()
{
    m_griddata_last.set_to_nan();
    m_griddata_save.set_to_nan();
    m_griddata_time_interpolated.set_to_nan();
}


void ClawPatch::copyFrom(ClawPatch *a_cp)
{
    m_griddata = a_cp->m_griddata;
}

// This is used by level_step.
double* ClawPatch::q()
{
    return m_griddata.dataPtr();
}

FArrayBox ClawPatch::newGrid()
{
    /* Create a grid based on the size of the existing grids */
    Box b = m_griddata.box();
    int fields = m_griddata.fields();
    FArrayBox A;
    A.define(b,fields);
    return A;
}

Box ClawPatch::dataBox()
{
    return m_griddata.box();
}

Box ClawPatch::areaBox()
{
    return m_area.box();
}

Box ClawPatch::edgeBox()
{
    return m_edge_lengths.box();
}

Box ClawPatch::nodeBox()
{
    return m_xp.box();
}


/* Return a pointer to either time interpolated data or regular grid data */
double* ClawPatch::q_time_sync(fclaw_bool time_interp)
{
    if (time_interp)
        return m_griddata_time_interpolated.dataPtr();
    else
        return m_griddata.dataPtr();
}

void ClawPatch::save_current_step()
{
    m_griddata_last = m_griddata; // Copy for time interpolation
}

double ClawPatch::mx()
{
    return m_mx;
}
double ClawPatch::my()
{
    return m_my;
}
double ClawPatch::mbc()
{
    return m_mbc;
}
double ClawPatch::meqn()
{
    return m_meqn;
}

double ClawPatch::dx()
{
    return m_dx;
}

double ClawPatch::dy()
{
    return m_dy;
}

double ClawPatch::xlower()
{
    return m_xlower;
}

double ClawPatch::ylower()
{
    return m_ylower;
}

double ClawPatch::xupper()
{
    return m_xupper;
}

double ClawPatch::yupper()
{
    return m_yupper;
}

double* ClawPatch::xp()
{
    return m_xp.dataPtr();
}

double* ClawPatch::yp()
{
    return m_yp.dataPtr();
}
double* ClawPatch::zp()
{
    return m_zp.dataPtr();
}

double* ClawPatch::xd()
{
    return m_xd.dataPtr();
}
double* ClawPatch::yd()
{
    return m_yd.dataPtr();
}
double* ClawPatch::zd()
{
    return m_zd.dataPtr();
}

double* ClawPatch::area()
{
    return m_area.dataPtr();
}

double* ClawPatch::xface_normals()
{
    return m_xface_normals.dataPtr();
}

double* ClawPatch::yface_normals()
{
    return m_yface_normals.dataPtr();
}

double* ClawPatch::xface_tangents()
{
    return m_xface_tangents.dataPtr();
}

double* ClawPatch::yface_tangents()
{
    return m_yface_tangents.dataPtr();
}

double* ClawPatch::surf_normals()
{
    return m_surf_normals.dataPtr();
}

double* ClawPatch::edge_lengths()
{
    return m_edge_lengths.dataPtr();
}

double* ClawPatch::curvature()
{
    return m_curvature.dataPtr();
}

int* ClawPatch::block_corner_count()
{
    return &m_block_corner_count[0];
}

/* ----------------------------------------------------
   Solver data and functions
   ---------------------------------------------------*/
// Wave propagation algorithms

void* ClawPatch::clawpack_patch_data(int id)
{
    return fclaw_package_get_data(m_package_data_ptr,id);
    // return m_clawpack_patch_data;
}

/* ----------------------------------------------------------------
   Time stepping routines
   ---------------------------------------------------------------- */


void ClawPatch::save_step()
{
    // Store a backup in case the CFL number is too large doesn't work out.
    m_griddata_save = m_griddata;
}

void ClawPatch::restore_step()
{
    m_griddata = m_griddata_save;
}

void ClawPatch::ghost_comm(double *qpack, int time_interp,
                           int packmode)
{
    int ierror;

    // Number of internal layers.  At least four are needed to
    // average fine grid ghost patches onto coarse grid (on-proc)
    // ghost cells.
    int mint = ClawPatch::pack_layers;

    int packarea = packmode/2;   // (0,1)/2 = 0;  (2,3)/2 = 1;

    /* This is computed twice - here, and in fclaw2d_clawpatch_ghost_packsize */
    int wg = (2 + m_mx)*(2 + m_my);  // Whole grid (one layer of ghost cells)
    int hole = (m_mx - 2*mint)*(m_my - 2*mint);  // Hole in center
    FCLAW_ASSERT(hole >= 0);

    int psize = (wg - hole)*(m_meqn + packarea);
    FCLAW_ASSERT(psize > 0);

    double *q = q_time_sync(time_interp);
    double *area = m_area.dataPtr();  // Might be NULL

    FCLAW2D_GHOST_PACK(&m_mx,&m_my,&m_mbc,&m_meqn,&mint,q,area,
                       qpack,&psize,&packmode,&mint,&ierror);

    if (ierror > 0)
    {
        fclaw_global_essentialf("ghost_pack (ClawPatch.cpp) : ierror = %d\n",ierror);
        exit(0);
    }
}

void ClawPatch::partition_pack(double *q)
{
    m_griddata.copyToMemory(q);
}

void ClawPatch::partition_unpack(double *q)
{
    /* We only copy to griddata, since all levels are time
       synchronized  when repartitioning */
    m_griddata.copyFromMemory(q);
}

/* ----------------------------------------------------------------
   Exchange/average/interpolate
   ---------------------------------------------------------------- */

void ClawPatch::set_boundary_to_value(const int& time_interp,
                                      double& value)
{
    double *q = q_time_sync(time_interp);
    FCLAW2D_SET_BOUNDARY_TO_VALUE(&m_mx, &m_my, &m_mbc, &m_meqn, q, &value);
}

void ClawPatch::set_corners_to_value(const int& time_interp,
                                     double& value)
{
    double *q = q_time_sync(time_interp);
    FCLAW2D_SET_CORNERS_TO_VALUE(&m_mx, &m_my, &m_mbc, &m_meqn, q, &value);
}


void ClawPatch::exchange_face_ghost(const int& a_iface,
                                    ClawPatch *neighbor_cp,
                                    int time_interp,
                                    fclaw2d_transform_data_t* transform_data)
{
    double *qthis = q_time_sync(time_interp);
    double *qneighbor = neighbor_cp->q_time_sync(time_interp);
    exchange_face_ghost_(m_mx,m_my,m_mbc,m_meqn,qthis,qneighbor,a_iface,
                         &transform_data);
}


void ClawPatch::exchange_corner_ghost(const int& a_corner, ClawPatch *cp_corner,
                                      int time_interp,
                                      fclaw2d_transform_data_t* transform_data)
{
    double *qthis = q_time_sync(time_interp);
    double *qcorner = cp_corner->q_time_sync(time_interp);

    exchange_corner_ghost_(m_mx, m_my, m_mbc, m_meqn, qthis, qcorner, a_corner,
                           &transform_data);

}

/* ----------------------------------------------------------------
   Multi-level operations
   ---------------------------------------------------------------- */
void ClawPatch::average_face_ghost(const int& a_idir,
                                   const int& a_iface_coarse,
                                   const int& a_p4est_refineFactor,
                                   const int& a_refratio,
                                   ClawPatch *neighbor_cp,
                                   fclaw_bool a_time_interp,
                                   const int& igrid,
                                   fclaw2d_transform_data_t* transform_data)
{
    double *qcoarse = q_time_sync(a_time_interp);

    double *qfine = neighbor_cp->m_griddata.dataPtr();

    /* These will be empty for non-manifolds cases */
    double *areacoarse = m_area.dataPtr();
    double *areafine = neighbor_cp->m_area.dataPtr();

    int manifold = m_manifold ? 1 : 0;
    average_face_ghost_(m_mx,m_my,m_mbc,m_meqn,
                        qcoarse,qfine,
                        areacoarse, areafine,
                        a_idir,a_iface_coarse,
                        a_p4est_refineFactor,a_refratio,igrid,
                        manifold, &transform_data);
}

void ClawPatch::interpolate_face_ghost(const int& a_idir,
                                       const int& a_iside,
                                       const int& a_p4est_refineFactor,
                                       const int& a_refratio,
                                       ClawPatch *neighbor_cp,
                                       fclaw_bool a_time_interp,
                                       const int& igrid,
                                       fclaw2d_transform_data_t* transform_data)
{
    double *qcoarse = q_time_sync(a_time_interp);
    double *qfine = neighbor_cp->m_griddata.dataPtr();

    interpolate_face_ghost_(m_mx,m_my,m_mbc,m_meqn,qcoarse,qfine,a_idir,a_iside,
                            a_p4est_refineFactor,a_refratio,igrid,
                            &transform_data);
}

//
void ClawPatch::average_corner_ghost(const int& a_coarse_corner,
                                     const int& a_refratio,
                                     ClawPatch *cp_corner,
                                     fclaw_bool a_time_interp,
                                     fclaw2d_transform_data_t* transform_data)
{

    double *qcoarse = q_time_sync(a_time_interp);
    double *qfine = cp_corner->m_griddata.dataPtr();
    double *areacoarse = this->m_area.dataPtr();
    double *areafine = cp_corner->m_area.dataPtr();

    int manifold = m_manifold ? 1 : 0;
    average_corner_ghost_(m_mx, m_my, m_mbc, m_meqn, a_refratio,
                          qcoarse, qfine,
                          areacoarse, areafine,
                          manifold,
                          a_coarse_corner,&transform_data);
}



void ClawPatch::interpolate_corner_ghost(const int& a_coarse_corner,
                                         const int& a_refratio,
                                         ClawPatch *cp_corner,
                                         fclaw_bool a_time_interp,
                                         fclaw2d_transform_data_t* transform_data)

{
    double *qcoarse = q_time_sync(a_time_interp);

    // qcorner is the finer level.
    double *qfine = cp_corner->m_griddata.dataPtr();

    interpolate_corner_ghost_(m_mx, m_my, m_mbc, m_meqn,
                              a_refratio, qcoarse, qfine,
                              a_coarse_corner,&transform_data);
}


/* ----------------------------------------------------------------
   Special case : Pillow grid ghost exchanges/average/interpolate
   ---------------------------------------------------------------- */

void ClawPatch::mb_exchange_block_corner_ghost(const int& a_corner,
                                               ClawPatch *cp_corner,
                                               int time_interp)
{
    double *qthis = q_time_sync(time_interp);
    double *qcorner = cp_corner->m_griddata.dataPtr();


    mb_exchange_block_corner_ghost_(m_mx, m_my, m_mbc, m_meqn, qthis, qcorner,
                                    a_corner, m_blockno);

}

// internal corners only a block boundaries.
void ClawPatch::mb_average_block_corner_ghost(const int& a_coarse_corner,
                                              const int& a_refratio,
                                              ClawPatch *cp_corner,
                                              fclaw_bool a_time_interp)
{
    // 'this' is the finer grid; 'cp_corner' is the coarser grid.
    double *qcoarse = q_time_sync(a_time_interp);


    double *areacoarse = this->m_area.dataPtr();
    double *areafine = cp_corner->m_area.dataPtr();
    double *qfine = cp_corner->m_griddata.dataPtr();

    mb_average_block_corner_ghost_(m_mx,m_my,m_mbc,m_meqn,
                                   a_refratio,qcoarse,qfine,
                                   areacoarse,areafine,
                                   a_coarse_corner,m_blockno);
}


void ClawPatch::mb_interpolate_block_corner_ghost(const int& a_coarse_corner,
                                                  const int& a_refratio,
                                                  ClawPatch *cp_corner,
                                                  fclaw_bool a_time_interp)

{
    double *qcoarse = q_time_sync(a_time_interp);

    /* qcorner is the finer level. */
    double *qfine = cp_corner->m_griddata.dataPtr();

    mb_interpolate_block_corner_ghost_(m_mx, m_my, m_mbc, m_meqn,
                                       a_refratio, qcoarse, qfine,
                                       a_coarse_corner, m_blockno);
}



/* ----------------------------------------------------------------
   Time interpolation
   ---------------------------------------------------------------- */
void ClawPatch::setup_for_time_interpolation(const double& alpha,
                                             const int& psize)
{
    /* Store time interpolated data that will be use in coarse grid
       exchanges */
    double *qlast = m_griddata_last.dataPtr();
    double *qcurr = m_griddata.dataPtr();
    double *qinterp = m_griddata_time_interpolated.dataPtr();
    int ierror;

    /* Do interpolation only on interior, since ghost cells in qcurr
       are invalid and will lead to floating point exceptions */
    FCLAW2D_TIMEINTERP_FORT(&m_mx,&m_my,&m_mbc,&m_meqn,&psize,
                            qcurr,qlast,qinterp,&alpha,&ierror);
    if (ierror > 0)
    {
        fclaw_global_essentialf("fclaw2d_timeinterp : incorrect pack size\n");
        exit(0);
    }
}

void ClawPatch::finegrid_neighbors(int y)
{
    m_has_finegrid_neighbors = y;
}

int ClawPatch::has_finegrid_neighbors()
{
    return m_has_finegrid_neighbors;
}


/* ----------------------------------------------------------------
   Mapped grids
   ---------------------------------------------------------------- */

void ClawPatch::set_block_corner_count(const int icorner, const int block_corner_count)
{
    m_block_corner_count[icorner] = block_corner_count;
}

void ClawPatch::setup_area_storage()
{
    int mx = m_mx;
    int my = m_my;
    int mbc = m_mbc;

    int ll[SpaceDim];
    int ur[SpaceDim];
    for (int idir = 0; idir < SpaceDim; idir++)
    {
        ll[idir] = -mbc;
    }
    ur[0] = mx + mbc + 1;
    ur[1] = my + mbc + 1;

    Box box_p(ll,ur);
    m_area.define(box_p,1);
}

void ClawPatch::setup_metric_storage()
{
    int mx = m_mx;
    int my = m_my;
    int mbc = m_mbc;

    int ll[SpaceDim];
    int ur[SpaceDim];
    for (int idir = 0; idir < SpaceDim; idir++)
    {
        ll[idir] = -mbc;
    }
    ur[0] = mx + mbc + 1;
    ur[1] = my + mbc + 1;

    Box box_p(ll,ur);   /* Store cell centered values here */

    /* Mesh cell centers of physical mesh */
    m_xp.define(box_p,1);
    m_yp.define(box_p,1);
    m_zp.define(box_p,1);
    m_surf_normals.define(box_p,3);
    m_curvature.define(box_p,3);

    /* Node centered values */
    for (int idir = 0; idir < SpaceDim; idir++)
    {
        ll[idir] = -mbc;
    }
    ur[0] = mx + mbc + 2;
    ur[1] = my + mbc + 2;
    Box box_d(ll,ur);

    m_xd.define(box_d,1);
    m_yd.define(box_d,1);
    m_zd.define(box_d,1);

    /* Face centered values */
    m_xface_normals.define(box_d,3);
    m_yface_normals.define(box_d,3);
    m_xface_tangents.define(box_d,3);
    m_yface_tangents.define(box_d,3);
    m_edge_lengths.define(box_d,2);
}

/* ----------------------------------------------------------------
   Output and diagnostics
   ---------------------------------------------------------------- */

int ClawPatch::size()
{
    /* Use this to create new data */
    return m_griddata.size();
}
