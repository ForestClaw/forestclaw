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
#include <fclaw2d_manifold_default.h>
#include <fclaw2d_manifold.h>

#include <ClawPatch.hpp>

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

void fclaw2d_clawpatch_setup_timeinterp(fclaw2d_domain_t* domain,
                                        fclaw2d_patch_t *this_patch,
                                        double alpha)
{
    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    cp->setup_for_time_interpolation(alpha);
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
#if 0
void fclaw2d_clawpatch_area_storage(fclaw2d_domain_t* domain,
                                    fclaw2d_patch_t* this_patch,
                                    int blockno,
                                    int patchno)
{
    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    cp->setup_area_storage();
}

void fclaw2d_clawpatch_manifold_storage(fclaw2d_domain_t* domain,
                                        fclaw2d_patch_t* this_patch,
                                        int blockno,
                                        int patchno)
{
    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    cp->setup_manifold_storage();
}
#endif

void fclaw2d_clawpatch_manifold_setup(fclaw2d_domain_t* domain,
                                      fclaw2d_patch_t* this_patch,
                                      int blockno,
                                      int patchno)
{
    fclaw2d_vtable_t vt;
    vt = fclaw2d_get_vtable(domain);

    /* vt.patch_manifold_setup_mesh(...) */
    fclaw2d_manifold_setup_mesh(domain,this_patch,blockno,patchno);

    /* vt.patch_manifold_compute_normals(...) */
    fclaw2d_manifold_compute_normals(domain,this_patch,blockno,patchno);
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
    const amr_options_t *gparms = get_domain_parms(domain);
    int level = this_patch->level;

    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    cp->define(this_patch->xlower,
               this_patch->ylower,
               this_patch->xupper,
               this_patch->yupper,
               blockno,
               level,
               gparms,
               build_mode);

}


/* This is called from fclaw2d_regrid_new_domain_setup (in fclaw2d_regrid.cpp)
   every time a new domain is created after regridding.  This is followed by
   an "iterate_adapted" step, in which data in these patches is copied, averaged or
   interpolated. */

void fclaw2d_clawpatch_build_cb(fclaw2d_domain_t *domain,
                                fclaw2d_patch_t *this_patch,
                                int blockno,
                                int patchno,
                                void *user)
{
    fclaw2d_vtable_t vt;

    fclaw2d_build_mode_t build_mode =  *((fclaw2d_build_mode_t*) user);
    const amr_options_t *gparms = get_domain_parms(domain);

    fclaw2d_clawpatch_define(domain,this_patch,blockno,patchno,build_mode);

    if (gparms->manifold)
    {
        if (build_mode != FCLAW2D_BUILD_FOR_GHOST_AREA_PACKED)
        {
            fclaw2d_manifold_compute_area(domain,this_patch,blockno,patchno);

            /* Don't need any more manifold info for ghost patches */
            if (build_mode == FCLAW2D_BUILD_FOR_UPDATE)
            {
                fclaw2d_clawpatch_manifold_setup(domain,this_patch,blockno,patchno);
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
    const amr_options_t *gparms = get_domain_parms(domain);

    fclaw2d_clawpatch_define(domain,coarse_patch,blockno,coarse_patchno,build_mode);

    if (gparms->manifold)
    {
        fclaw2d_manifold_average_area(domain,fine_patches,coarse_patch,
                                      blockno, coarse_patchno, fine0_patchno);
        fclaw2d_clawpatch_manifold_setup(domain,coarse_patch,blockno,
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

void fclaw2d_clawpatch_partition_pack_cb(fclaw2d_domain_t *domain,
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

void fclaw2d_clawpatch_partition_unpack_cb(fclaw2d_domain_t *domain,
                                           fclaw2d_patch_t *this_patch,
                                           int this_block_idx,
                                           int this_patch_idx,
                                           void *user)
{
    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    int patch_num = this_block->num_patches_before + this_patch_idx;
    double* patch_data = (double*) ((void**)user)[patch_num];

    /* First need to rebuild the patch */
    fclaw2d_patch_data_new(domain,this_patch);

    fclaw2d_build_mode_t build_mode = FCLAW2D_BUILD_FOR_UPDATE;
    fclaw2d_clawpatch_build_cb(domain,this_patch,this_block_idx,
                               this_patch_idx,(void*) &build_mode);

    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);

    /* Time interp is false, since we only partition when all grids
       are time synchronized */
    cp->partition_unpack(patch_data);

}
