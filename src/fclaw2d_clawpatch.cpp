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


#include <forestclaw2d.h>
#include "fclaw2d_defs.H"
#include "fclaw_options.h"
#include <fclaw_package.h>
#include <fclaw2d_vtable.h>
#include <amr_utils.h>

#include <fclaw2d_clawpatch.h>

class ClawPatch;

void fclaw2d_clawpatch_grid_data(fclaw2d_domain_t* domain,
                                 fclaw2d_patch_t* this_patch,
                                 int* mx, int* my, int* mbc,
                                 double* xlower, double* ylower,
                                 double* dx, double* dy)
{
    ClawPatch *cp = get_clawpatch(this_patch);
    *mx = cp->mx();
    *my = cp->my();
    *mbc = cp->mbc();
    *xlower = cp->xlower();
    *ylower = cp->ylower();
    *dx = cp->dx();
    *dy = cp->dy();
}

void fclaw2d_clawpatch_metric_data(fclaw2d_domain_t* domain,
                                   fclaw2d_patch_t* this_patch,
                                   double **xp, double **yp, double **zp,
                                   double **xd, double **yd, double **zd,
                                   double **area)
{
    /* Declared as a friend in ClawPatch.H */
    ClawPatch *cp = get_clawpatch(this_patch);
    *xp = cp->xp();
    *yp = cp->yp();
    *zp = cp->zp();
    *xd = cp->xd();
    *yd = cp->yd();
    *zd = cp->zd();
    *area = cp->area();
}

double* fclaw2d_clawpatch_get_area(fclaw2d_domain_t* domain,
                                   fclaw2d_patch_t* this_patch)
{
    ClawPatch *cp = get_clawpatch(this_patch);
    return cp->area();
}

double *fclaw2d_clawpatch_get_q(fclaw2d_domain_t* domain,
                                fclaw2d_patch_t* this_patch)
{
    ClawPatch *cp = get_clawpatch(this_patch);
    return cp->q();
}


double *fclaw2d_clawpatch_get_q_time_interp(fclaw2d_domain_t* domain,
                                            fclaw2d_patch_t* this_patch)
{
    ClawPatch *cp = get_clawpatch(this_patch);
    return cp->q_time_interp();
}


void fclaw2d_clawpatch_metric_data2(fclaw2d_domain_t* domain, fclaw2d_patch_t* this_patch,
                                    double **xnormals, double **ynormals,
                                    double **xtangents, double **ytangents,
                                    double **surfnormals,
                                    double **edgelengths, double **curvature)
{
    /* or just call the member functions? */
    ClawPatch *cp = get_clawpatch(this_patch);
    *xnormals    = cp->xface_normals();
    *ynormals    = cp->yface_normals();
    *xtangents   = cp->xface_tangents();
    *ytangents   = cp->yface_tangents();
    *surfnormals = cp->surf_normals();
    *curvature   = cp->curvature();
    *edgelengths = cp->edge_lengths();
}

void fclaw2d_clawpatch_soln_data(fclaw2d_domain_t* domain,
                                 fclaw2d_patch_t* this_patch,
                                 double **q, int* meqn)
{
    ClawPatch *cp = get_clawpatch(this_patch);
    *q = cp->q();
    *meqn = cp->meqn();
}

void fclaw2d_clawpatch_timesync_data(fclaw2d_domain_t* domain,
                                     fclaw2d_patch_t* this_patch,
                                     fclaw_bool time_interp,
                                     double **q, int* meqn)
{
    ClawPatch *cp = get_clawpatch(this_patch);
    *q = cp->q_time_sync(time_interp);
    *meqn = cp->meqn();
}


void fclaw2d_clawpatch_save_current_step(fclaw2d_domain_t* domain,
                                         fclaw2d_patch_t* this_patch)
{
    ClawPatch *cp = get_clawpatch(this_patch);
    cp->save_current_step();
}

int* fclaw2d_clawpatch_corner_count(fclaw2d_domain_t* domain,
                                   fclaw2d_patch_t* this_patch)
{
    ClawPatch *cp = get_clawpatch(this_patch);
    return cp->block_corner_count();
}

/* ------------------------------------------------------------------
   Repartition and rebuild new domains, or construct initial domain
 -------------------------------------------------------------------- */
#ifdef __cplusplus
extern "C" {
#if 0
}
#endif
#endif


void fclaw2d_clawpatch_define(fclaw2d_domain_t* domain,
                              fclaw2d_patch_t *this_patch,
                              int blockno, int patchno)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    int level = this_patch->level;

    ClawPatch *cp = new ClawPatch();
    cp->define(this_patch->xlower,
               this_patch->ylower,
               this_patch->xupper,
               this_patch->yupper,
               blockno,
               level,
               gparms);

    fclaw2d_patch_data_t *pdata = get_patch_data(this_patch);
    pdata->cp = cp;

    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    ++ddata->count_set_clawpatch;
}

void fclaw2d_clawpatch_build_cb(fclaw2d_domain_t *domain,
                                fclaw2d_patch_t *this_patch,
                                int this_block_idx,
                                int this_patch_idx,
                                void *user)
{
    fclaw2d_vtable_t vt;
    vt = fclaw2d_get_vtable(domain);

    fclaw2d_clawpatch_define(domain,this_patch,this_block_idx,this_patch_idx);

    if (vt.patch_setup != NULL)
    {
        vt.patch_setup(domain,this_patch,this_block_idx,this_patch_idx);
    }
}

void fclaw2d_clawpatch_pack_cb(fclaw2d_domain_t *domain,
                               fclaw2d_patch_t *this_patch,
                               int this_block_idx,
                               int this_patch_idx,
                               void *user)
{
    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    int patch_num = this_block->num_patches_before + this_patch_idx;
    double* patch_data = (double*) ((void**)user)[patch_num];

    pack_clawpatch(this_patch,patch_data);
}

void fclaw2d_clawpatch_unpack_cb(fclaw2d_domain_t *domain,
                                 fclaw2d_patch_t *this_patch,
                                 int this_block_idx,
                                 int this_patch_idx,
                                 void *user)
{
    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    int patch_num = this_block->num_patches_before + this_patch_idx;
    double* patch_data = (double*) ((void**)user)[patch_num];

    fclaw_bool time_interp = fclaw_false;
    unpack_clawpatch(domain,this_patch,this_block_idx,this_patch_idx,
                     patch_data,time_interp);
}

void delete_clawpatch(fclaw2d_domain_t* domain, fclaw2d_patch_t* this_patch)
{
    // We expect this ClawPatch to exist.
    fclaw2d_patch_data_t *pdata = get_patch_data(this_patch);
    delete pdata->cp;
    pdata->cp = NULL;

    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    ++ddata->count_delete_clawpatch;
}

void pack_clawpatch(fclaw2d_patch_t* this_patch,double* qdata)
{
    ClawPatch *cp = get_clawpatch(this_patch);
    cp->pack_griddata(qdata);
}

void unpack_clawpatch(fclaw2d_domain_t* domain, fclaw2d_patch_t* this_patch,
                      int this_block_idx, int this_patch_idx, double *qdata,
                      fclaw_bool time_interp)
{
    ClawPatch *cp = get_clawpatch(this_patch);
    if (time_interp)
        cp->unpack_griddata_time_interpolated(qdata);
    else
        cp->unpack_griddata(qdata);
}

size_t fclaw2d_clawpatch_pack_size(fclaw2d_domain_t* domain)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;
    size_t size = (2*mbc + mx)*(2*mbc+my)*meqn;
    return size*sizeof(double);
}


#ifdef __cplusplus
#if 0
{
#endif
}
#endif
