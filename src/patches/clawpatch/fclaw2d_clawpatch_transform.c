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

#ifndef REFINE_DIM
#define REFINE_DIM 2
#endif

#ifndef PATCH_DIM
#define PATCH_DIM 2
#endif

#include <forestclaw.h>
#include <fclaw_patch.h>
#include <fclaw_global.h>

#include <fclaw2d_clawpatch_transform.h>

#include <fclaw_clawpatch_options.h>

void fclaw2d_clawpatch_transform_init_data(fclaw_global_t* glob, 
                                           fclaw_patch_t* this_patch,
                                           int blockno, int patchno,
                                           fclaw_patch_transform_data_t* transform)
{
    /* Cell centered data */
    transform->based = 1;

    /* Store clawpatch_options in transform->user */
    transform->user = (void*) fclaw_clawpatch_get_options(glob);
}

void fclaw2d_clawpatch_face_transformation (int faceno, int rfaceno, int ftransform[])
{
    /* Defined in forestclaw2d.c */
    fclaw_patch_face_transformation (2, faceno, rfaceno, ftransform);
}

void fclaw2d_clawpatch_face_transformation_intra (int ftransform[])
{
    /* Defined in forestclaw2d.c */
    fclaw_patch_face_transformation_intra (2, ftransform);
}


/* Same size neighbor across a face */
void
FCLAW2D_CLAWPATCH_TRANSFORM_FACE (const int *i1, const int *j1,
                        int *i2, int *j2, fclaw_patch_transform_data_t** ptdata)
{
    fclaw_patch_transform_data_t *tdata = *ptdata;
    const fclaw_clawpatch_options_t *clawpatch_opt = 
        (fclaw_clawpatch_options_t*) tdata->user;

    int mx,my;
    if(clawpatch_opt->patch_dim == 2)
    {
        mx = clawpatch_opt->mx;
        my = clawpatch_opt->my;
    }
    else
    {
        mx = clawpatch_opt->mx;
        my = clawpatch_opt->my;
    }

    *i2 = *i1;
    *j2 = *j1;
    fclaw_patch_2d_transform_face (tdata->this_patch,
                                   tdata->neighbor_patch,
                                   tdata->transform,
                                   mx, my, 
                                   tdata->based, i2, j2);
}


/* Half size neighbor across a face */
void
FCLAW2D_CLAWPATCH_TRANSFORM_FACE_HALF (const int *i1, const int *j1,
                                       int i2[], int j2[],
                                       fclaw_patch_transform_data_t** ptdata)
{
    fclaw_patch_transform_data_t *tdata = *ptdata;
    const fclaw_clawpatch_options_t *clawpatch_opt = 
        (fclaw_clawpatch_options_t*) tdata->user;

    int mx,my;
    if(clawpatch_opt->patch_dim == 2)
    {
        mx = clawpatch_opt->mx;
        my = clawpatch_opt->my;
    }
    else
    {
        mx = clawpatch_opt->mx;
        my = clawpatch_opt->my;
    }

    i2[0] = *i1;
    j2[0] = *j1;
    fclaw_patch_2d_transform_face2 (tdata->this_patch,
                                    tdata->neighbor_patch,
                                    tdata->transform, 
                                    mx, my,
                                    tdata->based, i2, j2);
}


void
FCLAW2D_CLAWPATCH_TRANSFORM_CORNER (const int *i1, const int *j1,
                                    int *i2, int *j2,
                                    fclaw_patch_transform_data_t** ptdata)
{
    fclaw_patch_transform_data_t *tdata = *ptdata;
    const fclaw_clawpatch_options_t *clawpatch_opt = 
        (fclaw_clawpatch_options_t*) tdata->user;

    int mx,my;
    if(clawpatch_opt->patch_dim == 2)
    {
        mx = clawpatch_opt->mx;
        my = clawpatch_opt->my;
    }
    else
    {
        mx = clawpatch_opt->mx;
        my = clawpatch_opt->my;
    }

    *i2 = *i1;
    *j2 = *j1;
    if (tdata->block_iface >= 0)
    {
        /* block-face but not a block-corner */
#if 0
        FCLAW_ASSERT (tdata->block_iface < 4);
#endif
        fclaw_patch_2d_transform_face (tdata->this_patch,
                                       tdata->neighbor_patch, tdata->transform,
                                       mx, my,
                                       tdata->based, i2, j2);
    }
    else
    {
        /* Corner within a block or a block-block corner. For block-block
         * corners, we assume both patches lie in coordinate systems with the
         * same orientation. */
        FCLAW_ASSERT (tdata->block_iface == -1);
        fclaw_patch_2d_transform_corner (tdata->this_patch,
                                         tdata->neighbor_patch,
                                         tdata->icorner, tdata->is_block_corner,
                                         mx, my,
                                         tdata->based, i2, j2);
    }
}

void
FCLAW2D_CLAWPATCH_TRANSFORM_CORNER_HALF (const int *i1, const int *j1,
                                         int *i2, int *j2,
                                         fclaw_patch_transform_data_t** ptdata)
{
    fclaw_patch_transform_data_t *tdata = *ptdata;
    const fclaw_clawpatch_options_t *clawpatch_opt = 
        (fclaw_clawpatch_options_t*) tdata->user;

    int mx,my;
    if(clawpatch_opt->patch_dim == 2)
    {
        mx = clawpatch_opt->mx;
        my = clawpatch_opt->my;
    }
    else
    {
        mx = clawpatch_opt->mx;
        my = clawpatch_opt->my;
    }

    i2[0] = *i1;
    j2[0] = *j1;
    if (tdata->block_iface >= 0)
    {
        /* block-face but not a block-corner. */
        fclaw_patch_2d_transform_face2 (tdata->this_patch,
                                        tdata->neighbor_patch,
                                        tdata->transform, 
                                        mx, my,
                                        tdata->based, i2, j2);
    }
    else
    {
        /* Corner within a block or a block-block corner. For block-block
         * corners, we assume both patches lie in coordinate systems with the
         * same orientation. */
        FCLAW_ASSERT (tdata->block_iface == -1);
        fclaw_patch_2d_transform_corner2 (tdata->this_patch,
                                          tdata->neighbor_patch,
                                          tdata->icorner, tdata->is_block_corner,
                                          mx, my,
                                          tdata->based, i2, j2);
    }
}
