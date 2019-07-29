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

#include <fclaw2d_clawpatch_transform.h>

#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_patch.h>
#include <fclaw2d_global.h>

void fclaw2d_clawpatch_transform_init_data(fclaw2d_global_t* glob, 
                                           fclaw2d_patch_t* this_patch,
                                           int blockno, int patchno,
                                           fclaw2d_patch_transform_data_t* transform)
{
    /* Cell centered data */
    transform->based = 1;

    /* Nothing to do for transform->user */
}

void fclaw2d_clawpatch_face_transformation (int faceno, int rfaceno, int ftransform[])
{
    /* Defined in forestclaw2d.c */
    fclaw2d_patch_face_transformation (faceno, rfaceno, ftransform);
}

void fclaw2d_clawpatch_face_transformation_intra (int ftransform[])
{
    /* Defined in forestclaw2d.c */
    fclaw2d_patch_face_transformation_intra (ftransform);
}


/* Same size neighbor across a face */
void
FCLAW2D_CLAWPATCH_TRANSFORM_FACE (const int *i1, const int *j1,
                        int *i2, int *j2, fclaw2d_patch_transform_data_t** ptdata)
{
    fclaw2d_patch_transform_data_t *tdata = *ptdata;
    const fclaw2d_clawpatch_options_t *clawpatch_opt = 
             fclaw2d_clawpatch_get_options(tdata->glob);

    *i2 = *i1;
    *j2 = *j1;
    fclaw2d_patch_transform_face (tdata->this_patch,
                                  tdata->neighbor_patch,
                                  tdata->transform,
                                  clawpatch_opt->mx, 
                                  clawpatch_opt->my, 
                                  tdata->based, i2, j2);
}


/* Half size neighbor across a face */
void
FCLAW2D_CLAWPATCH_TRANSFORM_FACE_HALF (const int *i1, const int *j1,
                                       int i2[], int j2[],
                                       fclaw2d_patch_transform_data_t** ptdata)
{
    fclaw2d_patch_transform_data_t *tdata = *ptdata;
    const fclaw2d_clawpatch_options_t *clawpatch_opt = 
                fclaw2d_clawpatch_get_options(tdata->glob);

    i2[0] = *i1;
    j2[0] = *j1;
    fclaw2d_patch_transform_face2 (tdata->this_patch,
                                   tdata->neighbor_patch,
                                   tdata->transform, 
                                   clawpatch_opt->mx, 
                                   clawpatch_opt->my,
                                   tdata->based, i2, j2);
}


/* TODO: Extend this for a block-block corner */
void
FCLAW2D_CLAWPATCH_TRANSFORM_CORNER (const int *i1, const int *j1,
                                    int *i2, int *j2,
                                    fclaw2d_patch_transform_data_t** ptdata)
{
    fclaw2d_patch_transform_data_t *tdata = *ptdata;
    const fclaw2d_clawpatch_options_t *clawpatch_opt = 
                   fclaw2d_clawpatch_get_options(tdata->glob);

    *i2 = *i1;
    *j2 = *j1;
    if (tdata->block_iface >= 0)
    {
        /* block-face but not a block-corner */
#if 0
        FCLAW_ASSERT (tdata->block_iface < 4);
#endif
        fclaw2d_patch_transform_face (tdata->this_patch,
                                      tdata->neighbor_patch, tdata->transform,
                                      clawpatch_opt->mx, clawpatch_opt->my,
                                      tdata->based, i2, j2);
    }
    else
    {
        /* corner within a block */
        FCLAW_ASSERT (tdata->block_iface == -1);
        fclaw2d_patch_transform_corner (tdata->this_patch,
                                        tdata->neighbor_patch,
                                        tdata->icorner, tdata->is_block_corner,
                                        clawpatch_opt->mx, clawpatch_opt->my,
                                        tdata->based, i2, j2);
    }
    /* Done. */
    /* TODO: We need to permit that it's a block corner.  In this case,
     * call fclaw2d_patch_transform_corner with is_block_boundary = 1 */
}

/* TODO: Extend this for a block-block corner */
void
FCLAW2D_CLAWPATCH_TRANSFORM_CORNER_HALF (const int *i1, const int *j1,
                                         int *i2, int *j2,
                                         fclaw2d_patch_transform_data_t** ptdata)
{
    fclaw2d_patch_transform_data_t *tdata = *ptdata;
    const fclaw2d_clawpatch_options_t *clawpatch_opt = 
               fclaw2d_clawpatch_get_options(tdata->glob);

    i2[0] = *i1;
    j2[0] = *j1;
    if (tdata->block_iface >= 0)
    {
        /* block-face but not a block-corner. */
        fclaw2d_patch_transform_face2 (tdata->this_patch,
                                       tdata->neighbor_patch,
                                       tdata->transform, 
                                       clawpatch_opt->mx, 
                                       clawpatch_opt->my,
                                       tdata->based, i2, j2);
    }
    else
    {
        /* corner within a block */
        FCLAW_ASSERT (tdata->block_iface == -1);
        fclaw2d_patch_transform_corner2 (tdata->this_patch,
                                         tdata->neighbor_patch,
                                         tdata->icorner, tdata->is_block_corner,
                                         clawpatch_opt->mx, clawpatch_opt->my,
                                         tdata->based, i2, j2);
    }
    /* Done */
    /* TODO: We need to permit that it's a block corner.  In this case,
     * call fclaw2d_patch_transform_corner2 with is_block_boundary = 1 */
}
