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

#include <fclaw2d_transform.h>

/* Same size neighbor across a face */
void
FCLAW2D_TRANSFORM_FACE (const int *i1, const int *j1,
                        int *i2, int *j2, fclaw2d_transform_data_t** ptdata)
{
    fclaw2d_transform_data_t *tdata = *ptdata;
    *i2 = *i1;
    *j2 = *j1;
    fclaw2d_patch_transform_face (tdata->this_patch,
                                  tdata->neighbor_patch,
                                  tdata->transform,
                                  tdata->mx, tdata->my, tdata->based, i2, j2);
}

/* This works except for a block-block corner */
void
FCLAW2D_TRANSFORM_CORNER (const int *i1, const int *j1,
                          int *i2, int *j2,
                          fclaw2d_transform_data_t** ptdata)
{
    fclaw2d_transform_data_t *tdata = *ptdata;
    *i2 = *i1;
    *j2 = *j1;
    if (tdata->iface >= 0)
    {
        /* block-face but no block-corner */
        FCLAW_ASSERT (tdata->iface < 4);
        fclaw2d_patch_transform_face (tdata->this_patch,
                                      tdata->neighbor_patch, tdata->transform,
                                      tdata->mx, tdata->my,
                                      tdata->based, i2, j2);
    }
    else
    {
        /* corner within a block */
        FCLAW_ASSERT (tdata->iface == -1);
        fclaw2d_patch_transform_corner (tdata->this_patch,
                                        tdata->neighbor_patch, tdata->icorner,
                                        tdata->mx, tdata->my,
                                        tdata->based, i2, j2);
    }
}

/* Half size neighbor across a face */
void
FCLAW2D_TRANSFORM_FACE_HALF (const int *i1, const int *j1,
                             int i2[], int j2[],
                             fclaw2d_transform_data_t** ptdata)
{
    fclaw2d_transform_data_t *tdata = *ptdata;
    i2[0] = *i1;
    j2[0] = *j1;
    fclaw2d_patch_transform_face2 (tdata->this_patch,
                                   tdata->neighbor_patch,
                                   tdata->transform, tdata->mx, tdata->my,
                                   tdata->based, i2, j2);
}

/* This works except for a block-block corner */
void
FCLAW2D_TRANSFORM_CORNER_HALF (const int *i1, const int *j1,
                               int *i2, int *j2,
                               fclaw2d_transform_data_t** ptdata)
{
    fclaw2d_transform_data_t *tdata = *ptdata;
    i2[0] = *i1;
    j2[0] = *j1;
    if (tdata->iface >= 0)
    {
        /* block-face but no block-corner */
        FCLAW_ASSERT (tdata->iface < 4);
        fclaw2d_patch_transform_face2 (tdata->this_patch,
                                       tdata->neighbor_patch,
                                       tdata->transform, tdata->mx, tdata->my,
                                       tdata->based, i2, j2);
    }
    else
    {
        /* corner within a block */
        FCLAW_ASSERT (tdata->iface == -1);
        fclaw2d_patch_transform_corner2 (tdata->this_patch,
                                         tdata->neighbor_patch,
                                         tdata->icorner, tdata->mx, tdata->my,
                                         tdata->based, i2, j2);
    }
}

/* This obviously is bogus - but I have here it so I can get the
   rest of the code working.
   TODO: remove this function */
void
transform_corner_halfsize2_ (const int *i1, const int *j1,
                             int *i2, int *j2,
                             fclaw2d_transform_data_t** ptdata)
{


    fclaw2d_transform_data_t *tdata = *ptdata;

    printf("halfsize2 : we should not be calling this function\n");
    exit(0);

    int mbc = 2;
    int icorner = tdata->icorner;

    int ibc, jbc;

    for (jbc = 1; jbc <= mbc; jbc++)
    {
        for (ibc = 1; ibc <= mbc; ibc++)
        {
            /* These transformations don't depend on i1,j1 ... */
            if (icorner == 0)
            {
                *i2++ = tdata->mx + ibc;
                *j2++ = tdata->my + jbc;
            }
            else if (icorner == 1)
            {
                *i2++ = 1 - ibc;
                *j2++ = tdata->my + jbc;
            }
            else if (icorner == 2)
            {
                *i2++ = tdata->mx + ibc;
                *j2++ = 1 - jbc;
            }
            else if (icorner == 3)
            {
                *i2++ = 1 - ibc;
                *j2++ = 1 - jbc;
            }
        }
    }
}
