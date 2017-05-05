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
#include <fclaw2d_output.h>
#include <fclaw2d_vtk.h>
#include <fclaw2d_map.h>

#include <fclaw2d_clawpatch.h>

static void
fclaw2d_output_vtk_coordinate_cb (fclaw2d_global_t * glob,
                                  fclaw2d_patch_t * this_patch,
                                  int this_block_idx, int this_patch_idx,
                                  char *a)
{
    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(glob);
    fclaw2d_map_context_t *cont;

    int mx,my,mbc;
    double dx,dy,xlower,ylower;


    cont = glob->cont;

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

#if 0
    fclaw2d_clawpatch_metric_data(glob,this_patch,&xp,&yp,&zp,&xd,&yd,&zd,&area);

    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    const double xlower = cp->xlower ();
    const double ylower = cp->ylower ();
    const double dx = cp->dx ();
    const double dy = cp->dy ();
#endif

    /* Enumerate point coordinates in the patch */
    double *d = (double *) a;
    int i, j;
    double xpp,ypp,zpp;
    for (j = 0; j <= my; ++j)
    {
        const double y = ylower + j * dy;
        for (i = 0; i <= mx; ++i)
        {
            const double x = xlower + i * dx;
            if (gparms->manifold)
            {
                FCLAW2D_MAP_C2M(&cont,&this_block_idx,&x,&y,&xpp,&ypp,&zpp);
                *d++ = xpp;
                *d++ = ypp;
                *d++ = zpp;
            }
            else
            {
                *d++ = x;
                *d++ = y;
                *d++ = 0;
            }
        }
    }
}


static void
fclaw2d_output_vtk_value_cb (fclaw2d_global_t * glob,
                             fclaw2d_patch_t * this_patch,
                             int this_block_idx, int this_patch_idx,
                             char *a)
{
    double *q;
    double xlower,ylower,dx,dy;
    int mx,my,mbc,meqn;

    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    const int xlane = mx + 2 * mbc;
    const int ylane = my + 2 * mbc;

    // Enumerate equation data in the patch
    float *f = (float *) a;
    int i, j, k;
    for (j = 0; j < my; ++j)
    {
        for (i = 0; i < mx; ++i)
        {
            for (k = 0; k < meqn; ++k)
            {
                *f++ = (float) q[(k * ylane + j + mbc) * xlane + i + mbc];
            }
        }
    }
}

/*  ---------------------------------------------------------------------------
    Public interface
    --------------------------------------------------------------------------- */

void
fclaw2d_clawpatch_output_vtk (fclaw2d_global_t * glob, int iframe)
{
    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(glob);
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);


    char basename[BUFSIZ];
    snprintf (basename, BUFSIZ, "%s_frame_%04d", gparms->prefix, iframe);

    (void) fclaw2d_vtk_write_file (glob, basename,
                                   clawpatch_opt->mx, clawpatch_opt->my, 
                                   clawpatch_opt->meqn,
                                   gparms->vtkspace, gparms->vtkwrite,
                                   fclaw2d_output_vtk_coordinate_cb,
                                   fclaw2d_output_vtk_value_cb);
}