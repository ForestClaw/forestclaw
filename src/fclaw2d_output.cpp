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

#include <fclaw2d_output.h>
#include <fclaw2d_clawpatch.hpp>
#include <ClawPatch.H>
/* #include <amr_utils.H> */
#include <fclaw2d_forestclaw.h>
#include <fclaw2d_vtable.h>
#include <fclaw2d_vtk.h>
#include <fclaw2d_map.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

static void
cb_serial_output (fclaw2d_domain_t * domain,
                  fclaw2d_patch_t * this_patch,
                  int this_block_idx, int this_patch_idx,
                  void *user)
{
    fclaw2d_vtable_t vt;
    vt = fclaw2d_get_vtable(domain);

    int iframe = *((int *) user);
    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    int64_t patch_num =
        domain->global_num_patches_before +
        (int64_t) (this_block->num_patches_before + this_patch_idx);

    /* TODO Enable 64bit integers for global counters and indices */

    /* the user can also get this, but maybe we don't want the user
       to have access? */
    int level = this_patch->level;

    vt.patch_write_file(domain, this_patch, this_block_idx,
                        this_patch_idx, iframe, (int) patch_num,
                        level);
}

static
void fclaw2d_output_write_serial(fclaw2d_domain_t* domain,int iframe)
{
    fclaw2d_vtable_t vt;
    vt= fclaw2d_get_vtable(domain);
    /* BEGIN NON-SCALABLE CODE */
    /* Write the file contents in serial.
       Use only for small numbers of processors. */
    fclaw2d_domain_serialization_enter (domain);

    if (domain->mpirank == 0)
    {
        vt.write_header(domain,iframe);
    }

    fclaw2d_domain_iterate_patches (domain, cb_serial_output, (void *) &iframe);
    fclaw2d_domain_serialization_leave (domain);
    /* END OF NON-SCALABLE CODE */
}


static void
fclaw2d_output_vtk_coordinate_cb (fclaw2d_domain_t * domain,
                                  fclaw2d_patch_t * this_patch,
                                  int this_block_idx, int this_patch_idx,
                                  char *a)
{
    // Global parameters
    const amr_options_t *gparms = get_domain_parms (domain);
    fclaw2d_map_context_t *cont;
    const int mx = gparms->mx;
    const int my = gparms->my;

    cont = get_map_context(domain);

    // Patch specific parameters

    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    const double xlower = cp->xlower ();
    const double ylower = cp->ylower ();
    const double dx = cp->dx ();
    const double dy = cp->dy ();

    // Enumerate point coordinates in the patch
    double *d = (double *) a;
    int i, j;
    double xp,yp,zp;
    for (j = 0; j <= my; ++j)
    {
        const double y = ylower + j * dy;
        for (i = 0; i <= mx; ++i)
        {
            const double x = xlower + i * dx;
            if (gparms->manifold)
            {
                FCLAW2D_MAP_C2M(&cont,&this_block_idx,&x,&y,&xp,&yp,&zp);
                *d++ = xp;
                *d++ = yp;
                *d++ = zp;
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
fclaw2d_output_vtk_value_cb (fclaw2d_domain_t * domain,
                             fclaw2d_patch_t * this_patch,
                             int this_block_idx, int this_patch_idx,
                             char *a)
{
    // Global parameters
    const amr_options_t *gparms = get_domain_parms (domain);
    const int mx = gparms->mx;
    const int my = gparms->my;
    const int mbc = gparms->mbc;
    const int meqn = gparms->meqn;
    const int xlane = mx + 2 * mbc;
    const int ylane = my + 2 * mbc;

    // Patch specific parameters
    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    const double *q = cp->q ();

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

void
fclaw2d_output_write_vtk (fclaw2d_domain_t * domain, const char *basename)
{
    const amr_options_t *gparms = get_domain_parms (domain);

    (void) fclaw2d_vtk_write_file (domain, basename,
                                   gparms->mx, gparms->my, gparms->meqn,
                                   gparms->vtkspace, gparms->vtkwrite,
                                   fclaw2d_output_vtk_coordinate_cb,
                                   fclaw2d_output_vtk_value_cb);
}

void
fclaw2d_output_frame (fclaw2d_domain_t * domain, int iframe)
{
    fclaw2d_vtable_t vt;
    double time;
    vt = fclaw2d_get_vtable(domain);

    time = get_domain_time(domain);

    /* Record output time */
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data (domain);
    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_OUTPUT]);

    /* Output VTK file while we're at it */

    const amr_options_t *gparms = get_domain_parms (domain);
    if (gparms->vtkout & 2)
        if (gparms->serialout)
            fclaw_global_essentialf("Output Frame %4d  at time %16.8e (vtk,ascii)\n\n",
                                    iframe,time);
        else
            fclaw_global_essentialf("Output Frame %4d  at time %16.8e (vtk)\n\n",
                                    iframe,time);
    else
        fclaw_global_essentialf("Output Frame %d  at time %16.8e (ascii)\n\n",
                                iframe,time);
    if (gparms->vtkout & 2)
    {
        char basename[BUFSIZ];
        snprintf (basename, BUFSIZ, "%s_frame_%04d", gparms->prefix, iframe);
        fclaw2d_output_write_vtk (domain, basename);
    }

    if (gparms->serialout)
    {
        fclaw2d_output_write_serial(domain, iframe);
    }

    /* Record output time */
    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_OUTPUT]);
}


#ifdef __cplusplus
#if 0
{
#endif
}
#endif
