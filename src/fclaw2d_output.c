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


static void
cb_serial_output (fclaw2d_domain_t * domain,
                  fclaw2d_patch_t * this_patch,
                  int this_block_idx, int this_patch_idx,
                  void *user)
{
    fclaw2d_global_iterate_t* s = (fclaw2d_global_iterate_t*) user;

    int iframe = *((int *) s->user);
    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    int64_t patch_num =
        domain->global_num_patches_before +
        (int64_t) (this_block->num_patches_before + this_patch_idx);

    /* TODO Enable 64bit integers for global counters and indices */

    /* the user can also get this, but maybe we don't want the user
       to have access? */
    int level = this_patch->level;

    fclaw2d_patch_write_file(s->glob, this_patch, this_block_idx,
                             this_patch_idx, iframe, (int) patch_num,
                             level);
}

static
void fclaw2d_output_write_serial(fclaw2d_global_t* glob,int iframe)
{
    fclaw2d_domain_t *domain = glob->domain;

    /* BEGIN NON-SCALABLE CODE */
    /* Write the file contents in serial.
       Use only for small numbers of processors. */
    fclaw2d_domain_serialization_enter (domain);

    if (domain->mpirank == 0)
    {
        fclaw2d_patch_write_header(glob,iframe);
    }

    /* Write out each patch to fort.qXXXX */
    fclaw2d_global_iterate_patches (glob, cb_serial_output, (void *) &iframe);
    fclaw2d_domain_serialization_leave (domain);
    /* END OF NON-SCALABLE CODE */
}

void
fclaw2d_output_frame (fclaw2d_global_t * glob, int iframe)
{
    double time;

    time = glob->curr_time;

    /* Record output time */
    fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_OUTPUT]);

    /* Output VTK file while we're at it */

    const amr_options_t *gparms = fclaw2d_forestclaw_get_options(glob);
    if (gparms->vtkout & 2)
    {
        if (gparms->serialout)
            fclaw_global_essentialf("Output Frame %4d  at time %16.8e (vtk,ascii)\n\n",
                                    iframe,time);
        else
            fclaw_global_essentialf("Output Frame %4d  at time %16.8e (vtk)\n\n",
                                    iframe,time);
    }
    else if (gparms->serialout)
    {
        fclaw_global_essentialf("Output Frame %d  at time %16.8e (ascii)\n\n",
                                iframe,time);
    }
    else
    {
        fclaw_global_essentialf("Output Frame %d  at time %16.8e (none)\n\n",
                                iframe,time);
    }

    if (gparms->vtkout & 2)
    {
        char basename[BUFSIZ];
        snprintf (basename, BUFSIZ, "%s_frame_%04d", gparms->prefix, iframe);
        fclaw2d_output_write_vtk (glob, basename);
    }

    if (gparms->serialout)
    {
        fclaw2d_output_write_serial(glob, iframe);
    }

    if (gparms->tikzout)
    {
        fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_EXTRA3]);
        fclaw2d_output_write_tikz(glob,iframe);
        fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_EXTRA3]);
    }

    /* Record output time */
    fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_OUTPUT]);
}

