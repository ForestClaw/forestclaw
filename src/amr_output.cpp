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

#include "amr_forestclaw.H"
#include "amr_utils.H"
#include "clawpack_fort.H"
#include "amr_output.H"
#include <fclaw2d_vtk.h>

static
void cb_amrout(fclaw2d_domain_t *domain,
               fclaw2d_patch_t *this_patch,
               int this_block_idx,
               int this_patch_idx,
               void *user)
{
    int iframe = *((int *) user);
    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    int64_t patch_num =
      domain->global_num_patches_before +
      (int64_t) (this_block->num_patches_before + this_patch_idx);

    /* TODO Enable 64bit integers for global counters and indices */

    /* the user can also get this, but maybe we don't want the user
       to have access? */
    int level = this_patch->level;

    fclaw2d_output_functions_t* of = get_output_functions(domain);
    (of->f_patch_write_output)(domain,this_patch,this_block_idx,
                               this_patch_idx,iframe,(int) patch_num,level);
}


void amrout(fclaw2d_domain_t *domain, int iframe)
{
    // Get total number of patches
    int ngrids = domain->global_num_patches;
    fclaw2d_output_functions_t* of = get_output_functions(domain);

    /* BEGIN NON-SCALABLE CODE */
    /* Write the file contents in serial.
       Use only for small numbers of processors. */
    fclaw2d_domain_serialization_enter (domain);

    if (domain->mpirank == 0) {
        /* the header needs to be written by the first processor */
        (of->f_patch_write_header)(domain,iframe,ngrids);
    }
    fclaw2d_domain_iterate_patches(domain, cb_amrout, (void *) &iframe);

    fclaw2d_domain_serialization_leave (domain);
    /* END OF NON-SCALABLE CODE */
}

static void
amr_output_vtk_coordinate_cb (fclaw2d_domain_t * domain,
                              fclaw2d_patch_t * this_patch,
                              int this_block_idx, int this_patch_idx,
                              char *a)
{
    // In case this is needed by the setaux routine
    set_block_(&this_block_idx);

    // Global parameters
    const amr_options_t *gparms = get_domain_parms(domain);
    const int mx = gparms->mx;
    const int my = gparms->my;

    // Patch specific parameters
    ClawPatch *cp = get_clawpatch(this_patch);
    const double xlower = cp->xlower();
    const double ylower = cp->ylower();
    const double dx = cp->dx();
    const double dy = cp->dy();

    // Enumerate point coordinates in the patch
    double * d = (double *) a;
    int i, j;
    for (j = 0; j <= my; ++j)
    {
       const double y = ylower + j * dy;
       for (i = 0; i <= mx; ++i)
       {
           *d++ = xlower + i * dx;
           *d++ = y;
           *d++ = 0.;
       }
    }
}

static void
amr_output_vtk_value_cb (fclaw2d_domain_t * domain,
                         fclaw2d_patch_t * this_patch,
                         int this_block_idx, int this_patch_idx, char *a)
{
    // In case this is needed by the setaux routine
    set_block_(&this_block_idx);

    // Global parameters
    const amr_options_t *gparms = get_domain_parms(domain);
    const int mx = gparms->mx;
    const int my = gparms->my;
    const int mbc = gparms->mbc;
    const int meqn = gparms->meqn;
    const int xlane = mx + 2 * mbc;
    const int ylane = my + 2 * mbc;

    // Patch specific parameters
    ClawPatch *cp = get_clawpatch(this_patch);
    const double *q = cp->q();

    // Enumerate equation data in the patch
    float * f = (float *) a;
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
amr_output_write_vtk (fclaw2d_domain_t *domain, const char *basename)
{
    const amr_options_t *gparms = get_domain_parms(domain);

    fclaw2d_vtk_write_file (domain, basename,
                            gparms->mx, gparms->my, gparms->meqn,
                            gparms->vtkspace, gparms->vtkwrite,
                            amr_output_vtk_coordinate_cb,
                            amr_output_vtk_value_cb);
}
