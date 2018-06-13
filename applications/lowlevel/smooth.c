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

#include <fclaw_base.h>
#include <fclaw2d_convenience.h>

typedef struct fclaw_smooth
{
    sc_MPI_Comm mpicomm;
    int mpisize, mpirank;
    fclaw_app_t *a;

    /* parameters */
    int minlevel, maxlevel;
    int smooth_refine, smooth_level, coarsen_delay;

    /* mesh */
    fclaw2d_domain_t *domain;

    /* numerical data */
    double pxy[2];
    double radius, thickn;
}
fclaw_smooth_t;

static int
patch_overlap (fclaw2d_patch_t * patch,
               const double pxy[2], double rmin2, double rmax2)
{
    int i;
    int outside[2];
    double m, hw;
    double ssmin, ssmax;
    double center[2], fdist[2];

    hw = .5 * (patch->xupper - patch->xlower);
    center[0] = patch->xlower + hw;
    center[1] = patch->ylower + hw;

    for (i = 0; i < 2; ++i)
    {
        outside[i] = (fdist[i] = fabs (center[i] - pxy[i])) > hw;
    }

    ssmin = ssmax = 0.;
    for (i = 0; i < 2; ++i)
    {
        if (outside[i])
        {
            m = fdist[i] - hw;
            FCLAW_ASSERT (m >= 0.);
            ssmin += m * m;
        }
        m = fdist[i] + hw;
        ssmax += m * m;
    }
    return ssmin <= rmax2 && rmin2 <= ssmax;
}

static void
run_refine (fclaw_smooth_t * s)
{
    int lev;
    int ib, ip;
    double rmin2, rmax2;
    fclaw2d_block_t *block;
    fclaw2d_patch_t *patch;
    fclaw2d_domain_t *domain;

    /* prepare comparing geometric distance */
    rmin2 = s->radius - s->thickn;
    rmin2 = SC_SQR (rmin2);
    rmax2 = s->radius + s->thickn;
    rmax2 = SC_SQR (rmax2);

    /* loop over multiple initial refinements */
    for (lev = s->minlevel; lev < s->maxlevel; ++lev)
    {
        fclaw_global_productionf ("Initial refinement from level %d\n", lev);
        for (ib = 0; ib < s->domain->num_blocks; ++ib)
        {
            block = s->domain->blocks + ib;
            for (ip = 0; ip < block->num_patches; ++ip)
            {
                /* refine according to overlap with spherical ring */
                patch = block->patches + ip;
                if (patch_overlap (patch, s->pxy, rmin2, rmax2))
                {
                    /* we overlap and prompt refinement of this patch */
                    fclaw2d_patch_mark_refine (s->domain, ib, ip);
                }
            }
        }

        /* run internal mesh refinement */
        domain = fclaw2d_domain_adapt (s->domain);
        FCLAW_ASSERT (domain != s->domain);
        if (domain != NULL)
        {
            fclaw2d_domain_destroy (s->domain);
            s->domain = domain;
            domain = fclaw2d_domain_partition (s->domain, 0);
            FCLAW_ASSERT (domain != s->domain);
            if (domain != NULL)
            {
                fclaw2d_domain_destroy (s->domain);
                s->domain = domain;
                fclaw2d_domain_complete (s->domain);
            }
        }
        else
        {
            /* mesh did not change; initialization is done */
            break;
        }
    }
    fclaw_global_productionf ("Initial refinement to level %d\n", lev);
}

int
main (int argc, char **argv)
{
    fclaw_smooth_t smoo, *s = &smoo;

    s->a = fclaw_app_new (&argc, &argv, NULL);
    s->mpicomm = fclaw_app_get_mpi_size_rank (s->a, &s->mpisize, &s->mpirank);

    /* set parameters */
    s->minlevel = 3;
    s->maxlevel = 4;
    s->smooth_refine = 0;
    s->smooth_level = 0;
    s->coarsen_delay = 0;

    /* init numerical data */
    s->pxy[0] = .3;
    s->pxy[1] = .4;
    s->radius = .2;
    s->thickn = .05;

    /* create a new domain */
    s->domain = fclaw2d_domain_new_unitsquare (s->mpicomm, s->minlevel);
    fclaw2d_domain_set_refinement (s->domain, s->smooth_refine,
                                   s->smooth_level, s->coarsen_delay);

    /* run refinement cycles */
    run_refine (s);

    /* cleanup */
    fclaw2d_domain_destroy (s->domain);
    fclaw_app_destroy (s->a);
    return 0;
}
