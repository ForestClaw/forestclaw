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
#include <fclaw_convenience.h>

#define TBUFSIZ (2 * BUFSIZ)

typedef struct fclaw_smooth
{
    sc_MPI_Comm mpicomm;
    int mpisize, mpirank;
    fclaw_app_t *a;

    /* parameters */
    int minlevel, maxlevel;
    int smooth_refine, smooth_level, coarsen_delay;
    int write_vtk;
    double dt, finalt;
    double radius, thickn;

    /* derived data */
    double rmin2, rmax2;
    char prefix[BUFSIZ];

    /* mesh */
    fclaw_domain_t *domain;

    /* numerical data */
    int k;
    double time;
    double pxy[2];
    double vel[2];
}
fclaw_smooth_t;

#ifdef FCLAW_ENABLE_DEBUG

static int
domain_match (fclaw_domain_t * d1, fclaw_domain_t * d2)
{
    return 1;
}

#endif

static void
init_values (fclaw_smooth_t * s)
{
    double ss;

    /* prepare comparing geometric distance */
    ss = s->radius - s->thickn;
    s->rmin2 = SC_SQR (ss);
    ss = s->radius + s->thickn;
    s->rmax2 = SC_SQR (ss);
    snprintf (s->prefix, BUFSIZ, "sm_l%dL%d_s%dl%dc%d",
              s->minlevel, s->maxlevel,
              s->smooth_refine, s->smooth_level, s->coarsen_delay);
    fclaw_global_productionf ("Output prefix %s\n", s->prefix);

    /* init numerical data */
    s->pxy[0] = .4;
    s->pxy[1] = .3;
    s->vel[0] = .6;
    s->vel[1] = .8;
}

static int
patch_overlap (fclaw_patch_t * patch,
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
mark_patches (fclaw_smooth_t * s)
{
    int ib, ip;
    fclaw_block_t *block;
    fclaw_patch_t *patch;

    for (ib = 0; ib < s->domain->num_blocks; ++ib)
    {
        block = s->domain->blocks + ib;
        for (ip = 0; ip < block->num_patches; ++ip)
        {
            /* refine according to overlap with spherical ring */
            patch = block->patches + ip;
            if (patch_overlap (patch, s->pxy, s->rmin2, s->rmax2))
            {
                /* we overlap and prompt refinement of this patch */
                if (patch->level < s->maxlevel)
                {
                    fclaw_patch_mark_refine (s->domain, ib, ip);
                }
            }
            else
            {
                /* we coarsen if we do not overlap */
                if (patch->level > s->minlevel)
                {
                    fclaw_patch_mark_coarsen (s->domain, ib, ip);
                }
            }

        }
    }
}

static void
init_refine (fclaw_smooth_t * s)
{
    int lev;
    char basename[TBUFSIZ];
    fclaw_domain_t *domain;

    /* loop over multiple initial refinements */
    fclaw_global_infof ("Initial position %g %g\n", s->pxy[0], s->pxy[1]);
    for (lev = s->minlevel; lev <= s->maxlevel; ++lev)
    {
        fclaw_global_productionf ("Initial pseudo-level %d\n", lev);
        if (s->write_vtk)
        {
            snprintf (basename, TBUFSIZ, "%s_L%02d", s->prefix, lev);
            fclaw_domain_write_vtk (s->domain, basename);
        }

        /* run refinement indicator for all patches */
        mark_patches (s);

        /* run internal mesh refinement */
        domain = fclaw_domain_adapt (s->domain);
        FCLAW_ASSERT (domain != s->domain);
        if (domain != NULL)
        {
            FCLAW_ASSERT (domain_match (domain, s->domain));
            fclaw_domain_destroy (s->domain);
            s->domain = domain;
            domain = fclaw_domain_partition (s->domain, 0);
            FCLAW_ASSERT (domain != s->domain);
            if (domain != NULL)
            {
                FCLAW_ASSERT (domain_match (domain, s->domain));
                fclaw_domain_destroy (s->domain);
                s->domain = domain;
                fclaw_domain_complete (s->domain);
            }
        }
        else
        {
            /* mesh did not change; initialization is done */
            break;
        }
    }
    fclaw_global_productionf ("Initial pseudo-level %d\n", lev);
    if (s->write_vtk)
    {
        snprintf (basename, TBUFSIZ, "%s_L%02d", s->prefix, lev);
        fclaw_domain_write_vtk (s->domain, basename);
    }
}

static void
run_refine (fclaw_smooth_t * s)
{
    double nextt, deltat;
    char basename[TBUFSIZ];
    fclaw_domain_t *domain;

    /* initialize time stepping */
    s->k = 0;
    s->time = 0.;
    fclaw_global_productionf ("Run time %.3g step %d\n", s->time, s->k);
    if (s->write_vtk)
    {
        snprintf (basename, TBUFSIZ, "%s_K%05d", s->prefix, s->k);
        fclaw_domain_write_vtk (s->domain, basename);
    }

    /* run time loop */
    while (s->time < s->finalt)
    {
        /* adjust time step near final time */
        nextt = s->time + (deltat = s->dt);
        if (nextt > s->finalt - 1e-3 * deltat)
        {
            nextt = s->finalt;
            deltat = nextt - s->time;
        }

        /* advance position of sphere */
        s->pxy[0] += deltat * s->vel[0];
        s->pxy[1] += deltat * s->vel[1];
        fclaw_global_infof ("New position %g %g\n", s->pxy[0], s->pxy[1]);

        /* run refinement indicator for all patches */
        mark_patches (s);

        /* run internal mesh refinement */
        domain = fclaw_domain_adapt (s->domain);
        FCLAW_ASSERT (domain != s->domain);
        if (domain != NULL)
        {
            FCLAW_ASSERT (domain_match (domain, s->domain));
            fclaw_domain_destroy (s->domain);
            s->domain = domain;
            domain = fclaw_domain_partition (s->domain, 0);
            FCLAW_ASSERT (domain != s->domain);
            if (domain != NULL)
            {
                FCLAW_ASSERT (domain_match (domain, s->domain));
                fclaw_domain_destroy (s->domain);
                s->domain = domain;
                fclaw_domain_complete (s->domain);
            }
        }

        /* advance time */
        ++s->k;
        s->time = nextt;
        fclaw_global_productionf ("Run time %.3g step %d\n", s->time, s->k);
        if (s->write_vtk)
        {
            snprintf (basename, TBUFSIZ, "%s_K%05d", s->prefix, s->k);
            fclaw_domain_write_vtk (s->domain, basename);
        }
    }
}

int
main (int argc, char **argv)
{
    fclaw_smooth_t smoo, *s = &smoo;

    s->a = fclaw_app_new (&argc, &argv, NULL);
    s->mpicomm = fclaw_app_get_mpi_size_rank (s->a, &s->mpisize, &s->mpirank);

    /* set refinement parameters */
    s->minlevel = 2;
    s->maxlevel = 7;
    s->smooth_refine = 1;
    s->smooth_level = 6;
    s->coarsen_delay = 0;
    s->write_vtk = 1;

    /* set numerical parameters */
    s->dt = .02;
    s->finalt = .3;
    s->radius = .2;
    s->thickn = .0 * s->radius;

    /* init numerical data */
    init_values (s);

    /* create a new domain */
    s->domain = fclaw_domain_new_unitsquare (s->mpicomm, s->minlevel);
    fclaw_domain_set_refinement (s->domain, s->smooth_refine,
                                 s->smooth_level, s->coarsen_delay);

    /* run initial refinement loop */
    init_refine (s);

    /* run refinement/coarsening while sphere is moving */
    run_refine (s);

    /* cleanup */
    fclaw_domain_destroy (s->domain);
    fclaw_app_destroy (s->a);
    return 0;
}
