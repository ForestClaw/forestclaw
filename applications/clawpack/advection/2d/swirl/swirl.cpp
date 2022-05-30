/*
  Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun
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

#include "swirl_user.h"

#include "../all/advection_user.h"

typedef enum
{
    SWIRL_RAY_LINE,
    SWIRL_RAY_CIRCLE,
    SWIRL_RAY_TYPE_LAST
}
swirl_ray_type_t;

typedef struct swirl_ray
{
    swirl_ray_type_t rtype;
    double xy[2];
    union
    {
        struct
        {
            double vec[2];
        } line;
        struct
        {
            double radius;
        } circle;
    } r;
}
swirl_ray_t;

static int
intersect_ray (fclaw2d_domain_t * domain, fclaw2d_patch_t * patch,
               int blockno, int patchno, void *vray, double *integral,
               void *user)
{
    int i, ni;
    double corners[2][2];
    swirl_ray_t *ray = (swirl_ray_t *) vray;

    /*
     * This intersection routine takes the patch coordinate information directly.
     * For mappings of any kind, these would have to be applied here in addition.
     */

    /* assert that ray is a valid swirl_ray_t */
    FCLAW_ASSERT (ray != NULL);
    FCLAW_ASSERT (ray->rtype == SWIRL_RAY_LINE);        /* feel free to add circles */
    FCLAW_ASSERT (integral != NULL && *integral == 0.); /* documented precondition */

    /* just for demonstration purposes: not used in this example */
    FCLAW_ASSERT (user == NULL);

    if (fabs (ray->r.line.vec[0]) <= 1e-12 ||
        fabs (ray->r.line.vec[1]) <= 1e-12)
    {
        /* we cannot guarantee correct results for rays
           that run near parallel to patch boundaries */
        return 0;
    }

    /* store the patch corners in an indexable format */
    corners[0][0] = patch->xlower;
    corners[0][1] = patch->ylower;
    corners[1][0] = patch->xupper;
    corners[1][1] = patch->yupper;

    /* for stability we search in the dimension of the strongest component */
    i = (fabs (ray->r.line.vec[0]) <= fabs (ray->r.line.vec[1])) ? 1 : 0;
    ni = i ^ 1; /* not i */

    if (patchno >= 0)
    {
        /* We are at a leaf and the patch is a valid patch of the domain.
         * Based on the patch, the domain, the blockno and the information stored
         * in the swirl_ray_t we defined, we now have to set *integral to be the
         * contribution of this ray-patch combination to the ray integral.
         * We should return 1 (even though a leaf return value is ignored). */
        int j, nj;
        double t, shift;
        double hits[2][2];

        /* compute the coordinates of intersections with most orthogonal edges */
        t = (corners[0][i] - ray->xy[i]) / ray->r.line.vec[i];
        shift = ray->xy[ni] + t * ray->r.line.vec[ni];

        /* shift coordinate system to first hit */
        hits[0][0] = 0.;
        hits[0][1] = 0.;
#if 0
        corners[1][i] -= corners[0][i];
        corners[0][i] = 0.;
#endif
        corners[1][ni] -= shift;
        corners[0][ni] -= shift;

        /* compute second hit in shifted coordinate system */
        t = (corners[1][i] - corners[0][i]) / ray->r.line.vec[i];
        hits[1][i] = (corners[1][i] - corners[0][i]);
        hits[1][ni] = t * ray->r.line.vec[ni];

        /* compute the actual hit coordinates */
        t = 0.;
        for (j = 0; j < 2; j++)
        {
            nj = j ^ 1;
            /* check if we hit the faces parallel to the search dimension */
            if (hits[j][ni] < corners[0][ni])
            {
                t = (corners[0][ni] - hits[j][ni]) /
                    (hits[nj][ni] - hits[j][ni]);
                hits[j][ni] = corners[0][ni];
            }
            else if (hits[j][ni] > corners[1][ni])
            {
                t = (corners[1][ni] - hits[j][ni]) /
                    (hits[nj][ni] - hits[j][ni]);
                hits[j][ni] = corners[1][ni];
            }
            else
            {
                /* if we get here we hit a face orthogonal to the search dimension */
                continue;
            }

            if (t < 0. || t > 1.)
            {
                return 0;       /* if t is not in [0, 1], both hits lie outside the patch */
            }
            hits[j][i] = hits[j][i] * (1 - t) + hits[nj][i] * t;
        }

        /* compute the distance of the two hits and store it in integral */
        hits[0][0] = hits[0][0] - hits[1][0];
        hits[0][1] = hits[0][1] - hits[1][1];
        *integral = sqrt (hits[0][0] * hits[0][0] + hits[0][1] * hits[0][1]);
        return 1;
    }
    else
    {
        /* We are not at a leaf and the patch is an artificial patch containing all
         * standard patch information except for the pointer to the next patch and
         * user-data of any kind.
         * Only FCLAW2D_PATCH_CHILDID and FCLAW2D_PATCH_ON_BLOCK_FACE_* flags are set.
         * Based on this, we run a test to test whether the ray and the patch intersect.
         *
         * We return 0 if we are certain that the ray does not intersect any
         * descendant of this patch.
         * We return 1 if the test concludes that the ray may intersect the patch.
         * This test may be overinclusive (a false positive) to optimize for speed.
         *
         * The purpose of this test is to remove irrelevant ancestor
         * patch-ray-combinations early on to avoid unnecessary computation.
         * We do not need to assign to the integral value for ancestor patches. */
        int j, is_left;
        double t, rayatt;

        is_left = 0;
        for (j = 0; j < 2; j++)
        {
            /* compute the coordinates of the ray when it intersects the patch
             * interval of the search dimension */
            t = (corners[j][i] - ray->xy[i]) / ray->r.line.vec[i];
            rayatt = ray->xy[ni] + t * ray->r.line.vec[ni];

            if (rayatt < corners[0][ni])
            {
                is_left++;      /* we pass the patch to the left (or bottom) */
            }
            else if (rayatt <= corners[1][ni])
            {
                return 1;       /* hits the edge between corners[0][ni] and corners[1][ni] */
            }
        }
        /* If we reach this point, there was no direct hit in the search dimension.
         * If there is one point on each side of the patch, there will be a hit in
         * the other dimension, else there will be none. */
        return (is_left == 1);
    }
}

static void
print_integrals (fclaw2d_domain_t * domain, sc_array_t * rays,
                 sc_array_t * integrals)
{
    if (domain->mpirank == 0)
    {
        int i;
        double *integral;
        swirl_ray_t *ray;

        /* replace with proper logging as desired */
        fprintf (stderr, "Results of the domain integration:\n");
        for (i = 0; i < (int) rays->elem_count; i++)
        {
            ray = (swirl_ray_t *) sc_array_index_int (rays, i);
            integral = (double *) sc_array_index_int (integrals, i);
            fprintf (stderr,
                     "Ray %d: [%2.5f,%2.5f] + t * [%2.5f,%2.5f] integral %g\n",
                     i, ray->xy[0], ray->xy[1], ray->r.line.vec[0],
                     ray->r.line.vec[1], *integral);
        }
    }
}

static
fclaw2d_domain_t* create_domain(sc_MPI_Comm mpicomm, fclaw_options_t* gparms)
{
    /* Mapped, multi-block domain */
    p4est_connectivity_t     *conn = NULL;
    fclaw2d_domain_t         *domain;
    fclaw2d_map_context_t    *cont = NULL;

    /* Map unit square to disk using mapc2m_disk.f */
    gparms->manifold = 0;
    conn = p4est_connectivity_new_unitsquare();
    cont = fclaw2d_map_new_nomap();

    domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);
    fclaw2d_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw2d_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);
    return domain;
}

static
void run_program(fclaw2d_global_t* glob, sc_array_t *rays, sc_array_t *integrals)
{
    const user_options_t           *user_opt;

    /* ---------------------------------------------------------------
       Set domain data.
       --------------------------------------------------------------- */
    fclaw2d_domain_data_new(glob->domain);

    user_opt = swirl_get_options(glob);

    /* Initialize virtual table for ForestClaw */
    fclaw2d_vtables_initialize(glob);

    /* Initialize virtual tables for solvers */
    if (user_opt->claw_version == 4)
    {
        fc2d_clawpack46_solver_initialize();
    }
    else if (user_opt->claw_version == 5)
    {
        fc2d_clawpack5_solver_initialize();
    }

    swirl_link_solvers(glob);

    /* ---------------------------------------------------------------
       Run
       --------------------------------------------------------------- */
    fclaw2d_initialize(glob);
    fclaw2d_run(glob);

    /* For convenience, we integrate after the run of the solver has finished.
     * In practice, it may be of interest to integrate several times during the
     * run for different time steps.  To implement this, move the following
     * call into a repeated diagnostic steps and output the integral values.
     */
    fclaw2d_domain_integrate_rays (glob->domain,
                                   intersect_ray, rays, integrals, NULL);
    print_integrals (glob->domain, rays, integrals);

    fclaw2d_finalize(glob);
}

static sc_array_t *
swirl_rays_new (int nlines)
{
    int i;
    swirl_ray_t *ray;
    sc_array_t *a = sc_array_new (sizeof (swirl_ray_t));

    /* add a couple straight rays */
    FCLAW_ASSERT (nlines >= 0);
    for (i = 0; i < nlines; ++i)
    {
        ray = (swirl_ray_t *) sc_array_push (a);
        ray->rtype = SWIRL_RAY_LINE;
        ray->xy[0] = 0.5;
        ray->xy[1] = 0.5;
        /* we add 0.1, since the intersection callback does not guarantee exact
         * results for axis-parallel rays, which may lie on a patch boundary */
        ray->r.line.vec[0] = cos ((i + 0.1) * 2 * M_PI / nlines);
        ray->r.line.vec[1] = sin ((i + 0.1) * M_PI / nlines);
    }

    /* add no circles yet */
    return a;
}

static sc_array_t *
swirl_integrals_new (int nlines)
{
    FCLAW_ASSERT (nlines >= 0);
    return sc_array_new_count (sizeof (double), nlines);
}

int
main (int argc, char **argv)
{
    /* number of rays that are straight lines.  No circles yet */
    const int swirl_nlines = 8;

    int first_arg;
    sc_array_t *rays;
    sc_array_t *integrals;
    fclaw_app_t *app;
    fclaw_exit_type_t vexit;

    /* Options */
    sc_options_t                *options;
    user_options_t              *user_opt;
    fclaw_options_t             *fclaw_opt;
    fclaw2d_clawpatch_options_t *clawpatch_opt;
    fc2d_clawpack46_options_t   *claw46_opt;
    fc2d_clawpack5_options_t    *claw5_opt;

    fclaw2d_global_t            *glob;
    fclaw2d_domain_t            *domain;
    sc_MPI_Comm mpicomm;

    int retval;

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, NULL);

    /* Create new options packages */
    fclaw_opt =                   fclaw_options_register(app,"fclaw_options.ini");
    clawpatch_opt =   fclaw2d_clawpatch_options_register(app,"fclaw_options.ini");
    claw46_opt =        fc2d_clawpack46_options_register(app,"fclaw_options.ini");
    claw5_opt =          fc2d_clawpack5_options_register(app,"fclaw_options.ini");
    user_opt =                    swirl_options_register(app,"fclaw_options.ini");

    /* Read configuration file(s) and command line, and process options */
    options = fclaw_app_get_options (app);
    retval = fclaw_options_read_from_file(options);
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    /* Setup some rays to integrate along/around */
    rays = swirl_rays_new (swirl_nlines);
    integrals = swirl_integrals_new (swirl_nlines);

    /* Run the program */
    if (!retval & !vexit)
    {
        /* Options have been checked and are valid */

        mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);
        domain = create_domain(mpicomm, fclaw_opt);

        /* Create global structure which stores the domain, timers, etc */
        glob = fclaw2d_global_new();
        fclaw2d_global_store_domain(glob, domain);

        /* Store option packages in glob */
        fclaw2d_options_store           (glob, fclaw_opt);
        fclaw2d_clawpatch_options_store (glob, clawpatch_opt);
        fc2d_clawpack46_options_store   (glob, claw46_opt);
        fc2d_clawpack5_options_store    (glob, claw5_opt);
        swirl_options_store             (glob, user_opt);

        run_program(glob, rays, integrals);
        fclaw2d_global_destroy(glob);
    }

    sc_array_destroy (rays);
    sc_array_destroy (integrals);
    fclaw_app_destroy (app);

    return 0;
}
