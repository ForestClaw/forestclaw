/*
  Copyright (c) 2012-2023 Carsten Burstedde, Donna Calhoun, Scott Aiton
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
#include <fclaw2d_rays.h>

#include "../all/advection_user.h"

/* A #define of 0 enables a set of rays that is the same as in the swirl
 * example parallel to this application directory.
 * If #define is 1, use an ad-hoc setup local to this file.
 */
#if 1
#define STAR_OF_RAYS
#endif

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
            /* set to 0 if ray is approximately parallel to the x axis,
               set to 1 if ray is approximately parallel to the y axis,
               and to 2 otherwise.  Will only work for Cartesian setup.
             */
            int parallel;
            /* set to 0 if the ray is more parallel to x than y axis. */
            int dominant;
            double vec[2];
        } line;
        struct
        {
            double radius;
        } circle;
    } r;
}
swirl_ray_t;

#ifdef STAR_OF_RAYS
const int swirl_nlines = 8;
#else
const int swirl_nlines = 3;
#endif

/* Virtual function for setting rays */
static void
swirl_allocate_and_define_rays (fclaw_global_t * glob,
                                fclaw2d_ray_t ** rays, int *num_rays)
{
    int i;

    *num_rays = swirl_nlines;

    /* We let the user allocate an array of rays, although what is inside the
       generic ray type is left opaque. This is destroy in matching FREE,
       below. */

    *rays = fclaw2d_ray_allocate_rays(*num_rays);
    fclaw2d_ray_t *ray_vec = *rays;
    for (i = 0; i < swirl_nlines; ++i)
    {
#ifndef STAR_OF_RAYS
        double dth;
#endif
        swirl_ray_t *sr = (swirl_ray_t *) FCLAW_ALLOC (swirl_ray_t, 1);
        sr->rtype = SWIRL_RAY_LINE;

        /* we define all vecs to have norm 1 to make integration easier */
#ifdef STAR_OF_RAYS
        sr->xy[0] = 0.5;
        sr->xy[1] = 0.5;
        sr->r.line.vec[0] = cos ((i) * 2 * M_PI / swirl_nlines);
        sr->r.line.vec[1] = sin ((i) * 2 * M_PI / swirl_nlines);
#else
        /* End points are on a semi-circle in x>0,y>0 quad */
        FCLAW_ASSERT(swirl_nlines >= 2);
        sr->xy[0] = 0; /* -0.1; */
        sr->xy[1] = 0; /* -0.1; */
        dth = M_PI/(2*swirl_nlines);
        sr->r.line.vec[0] = cos ((i + 0.5) * dth);
        sr->r.line.vec[1] = sin ((i + 0.5) * dth);
#endif

        /* Determine whether a ray is axis parallel */
        if (fabs (sr->r.line.vec[1]) <= 1e-12) {
          sr->r.line.parallel = 0;
        }
        else if (fabs (sr->r.line.vec[0]) <= 1e-12) {
          sr->r.line.parallel = 1;
        }
        else {
          sr->r.line.parallel = 2;
        }
        sr->r.line.dominant =
          fabs (sr->r.line.vec[0]) <= fabs (sr->r.line.vec[1]) ? 1 : 0;
        /* check sanity of computation of parallel and dominant */
        FCLAW_ASSERT(sr->r.line.parallel == 2 ||
                     sr->r.line.parallel == sr->r.line.dominant);

        /* Assign ray to diagnostics item */
        fclaw2d_ray_t *ray = &ray_vec[i];
        fclaw2d_ray_set_ray (ray, i + 1, sr);
    }
}

static
void swirl_deallocate_rays(fclaw_global_t *glob,
                           fclaw2d_ray_t** rays,
                           int* num_rays)
{
    int i;
    fclaw2d_ray_t *ray_vec = *rays;

    for(i = 0; i < *num_rays; i++)
    {
        /* Retrieve rays set above and deallocate them */
        int id;
        fclaw2d_ray_t *ray = &ray_vec[i];
        swirl_ray_t *rs = (swirl_ray_t*) fclaw2d_ray_get_ray(ray, &id);
        FCLAW_ASSERT (rs != NULL);
        FCLAW_FREE (rs);
        rs = NULL;
    }
    /* Match FCLAW_ALLOC, above */
    *num_rays = fclaw2d_ray_deallocate_rays(rays);
}

/** This function checks if a linear ray intersects a patch.
 * In case the ray does not intersect the patch it returns 0.
 * In case the ray does intersect the patch it returns 1, dt contains the
 * length of the intersection (as long as the vec has norm 1) and rayni
 * contains the not-i-coordinates of the ray for both its intersections with
 * the patch-boundary in the search dimension i.
 */
static int
intersect_patch (fclaw_patch_t *patch, swirl_ray_t *swirl_ray,
                 int i, int *untrustworthy, double *dt, double rayni[2])
{
    int ni, j, isleft, iscenter;
    double corners[2][2], h = 0.;

    FCLAW_ASSERT (untrustworthy != NULL);
    FCLAW_ASSERT (dt != NULL);

    /* store the patch corners in an indexable format */
    corners[0][0] = patch->d2->xlower;
    corners[0][1] = patch->d2->ylower;
    corners[1][0] = patch->d2->xupper;
    corners[1][1] = patch->d2->yupper;

    /* compute the coordinates of intersections with most orthogonal edges */
    ni = i ^ 1;
    rayni[0] = swirl_ray->xy[ni] + swirl_ray->r.line.vec[ni] *
                                   ((corners[0][i] - swirl_ray->xy[i]) / swirl_ray->r.line.vec[i]);
    rayni[0] -= corners[0][ni]; /* shift coordinate system to lower left corner */

    /* compute second hit in shifted coordinate system */
    *dt = (corners[1][i] - corners[0][i]) / swirl_ray->r.line.vec[i];
    rayni[1] = rayni[0] + *dt * swirl_ray->r.line.vec[ni];

    if ((j = swirl_ray->r.line.parallel) != 2)
    {
        /* this ray is parallel to one axis:
           we will not integrate ray-patch-intersections that are closer than
           h to the patch faces in dimension ni, because we cannot guarantee
           reliable results in this case. */
        FCLAW_ASSERT (i == j);
        h = 2e-12 * (corners[1][j ^ 1] - corners[0][j ^ 1]);
    }

    isleft = iscenter = 0;
    for (j = 0; j < 2; j++)
    {
        if (rayni[j] < 0.)
        {
            ++isleft;       /* we pass the patch to the left (or bottom) */
        }
        else if (rayni[j] <= (corners[1][ni] - corners[0][ni]))
        {
            ++iscenter;     /* hits the edge between corners[0][ni] and corners[1][ni] */
            if (swirl_ray->r.line.parallel != 2)
            {
                if (rayni[j] < h ||
                    rayni[j] > (corners[1][ni] - corners[0][ni]) - h) {
                    /* the ray hits the patch too close to a face */
                    *untrustworthy = 1;
                    return 0;
                }
            }
        }
    }

    /* verify that we caught all axis parallel rays that intersect the patch
     * boundary */
    FCLAW_ASSERT(swirl_ray->r.line.parallel == 2 || iscenter != 1);

    /* If the ray misses the patch there will be two hits on the same side
     * of the patch. */
    return isleft == 1 || iscenter;
}

static int
swirl_intersect_ray (fclaw_domain_t *domain, fclaw_patch_t *patch,
                     int blockno, int patchno, void *ray, double *integral,
                     void *user)
{
    int i, id;
    double dt, rayni[2];

    /* assert that ray is a valid swirl_ray_t */
    fclaw2d_ray_t *fclaw_ray = (fclaw2d_ray_t *) ray;

    swirl_ray_t *swirl_ray = (swirl_ray_t*) fclaw2d_ray_get_ray(fclaw_ray,&id);
    FCLAW_ASSERT(swirl_ray != NULL);
    FCLAW_ASSERT(swirl_ray->rtype == SWIRL_RAY_LINE); /* Circles not there yet. */
    FCLAW_ASSERT (integral != NULL && *integral == 0.); /* documented precondition */

    /*
     * This intersection routine takes the patch coordinate information directly.
     * For mappings of any kind, these would have to be applied here in addition.
     * We cannot guarantee correct results for rays that run near parallel to
     * any coordinate axis, because they might cause undefined behaviour.
     * We propose to rotate axis-parallel rays such that all entries of the
     * vector are at least 1e-12.
     */

    /* for stability we search in the dimension of the strongest component */
    i = swirl_ray->r.line.dominant;

    if (patchno >= 0)
    {
        /* We are at a leaf and the patch is a valid patch of the domain.
         * Based on the patch, the domain, the blockno and the information stored
         * in the swirl_ray_t we defined, we now have to set *integral to be the
         * contribution of this ray-patch combination to the ray integral.
         * We should return 1 (even though a leaf return value is ignored). */
        if (!intersect_patch (patch, swirl_ray, i, &fclaw_ray->untrustworthy,
                              &dt, rayni)) {
            /* We do not have an intersection, the return value is ignored. */
            return 1;
        }

        int ni = i ^ 1;
        int mx, my, mbc, meqn, sol_rows, j, k, klower, kupper, k_current;
        double xlower, ylower, dx, dy, *sol, tstep,
               nilower, niupper, raynilower, rayniupper;
        fclaw_global_t *glob = (fclaw_global_t *) user;
        FCLAW_ASSERT(glob != NULL);

        /* Obtain cell indices of the hits. */
        fclaw2d_clawpatch_grid_data(glob, patch, &mx, &my, &mbc,
                                    &xlower, &ylower, &dx, &dy);
        fclaw_clawpatch_soln_data(glob, patch, &sol, &meqn);
        sol_rows = mx + 2 * mbc;

        /* Make mx, my, dx, dy indexable. */
        int m[2] = {mx, my};
        double d[2] = {dx, dy};

        /* We will apply line-tracing to follow the ray through the cell-grid.
         * For this we divide dt in m[i] equal steps tstep to end up with the
         * values of the ray on the grid-boundaries. */
        tstep = dt/m[i];
        rayni[1] = rayni[0];
        for(j = 0; j < m[i]; j++) {
            /* Compute the range of the ni-coordinate of the next tstep and the
             * corresponding indices in the mx x my cell-grid. We work on a
             * coordinate system with (0,0) in the lower left patch corner. */
            rayni[0] = rayni[1];
            rayni[1] += tstep * swirl_ray->r.line.vec[ni];
            raynilower = SC_MIN(rayni[0], rayni[1]);
            rayniupper = SC_MAX(rayni[0], rayni[1]);
            klower = floor(raynilower / d[ni]);
            kupper = floor(rayniupper / d[ni]);

            if(kupper < 0 || klower >= m[ni]) {
                /* The current cell does not belong to the patch. */
                continue;
            }

            if(klower == kupper) {
                /* The solution is constant on every cell, so we are only
                 * interested in the length of the intersection.
                 * Since we demanded the ray-vector to have length 1, we can
                 * obtain it directly from the difference fabs(tstep) in t.
                 * In other settings one might compute the actual coordinates
                 * of the intersection of the ray with the cell boundaries
                 * and apply more advanced numerical integration schemes here. */
                *integral += sol[(mbc + klower) * sol_rows + (mbc + j)] * fabs(tstep);
            } else {
                for(k = klower; k < kupper; k++) {
                    /* compute the ni-range of the current cell. All values
                     * outside the patch are ignored */
                    nilower = (k == klower) ? SC_MAX(raynilower, 0.) : k * d[ni];
                    niupper = (k + 1 == kupper) ? SC_MIN(rayniupper, m[ni] * d[ni])
                            : (k + 1) * d[ni];
                    dt = tstep * (niupper - nilower) / (rayni[1] - rayni[0]);

                    /* We do not want to use k as the cell-index, because it may
                     * be afflicted by rounding errors. Instead, we use the
                     * weighted mean of niupper and nilower to compute
                     * the current cell-index. */
                    k_current = floor((niupper + nilower) / (2. * d[ni]));
                    *integral += sol[(mbc + k_current) * sol_rows + (mbc + j)] * fabs(dt);
                }
            }
        }
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
        return intersect_patch (patch, swirl_ray, i, &fclaw_ray->untrustworthy,
                                &dt, rayni);
    }
}

void swirl_initialize_rays(fclaw_global_t* glob)
{
    /* Set up rays */
    fclaw2d_ray_vtable_t* rays_vt = fclaw2d_ray_vt(glob);

    rays_vt->allocate_and_define = swirl_allocate_and_define_rays;
    rays_vt->deallocate = swirl_deallocate_rays;

    rays_vt->integrate = swirl_intersect_ray;
}

static
fclaw_domain_t* create_domain(sc_MPI_Comm mpicomm, fclaw_options_t* gparms)
{
    /* Mapped, multi-block domain */
    p4est_connectivity_t     *conn = NULL;
    fclaw_domain_t         *domain;
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
void run_program(fclaw_global_t* glob)
{
    const user_options_t           *user_opt;

    /* ---------------------------------------------------------------
       Set domain data.
       --------------------------------------------------------------- */
    fclaw_domain_data_new(glob->domain);

    user_opt = swirl_get_options(glob);

    /* Initialize virtual table for ForestClaw */
    fclaw2d_vtables_initialize(glob);

    /* Initialize virtual tables for solvers */
    if (user_opt->claw_version == 4)
    {
        fc2d_clawpack46_solver_initialize(glob);
    }
    else if (user_opt->claw_version == 5)
    {
        fc2d_clawpack5_solver_initialize(glob);
    }

    swirl_link_solvers(glob);
    swirl_initialize_rays(glob);

    /* ---------------------------------------------------------------
       Run
       --------------------------------------------------------------- */
    fclaw2d_initialize(glob);
    fclaw2d_run(glob);

    fclaw2d_finalize(glob);
}

int
main (int argc, char **argv)
{
    int first_arg;
    fclaw_app_t *app;
    fclaw_exit_type_t vexit;

    /* Options */
    user_options_t              *user_opt;
    fclaw_options_t             *fclaw_opt;
    fclaw_clawpatch_options_t *clawpatch_opt;
    fc2d_clawpack46_options_t   *claw46_opt;
    fc2d_clawpack5_options_t    *claw5_opt;

    fclaw_global_t            *glob;
    fclaw_domain_t            *domain;
    sc_MPI_Comm mpicomm;

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, NULL);

    /* Create new options packages */
    fclaw_opt =                   fclaw_options_register(app,  NULL,        "fclaw_options.ini");
    clawpatch_opt =   fclaw_clawpatch_options_register_2d(app, "clawpatch",  "fclaw_options.ini");
    claw46_opt =        fc2d_clawpack46_options_register(app, "clawpack46", "fclaw_options.ini");
    claw5_opt =          fc2d_clawpack5_options_register(app, "clawpack5",  "fclaw_options.ini");
    user_opt =                    swirl_options_register(app,               "fclaw_options.ini");

    /* Read configuration file(s) and command line, and process options */
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    /* Run the program */
    if (!vexit)
    {
        /* Options have been checked and are valid */

        mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);
        domain = create_domain(mpicomm, fclaw_opt);

        /* Create global structure which stores the domain, timers, etc */
        glob = fclaw_global_new();
        fclaw_global_store_domain(glob, domain);

        /* Store option packages in glob */
        fclaw_options_store           (glob, fclaw_opt);
        fclaw_clawpatch_options_store (glob, clawpatch_opt);
        fc2d_clawpack46_options_store   (glob, claw46_opt);
        fc2d_clawpack5_options_store    (glob, claw5_opt);
        swirl_options_store             (glob, user_opt);

        run_program(glob);
        fclaw_global_destroy(glob);
    }

    fclaw_app_destroy (app);

    return 0;
}
