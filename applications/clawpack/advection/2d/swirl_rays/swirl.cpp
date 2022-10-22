/*
  Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#if 0
const int swirl_nlines = 3;
#else
const int swirl_nlines = 8;
#endif

/* Virtual function for setting rays */
static
void swirl_allocate_and_define_rays(fclaw2d_global_t *glob,
                                    fclaw2d_ray_t** rays,
                                    int* num_rays)
{
    *num_rays = swirl_nlines;

    /* We let the user allocate an array of rays, although what is inside the
       generic ray type is left opaque. This is destroy in matching FREE,
       below. */

    *rays = fclaw2d_ray_allocate_rays(*num_rays);
    fclaw2d_ray_t *ray_vec = *rays;
    for (int i = 0; i < swirl_nlines; ++i)
    {
        swirl_ray_t *sr = (swirl_ray_t*) FCLAW_ALLOC(swirl_ray_t,1);
        sr->rtype = SWIRL_RAY_LINE;

#if 0
        /* End points are on a semi-circle in x>0,y>0 quad */
        FCLAW_ASSERT(swirl_nlines >= 2);
        sr->xy[0] = 0; //-0.1;
        sr->xy[1] = 0; //-0.1;
        double R = 2.0;
        double dth = M_PI/(2*swirl_nlines);
        sr->r.line.vec[0] = R*cos ((i+0.5) * dth);
        sr->r.line.vec[1] = R*sin ((i+0.5) * dth);
#else
        sr->xy[0] = 0.5;
        sr->xy[1] = 0.5;
        /* we add 0.1, since the intersection callback does not guarantee exact
         * results for axis-parallel rays */
        sr->r.line.vec[0] = cos ((i + 0.1) * 2 * M_PI / swirl_nlines);
        sr->r.line.vec[1] = sin ((i + 0.1) * M_PI / swirl_nlines);
#endif

        fclaw2d_ray_t *ray = &ray_vec[i];
        int id = i + 1;
        fclaw2d_ray_set_ray(ray,id, sr);
    }
}

static
void swirl_deallocate_rays(fclaw2d_global_t *glob,
                           fclaw2d_ray_t** rays,
                           int* num_rays)
{
    fclaw2d_ray_t *ray_vec = *rays;
    for(int i = 0; i < *num_rays; i++)
    {
        /* Retrieve rays set above and deallocate them */
        int id;
        fclaw2d_ray_t *ray = &ray_vec[i];
        swirl_ray_t *rs = (swirl_ray_t*) fclaw2d_ray_get_ray(ray,&id);
        FCLAW_ASSERT(rs != NULL);
        FCLAW_FREE(rs);
        rs = NULL;
    }
    /* Match FCLAW_ALLOC, above */
    *num_rays = fclaw2d_ray_deallocate_rays(rays);
}

static int
swirl_intersect_ray (fclaw2d_domain_t * domain, fclaw2d_patch_t * patch,
               int blockno, int patchno, void *ray, double *integral,
               void *user)
{
    int i, ni;
    double corners[2][2];

    /* assert that ray is a valid swirl_ray_t */
    fclaw2d_ray_t *fclaw_ray = (fclaw2d_ray_t *) ray;

    int id;
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

    if (fabs (swirl_ray->r.line.vec[0]) <= 1e-12 ||
        fabs (swirl_ray->r.line.vec[1]) <= 1e-12)
    {
        /* we cannot guarantee correct results for rays
         * that run near parallel to coordinate axis */
        return 0;
    }

    /* store the patch corners in an indexable format */
    corners[0][0] = patch->xlower;
    corners[0][1] = patch->ylower;
    corners[1][0] = patch->xupper;
    corners[1][1] = patch->yupper;

    /* for stability we search in the dimension of the strongest component */
    i = (fabs (swirl_ray->r.line.vec[0]) <= fabs (swirl_ray->r.line.vec[1])) ? 1 : 0;
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
        t = (corners[0][i] - swirl_ray->xy[i]) / swirl_ray->r.line.vec[i];
        shift = swirl_ray->xy[ni] + t * swirl_ray->r.line.vec[ni];

        /* shift coordinate system to the first hit */
        hits[0][0] = 0.;
        hits[0][1] = 0.;
        corners[1][ni] -= shift;
        corners[0][ni] -= shift;

        /* compute second hit in shifted coordinate system */
        hits[1][i] = corners[1][i] - corners[0][i];
        t = hits[1][i] / swirl_ray->r.line.vec[i];
        hits[1][ni] = t * swirl_ray->r.line.vec[ni];

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
            t = (corners[j][i] - swirl_ray->xy[i]) / swirl_ray->r.line.vec[i];
            rayatt = swirl_ray->xy[ni] + t * swirl_ray->r.line.vec[ni];

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

void swirl_initialize_rays(fclaw2d_global_t* glob)
{
    /* Set up rays */
    fclaw2d_ray_vtable_t* rays_vt = fclaw2d_ray_vt(glob);

    rays_vt->allocate_and_define = swirl_allocate_and_define_rays;
    rays_vt->deallocate = swirl_deallocate_rays;

    rays_vt->integrate = swirl_intersect_ray;
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
void run_program(fclaw2d_global_t* glob)
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
    fclaw_opt =                   fclaw_options_register(app,  NULL,        "fclaw_options.ini");
    clawpatch_opt =   fclaw2d_clawpatch_options_register(app, "clawpatch",  "fclaw_options.ini");
    claw46_opt =        fc2d_clawpack46_options_register(app, "clawpack46", "fclaw_options.ini");
    claw5_opt =          fc2d_clawpack5_options_register(app, "clawpack5",  "fclaw_options.ini");
    user_opt =                    swirl_options_register(app,               "fclaw_options.ini");  

    /* Read configuration file(s) and command line, and process options */
    options = fclaw_app_get_options (app);
    retval = fclaw_options_read_from_file(options);
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

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

        run_program(glob);
        fclaw2d_global_destroy(glob);
    }

    fclaw_app_destroy (app);

    return 0;
}
