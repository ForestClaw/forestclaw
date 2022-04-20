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

typedef enum {
  SWIRL_RAY_LINE,
  SWIRL_RAY_CIRCLE,
  SWIRL_RAY_TYPE_LAST
}
swirl_ray_type_t;

typedef struct swirl_ray
{
  swirl_ray_type_t    rtype;
  double              xy[2];
  union {
    struct {
      double              vec[2];
    } line;
    struct {
      double              radius;
    } circle;
  } r;
}
swirl_ray_t;

static int
intersect_ray (fclaw2d_domain_t *domain, fclaw2d_patch_t * patch,
               int blockno, int patchno, void *ray, double *integral)
{
  swirl_ray_t *swirl_ray;

  /* assert that ray is a valid swirl_ray_t */
  swirl_ray = (swirl_ray_t *) ray;
  FCLAW_ASSERT(swirl_ray != NULL);

  if(patchno >= 0) {
    /* We are at a leaf and the patch is a valid patch of the domain.
     * Based on the patch, the domain, the blockno and the information stored
     * in the swirl_ray_t we defined, we now have to set *integral to be the
     * contribution of this ray-patch combination to the ray integral.
     * Additionaly, we should return 1, if ray and patch do intersect,
     * and 0 otherwise. */
    *integral = swirl_ray->r.line.vec[0];
    return 1;
  } else {
    /* We are not at a leaf and the patch is an artificial patch containing all
     * standard patch information except for the pointer to the next patch and
     * user-data of any kind. Only the FCLAW2D_PATCH_CHILDID and the
     * FCLAW2D_PATCH_ON_BLOCK_FACE_* flags are set.
     * Based on this, we now can run a test to check if the ray and the patch
     * intersect.
     * We return 0, if we know that the ray does not intersect any descendant
     * of this patch.
     * We return 1, if the test concludes that the ray may intersect the patch.
     * This test may be overinclusive.
     * Its purpose is to remove irrelevant patch-ray-combinations early on to
     * avoid unnecessary computations. */
    return 1;
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
    /* We integrate after the run of the solver finished.
     * It may be of interest to integrate several times during the run for
     * different time steps */
    fclaw2d_domain_integrate_rays(glob->domain, intersect_ray, rays, integrals);
    fclaw2d_finalize(glob);
}

static int           nlines = 3;

static sc_array_t *
swirl_rays_new (void)
{
  int                 i;
  swirl_ray_t        *ray;
  sc_array_t         *a = sc_array_new (sizeof (swirl_ray_t));

  /* add a couple straight rays */
  for (i = 0; i < nlines; ++i) {
    ray = (swirl_ray_t *) sc_array_push (a);
    ray->rtype = SWIRL_RAY_LINE;
    ray->xy[0] = 0.;
    ray->xy[1] = 0.;
    ray->r.line.vec[0] = cos (i * M_PI / nlines);
    ray->r.line.vec[1] = sin (i * M_PI / nlines);
  }

  /* add no circles yet */
  return a;
}

static sc_array_t *
swirl_integrals_new(void)
{
  sc_array_t         *a = sc_array_new (sizeof (double));
  sc_array_resize (a, nlines);
  return a;
}

int
main (int argc, char **argv)
{
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
    rays = swirl_rays_new ();
    integrals = swirl_integrals_new();

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
