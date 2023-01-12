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

#include "filament/filament_user.h"
#include "swirl/swirl_user.h"
#include "user.h"

#include "../all/advection_user.h"

#include <fclaw_filesystem.h>
#include <fclaw_base.h>

#include <fclaw2d_forestclaw.h>
#include <fclaw2d_global.h>

static
void filament_initialize(fclaw2d_global_t* glob)
{
    fclaw2d_set_global_context(glob);

    filament_options_t             *user;

    user = (filament_options_t*) filament_get_options(glob);


    /* ---------------------------------------------------------------
       Set domain data.
       --------------------------------------------------------------- */
    fclaw2d_domain_data_new(glob->domain);

    /* Initialize virtual table for ForestClaw */
    fclaw2d_vtables_initialize(glob);

    if (user->claw_version == 4)
    {
      fc2d_clawpack46_solver_initialize(glob);
    }
    else if (user->claw_version == 5)
    {
      fc2d_clawpack5_solver_initialize(glob);
    }

    filament_link_solvers(glob);

    fclaw2d_initialize(glob);

    fclaw2d_clear_global_context(glob);
}

static
void filament_finalize(fclaw2d_global_t* glob)
{
    fclaw2d_set_global_context(glob);

    fclaw2d_problem_setup(glob);
    fclaw2d_finalize(glob);

    fclaw2d_clear_global_context(glob);
}
static
void swirl_initialize(fclaw2d_global_t* glob)
{
    fclaw2d_set_global_context(glob);

    const swirl_options_t           *swirl_opt;

    /* ---------------------------------------------------------------
       Set domain data.
       --------------------------------------------------------------- */
    fclaw2d_domain_data_new(glob->domain);

    swirl_opt = swirl_get_options(glob);

    /* Initialize virtual table for ForestClaw */
    fclaw2d_vtables_initialize(glob);

    /* Initialize virtual tables for solvers */
    if (swirl_opt->claw_version == 4)
    {
        fc2d_clawpack46_solver_initialize(glob);
    }
    else if (swirl_opt->claw_version == 5)
    {
        fc2d_clawpack5_solver_initialize(glob);
    }

    swirl_link_solvers(glob);

    /* ---------------------------------------------------------------
       Run
       --------------------------------------------------------------- */
    fclaw2d_initialize(glob);

    fclaw2d_clear_global_context(glob);
}
static
void swirl_finalize(fclaw2d_global_t* glob)
{
    fclaw2d_set_global_context(glob);

    fclaw2d_problem_setup(glob);
    fclaw2d_finalize(glob);

    fclaw2d_clear_global_context(glob);
}

typedef struct overlap_prodata
{
  double              myvalue;
  int                 isset;
}
overlap_prodata_t;

typedef struct overlap_point
{
  size_t              lnum;
  double              xy[2];
  overlap_prodata_t   prodata;
}
overlap_point_t;

typedef struct overlap_consumer
{
  fclaw2d_global_t   *glob;
  fclaw2d_domain_t   *domain;
  sc_array_t         *query_points;
  size_t              cell_idx;
  int                 num_cells_in_patch;
}
overlap_consumer_t;

static
void add_cell_centers (fclaw2d_domain_t * domain, fclaw2d_patch_t * patch,
                       int blockno, int patchno, void *user)
{
    overlap_point_t *op;
    int mx, my, mbc, i, j;
    double xlower, ylower, dx, dy;

    /* assert that user is a valid overlap_consumer_t */
    overlap_consumer_t *c = (overlap_consumer_t *) user;
    FCLAW_ASSERT (c != NULL);
    FCLAW_ASSERT (c->domain != NULL);
    FCLAW_ASSERT (c->query_points != NULL);

    /* create one query point for every cell in the patch */
    fclaw2d_clawpatch_grid_data (c->glob, patch, &mx, &my, &mbc,
                                 &xlower, &ylower, &dx, &dy);
    for (i = 0; i < mx; i++)
    {
        for (j = 0; j < my; j++)
        {
            /* initialize the query-point corresponding to this cell */
            op = (overlap_point_t *) sc_array_index (c->query_points,
                                                     c->cell_idx);
            op->lnum = c->cell_idx++;   /* local index of the cell */
            op->prodata.isset = 0;
            op->prodata.myvalue = -1;

            /* choose the middle point of the cell */
            op->xy[0] = xlower + (2 * i + 1) * dx / 2.;
            op->xy[1] = ylower + (2 * j + 1) * dy / 2.;

            /* Here, we would need to apply the mapping of the consumer domain
             * (in this example the swirl domain) to the point. */
        }
    }
}

static
void create_query_points (overlap_consumer_t * c)
{
    /* We create a process-local set of query points, for which we want to
     * obtain interpolation data from the producer side. We query the
     * center-point of every local cell. */
    c->query_points = sc_array_new_count (sizeof (overlap_point_t),
                                          c->domain->local_num_patches *
                                          c->num_cells_in_patch);
    c->cell_idx = 0;
    fclaw2d_domain_iterate_patches (c->domain, add_cell_centers, c);

    /* verify that we created as many query_points as expected */
    FCLAW_ASSERT (c->cell_idx ==
                  (size_t) c->domain->local_num_patches *
                  c->num_cells_in_patch);
}

static
int overlap_interpolate (fclaw2d_domain_t * domain, fclaw2d_patch_t * patch,
                         int blockno, int patchno, void *point, void *user)
{
    overlap_point_t *op;
    double tol;
    int consumer_side;

    /* assert that we got passed a valid overlap_point_t */
    FCLAW_ASSERT (point != NULL);
    op = (overlap_point_t *) point;

    /* check, if we are on the consumer or the producer side */
    consumer_side = domain_is_meta (domain);

    /* set tolerances */
    if (consumer_side)
    {
        /* we are on the consumer side and can only rely on basic domain
         * information. We do a stricter interpolation test in order to not lose
         * the accepted points on the producer side */
        tol = 0.5 * SC_1000_EPS;
    }
    else
    {
        tol = SC_1000_EPS;
    }

    /* We check if the query point intersects the patch.
     * It is important to take the underlying geometries of the consumer and
     * the producer side under consideration.
     * Here, we would need to apply the inverse mapping of the producer side
     * (in this example the filament side) to the query point first, so that it
     * lies in the same coordinate system as the patch-boundaries. */
    if ((op->xy[0] < patch->xlower - tol
         || op->xy[0] > patch->xupper + tol)
        || (op->xy[1] < patch->ylower - tol
            || op->xy[1] > patch->yupper + tol))
    {
        return 0;
    }

    fclaw_debugf ("Found point [%f,%f] in patch [%f,%f]x[%f,%f].\n",
                  op->xy[0], op->xy[1], patch->xlower, patch->xupper,
                  patch->ylower, patch->yupper);

    /* Although the point is located within a certain tolerance of the patch,
     * it may still lie outside of the [0,1]x[0,1]-block on which the domain is
     * defined. */
    op->xy[0] = SC_MAX (op->xy[0], 0.);
    op->xy[0] = SC_MIN (op->xy[0], 1.);
    op->xy[1] = SC_MAX (op->xy[1], 0.);
    op->xy[1] = SC_MIN (op->xy[1], 1.);

    /* update interpolation data */
    if (patchno >= 0)
    {
        /* we are on a leaf on the producer side */
        if (!op->prodata.isset)
        {
            /* We update the query point by setting its interpolation data.
             * We increment isset, to keep track of how many local patches
             * contributed to the points interpolation data.
             * We set myvalue to the mpirank of the domain for simplicity and
             * easy verification of the interpolation data at the end.
             * In practice, one should implement actual interpolation and store
             * the resulting data in the point-struct here. */
            op->prodata.isset++;
            op->prodata.myvalue = (double) domain->mpirank;
            fclaw_debugf
                ("Setting interpolation data of point [%f,%f] to %f.\n",
                 op->xy[0], op->xy[1], op->prodata.myvalue);
        }
    }

    return 1;
}

static
void output_query_points (overlap_consumer_t * c)
{
    size_t iz, npz;
    overlap_point_t *op;

    npz = c->query_points->elem_count;
    for (iz = 0; iz < npz; iz++)
    {
        op = (overlap_point_t *) sc_array_index (c->query_points, iz);
        printf
            ("Query point %ld on process %d is [%f,%f] and has interpolation data %f.\n",
             iz, c->domain->mpirank, op->xy[0], op->xy[1],
             op->prodata.myvalue);

    }
}

int
main (int argc, char **argv)
{
    fclaw_app_t *app;
    int first_arg;
    fclaw_exit_type_t vexit;
    overlap_consumer_t consumer, *c = &consumer;

    /* Options */
    sc_options_t                *options;

    filament_options_t          *filament_user_opt;
    fclaw_options_t             *filament_fclaw_opt;
    fclaw2d_clawpatch_options_t *filament_clawpatch_opt;
    fc2d_clawpack46_options_t   *filament_claw46_opt;
    fc2d_clawpack5_options_t    *filament_claw5_opt;

    swirl_options_t             *swirl_user_opt;
    fclaw_options_t             *swirl_fclaw_opt;
    fclaw2d_clawpatch_options_t *swirl_clawpatch_opt;
    fc2d_clawpack46_options_t   *swirl_claw46_opt;
    fc2d_clawpack5_options_t    *swirl_claw5_opt;

    fclaw2d_global_t         *filament_glob;
    fclaw2d_domain_t         *filament_domain;

    fclaw2d_global_t         *swirl_glob;
    fclaw2d_domain_t         *swirl_domain;

    sc_MPI_Comm mpicomm;

    int retval;

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, NULL);

    /* Register packages */
    filament_fclaw_opt                    = fclaw_options_register(app,  "filament",           "fclaw_options.ini");
    filament_clawpatch_opt    = fclaw2d_clawpatch_options_register(app, "filament-clawpatch",  "fclaw_options.ini");
    filament_claw46_opt         = fc2d_clawpack46_options_register(app, "filament-clawpack46", "fclaw_options.ini");
    filament_claw5_opt           = fc2d_clawpack5_options_register(app, "filament-clawpack5",  "fclaw_options.ini");
    filament_user_opt =                  filament_options_register(app, "filament-user",       "fclaw_options.ini");  

    swirl_fclaw_opt =                   fclaw_options_register(app, "swirl",            "fclaw_options.ini");
    swirl_clawpatch_opt =   fclaw2d_clawpatch_options_register(app, "swirl-clawpatch",  "fclaw_options.ini");
    swirl_claw46_opt =        fc2d_clawpack46_options_register(app, "swirl-clawpack46", "fclaw_options.ini");
    swirl_claw5_opt =          fc2d_clawpack5_options_register(app, "swirl-clawpack5",  "fclaw_options.ini");
    swirl_user_opt =                    swirl_options_register(app, "swirl-user",       "fclaw_options.ini");  

    /* Read configuration file(s) */
    options = fclaw_app_get_options (app);
    retval = fclaw_options_read_from_file(options);
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    if (!retval & !vexit)
    {
        /* Options have been checked and are valid */

        mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);

        /* Domains */
        filament_domain = filament_create_domain(mpicomm, filament_fclaw_opt, filament_user_opt,filament_clawpatch_opt);
        swirl_domain = swirl_create_domain(mpicomm, swirl_fclaw_opt);
            

        /* Globs */
        filament_glob = fclaw2d_global_new();
        fclaw2d_global_store_domain(filament_glob, filament_domain);

        fclaw2d_options_store            (filament_glob, filament_fclaw_opt);
        fclaw2d_clawpatch_options_store  (filament_glob, filament_clawpatch_opt);
        fc2d_clawpack46_options_store    (filament_glob, filament_claw46_opt);
        fc2d_clawpack5_options_store     (filament_glob, filament_claw5_opt);
        filament_options_store           (filament_glob, filament_user_opt);

        swirl_glob = fclaw2d_global_new();
        fclaw2d_global_store_domain(swirl_glob, swirl_domain);

        fclaw2d_options_store           (swirl_glob, swirl_fclaw_opt);
        fclaw2d_clawpatch_options_store (swirl_glob, swirl_clawpatch_opt);
        fc2d_clawpack46_options_store   (swirl_glob, swirl_claw46_opt);
        fc2d_clawpack5_options_store    (swirl_glob, swirl_claw5_opt);
        swirl_options_store             (swirl_glob, swirl_user_opt);

        /* initialize */
        filament_initialize(filament_glob);
        swirl_initialize(swirl_glob);

        /* compute process-local query points on the consumer side */
        c->glob = swirl_glob;
        c->domain = swirl_glob->domain;
        c->num_cells_in_patch =
            swirl_clawpatch_opt->mx * swirl_clawpatch_opt->my;
        create_query_points (c);

        /* obtain interpolation data of the points from the producer side */
        fclaw2d_overlap_exchange (filament_glob->domain, c->query_points,
                                  overlap_interpolate, NULL);

        /* output the interpolation data for all query points */
        output_query_points (c);

        /* run */
        fclaw2d_global_t *globs[2];
        globs[0] = filament_glob;
        globs[1] = swirl_glob;
        user_run (globs, 2);

        /* finalzie */
        filament_finalize(filament_glob);
        swirl_finalize(swirl_glob);


        /* destroy */
        sc_array_destroy (c->query_points);
        fclaw2d_global_destroy (filament_glob);
        fclaw2d_global_destroy (swirl_glob);

    }

    fclaw_app_destroy (app);

    return 0;
}
