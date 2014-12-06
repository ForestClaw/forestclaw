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

#include <fclaw2d_clawpack.H>
#include <fclaw2d_map.h>
#include <p4est_connectivity.h>

#include <amr_forestclaw.H>
#include <amr_utils.H>
#include <fclaw2d_options.h>
#include <fclaw2d_map_query.h>

#include "torus_user.H"

static int
torus_checkparms (int example, int lp)
{
    if (example < 1 || example > 4) {
        fclaw2d_global_log (lp, "Option --example must be 1, 2, 3 or 4\n");
        return -1;
    }

    return 0;
}

int
main (int argc, char **argv)
{
  int                   lp;
  int                   retval;
  sc_MPI_Comm           mpicomm;
  sc_options_t          *options;
  p4est_connectivity_t  *conn = NULL;
  fclaw2d_map_context_t *cont = NULL, *brick = NULL;
  fclaw2d_domain_t	*domain;
  amr_options_t         samr_options, *gparms = &samr_options;
  fclaw2d_clawpack_parms_t  sclawpack_parms, *clawpack_parms = &sclawpack_parms;

  int example;
  double pi = M_PI;

  /* Mapping variables */
  double rotate[2];
  double alpha;  /* Ratio of torus inner radius to outer radius */
  double lat[2];
  double longitude[2];
  int mi, mj, a,b;
  fclaw2d_map_data_t smap_data, *map_data = &smap_data;

  int trapfpe;
  int mpi_debug;

  lp = SC_LP_PRODUCTION;
  mpicomm = sc_MPI_COMM_WORLD;
  fclaw_mpi_init (&argc, &argv, mpicomm, lp);

  options = sc_options_new (argv[0]);

  /* Register local options */
  sc_options_add_int (options, 0, "example", &example, 0,
                      "[Example] 1 = cart; 2 = torus; 3 = lat-long; 4 = annulus");

  sc_options_add_switch (options, 0, "mpi_debug", &mpi_debug,
                      "[Init] Start MPI debug session to attach processes to gdb [F]");

  sc_options_add_switch (options, 0, "trapfpe", &trapfpe,
                      "[Init] Trap floating point exceptions [T]");

  /* [forestclaw] Register general ForestClaw options */
  fclaw2d_register_options(options,gparms);
  fclaw2d_read_options_from_file(options);
  fclaw2d_postprocess_parms(gparms);  /* post-process array input */

  /* [Mapping] Register general mapping data */
  fclaw2d_register_map_data(options,map_data); /* sets default values */
  fclaw2d_read_options_from_file(options);  /* Read options from fclaw2d_defaults.ini */
  mi = map_data->mi;
  mj = map_data->mj;
  rotate[0] = pi*map_data->theta/180.0;
  rotate[1] = pi*map_data->phi/180.0;

  /* [clawpack] Register solver options */
  clawpack46_register_options(options,clawpack_parms); /* Change name space? */
  clawpack46_read_options_from_file(options);
  clawpack46_postprocess_parms(clawpack_parms);  /* post-process array input */

  /* TODO: would it make sense to provide this code up to and excluding
   *       torus_checkparms into a reusable function? */
  /* Hm... maybe.  But such a reusable function would have to be specific for
     each example, I think. And it isn't clear how it would be get re-used. */

  /* This prints out all the current options */
  fclaw2d_parse_command_line (options,argc, argv, lp);

  fclaw2d_postprocess_parms(gparms); /* post-process array input */
  clawpack46_postprocess_parms(clawpack_parms); /* post-process array input */

  /* Check final state of parameters.
   * The call to amr_checkparms2 also checks for a --help message.
   * We are exploiting C short-circuit boolean evaluation.
   */
  retval = fclaw2d_checkparms (options, gparms, lp);
  retval = retval || clawpack46_checkparms(options,clawpack_parms,gparms,lp);
  retval = retval || torus_checkparms (example, lp);  /* No more to check here */
  if (!retval) {
     /* the do-the-work block. TODO: put everything below into a function */

      if (trapfpe == 1)
      {
          printf("Enabling floating point traps\n");
          feenableexcept(FE_INVALID);
      }

      if (mpi_debug == 1)
      {
          fclaw2d_mpi_debug();
      }

      /* ---------------------------------------------------------------
         Mapping geometry
         --------------------------------------------------------------- */
      switch (example)
      {
      case 1:
          a = 1;
          b = 1;
          conn = p4est_connectivity_new_brick(mi,mj,a,b);
          brick = fclaw2d_map_new_brick(conn,mi,mj);
          cont = fclaw2d_map_new_cart(brick,map_data->scale,map_data->shift,rotate);
          break;
      case 2:
          /* Generally, we should have mj \approx \alpha*mi */
          /* where alpha is the ratio of inner radius to outer radius */
          a = 1;
          b = 1;
          alpha = 0.4;
          mj = alpha*mi;
          if (mj == 0)
          {
              mi = 1;
              mj = 1;
          }
          conn = p4est_connectivity_new_brick(mi,mj,a,b);
          brick = fclaw2d_map_new_brick(conn,mi,mj);
          cont = fclaw2d_map_new_torus(brick,map_data->scale,map_data->shift,rotate,alpha);
          break;
      case 3:
          /* Lat-long example */
          a = 1;
          b = 0;
          longitude[0] = 0; /* x-coordinate */
          longitude[1] = 360;  /* if a == 1, long[1] will be computed as [0] + 180 */
          lat[0] = -50;  /* y-coordinate */
          lat[1] = 50;
          alpha = (lat[1]-lat[0])/180;
          mj = alpha*mi/2.0;
          if (mj == 0)
          {
              mi = 1;
              mj = 1;
          }
          conn = p4est_connectivity_new_brick(mi,mj,a,b);
          brick = fclaw2d_map_new_brick(conn,mi,mj);
          cont = fclaw2d_map_new_latlong(brick,map_data->scale,map_data->shift,
                                         rotate,lat,longitude,a,b);
          break;
      case 4:
          /* Annulus */
          a = 1;
          b = 0;
          alpha = 0.4;  /* Inner radius */
          mj = (1-alpha)/(1+alpha)*mi/pi;
          if (mj == 0)
          {
              mi = 1;
              mj = 1;
          }
          conn = p4est_connectivity_new_brick(mi,mj,a,b);
          brick = fclaw2d_map_new_brick(conn,mi,mj);
          cont = fclaw2d_map_new_annulus(brick,map_data->scale,map_data->shift,
                                         rotate,alpha);
          break;

      default:
          SC_ABORT_NOT_REACHED (); /* must be checked in torus_checkparms */
      }

      domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);


      /* ---------------------------------------------------------- */
      if (gparms->verbosity > 0)
      {
          fclaw2d_domain_list_levels(domain, lp);
          fclaw2d_domain_list_neighbors(domain, lp);
      }

      /* ---------------------------------------------------------------
         Set domain data.
         --------------------------------------------------------------- */
      init_domain_data(domain);

      set_domain_parms(domain,gparms);
      set_clawpack_parms(domain,clawpack_parms);

      link_problem_setup(domain,torus_problem_setup);

      torus_link_solvers(domain);

      link_regrid_functions(domain,
                            torus_patch_tag4refinement,
                            torus_patch_tag4coarsening);
      amrinit(&domain);
      amrrun(&domain);
      amrreset(&domain);

      /* Destroy mapping functions */
      fclaw2d_map_destroy(cont);
  } /* this bracket ends the do_the_work block */

  sc_options_destroy (options);         /* this could be moved up */
  fclaw2d_options_destroy (gparms);
  fclaw2d_clawpack_parms_delete(clawpack_parms);

  fclaw_mpi_finalize ();

  return 0;
}
