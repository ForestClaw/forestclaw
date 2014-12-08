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
#include <fclaw2d_map_query.h>

#include <amr_forestclaw.H>
#include <amr_utils.H>
#include <fclaw_options.h>

#include <p4est_connectivity.h>

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
  int                      lp;
  sc_MPI_Comm              mpicomm;
  sc_options_t             *options;
  p4est_connectivity_t     *conn = NULL;
  fclaw2d_domain_t	   *domain;

  fclaw2d_map_context_t    *cont = NULL, *brick = NULL;
  amr_options_t             samr_options, *gparms = &samr_options;
  fclaw2d_clawpack_parms_t  sclawpack_parms, *clawpack_parms = &sclawpack_parms;
  fclaw2d_map_data_t        smap_data, *map_data = &smap_data;

  int example, retval;
  double pi = M_PI;

  /* Mapping variables */
  double rotate[2];
  double alpha, beta;  /* Ratio of torus inner radius to outer radius */
  const char* latitude_string, *longitude_string;
  double *latitude, *longitude;
  int mi, mj, a,b;

  lp = SC_LP_PRODUCTION;
  mpicomm = sc_MPI_COMM_WORLD;
  fclaw_mpi_init (&argc, &argv, mpicomm, lp);

  options = sc_options_new (argv[0]);



  /* [main] Register parameters local to the example.  This could be put into
     a function,  which might make things more readable.  One slight hassle is
     that the strings needed to read arrays have to then be passed back so they
     can be converted to numerical arrays after the file and command line have
     been read. */

  sc_options_add_int (options, 0, "main:example", &example, 0,
                      "[main] 1 = cart; 2 = torus; 3 = lat-long; 4 = annulus [2]");

  sc_options_add_double (options, 0, "main:alpha", &alpha, 0.4,
                         "[main] Ratio r/R, r=outer radius, R=inner radius " \
                         "(used for torus) [0.4]");

  fclaw_options_add_double_array(options, 0, "main:latitude", &latitude_string,
                                 "-50 50", &latitude, 2,
                                 "[main] Latitude range (degrees) [-50 50]");

  fclaw_options_add_double_array(options, 0, "main:longitude", &longitude_string,
                                 "0 360", &longitude, 2,
                                 "[main] Longitude range (degrees) [0 360]");

  sc_options_add_double (options, 0, "main:beta", &beta, 0.4,
                         "[main] Inner radius of annulus [0.4]");

  /* [Options] Register general ForestClaw options */
  fclaw_options_register(options,gparms);

  /* [mapping] Register general mapping data */
  fclaw2d_register_map_data(options,map_data);

  /* [clawpack46] Register solver options */
  clawpack46_register_options(options,clawpack_parms);

  /* Read options from fclaw_options.ini */
  retval = fclaw_options_read_from_file(options, lp);

  /* Parse command line */
  retval = retval || fclaw_options_parse_command_line (options,argc, argv, lp);

  /* post-process array input */
  fclaw_options_postprocess(gparms);
  clawpack46_postprocess_parms(clawpack_parms);
  fclaw2d_options_postprocess_map_data(map_data);

  if (example == 3)
  {
      fclaw_options_convert_double_array (latitude_string, &latitude,2);
      fclaw_options_convert_double_array (longitude_string, &longitude,2);
  }

  /* Check final state of parameters.  Return from help message, if necessary. */
  retval = retval || fclaw_options_check (options, gparms, lp);
  retval = retval || clawpack46_checkparms(options,clawpack_parms,gparms,lp);
  retval = retval || torus_checkparms (example, lp);  /* Nothing more to check here */
  if (!retval)
  {
      fclaw_options_print_summary(options,lp);

      if (gparms->trapfpe == 1)
      {
          printf("Enabling floating point traps\n");
          feenableexcept(FE_INVALID);
      }

      if (gparms->mpi_debug == 1)
      {
          fclaw2d_mpi_debug();
      }

      /* ---------------------------------------------------------------
         Mapping geometry
         --------------------------------------------------------------- */
      mi = map_data->mi;
      mj = map_data->mj;
      rotate[0] = pi*map_data->theta/180.0;
      rotate[1] = pi*map_data->phi/180.0;
      a = map_data->periodic_x;
      b = map_data->periodic_y;

      switch (example)
      {
      case 1:
          conn = p4est_connectivity_new_brick(mi,mj,a,b);
          brick = fclaw2d_map_new_brick(conn,mi,mj);
          cont = fclaw2d_map_new_cart(brick,map_data->scale,map_data->shift,rotate);
          break;
      case 2:
          conn = p4est_connectivity_new_brick(mi,mj,a,b);
          brick = fclaw2d_map_new_brick(conn,mi,mj);
          cont = fclaw2d_map_new_torus(brick,map_data->scale,map_data->shift,rotate,alpha);
          break;
      case 3:
          /* Lat-long example */
          conn = p4est_connectivity_new_brick(mi,mj,a,b);
          brick = fclaw2d_map_new_brick(conn,mi,mj);
          cont = fclaw2d_map_new_latlong(brick,map_data->scale,map_data->shift,
                                         rotate,latitude,longitude,a,b);
          break;
      case 4:
          /* Annulus */
          conn = p4est_connectivity_new_brick(mi,mj,a,b);
          brick = fclaw2d_map_new_brick(conn,mi,mj);
          cont = fclaw2d_map_new_annulus(brick,map_data->scale,map_data->shift,
                                         rotate,beta);
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

      fclaw2d_map_destroy(cont);
  }

  /* Destroy arrays used in options  */
  if (example == 3)
  {
      fclaw_options_destroy_array((void*) latitude);
      fclaw_options_destroy_array((void*) longitude);
  }

  fclaw2d_map_destroy_arrays(map_data);
  sc_options_destroy (options);
  fclaw_options_destroy_arrays (gparms);
  fclaw2d_clawpack_parms_delete(clawpack_parms);

  fclaw_mpi_finalize ();

  return 0;
}
