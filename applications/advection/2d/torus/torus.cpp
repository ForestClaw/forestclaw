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

#include <forestclaw2d.h>
#include <fclaw_options.h>

#include <fclaw2d_clawpack.H>
#include <fclaw2d_map.h>
#include <fclaw2d_map_query.h>

#include <amr_forestclaw.H>
#include <amr_utils.H>
#include <fclaw_options.h>

#include <p4est_connectivity.h>

#include "torus_user.H"

static int
torus_checkparms (int example)
{
    if (example < 1 || example > 4) {
        fclaw_global_essentialf ("Option --example must be 1, 2, 3 or 4\n");
        return -1;
    }
    return 0;
}

typedef struct user_options
{
    int example;
    double alpha;

    const char* latitude_string;
    double *latitude;

    const char* longitude_string;
    double *longitude;

    double beta;
} user_options_t;



int
main (int argc, char **argv)
{
  fclaw_app_t *app;
  int first_arg;
  fclaw_exit_type_t vexit;

  sc_MPI_Comm              mpicomm;
  sc_options_t             *options;
  p4est_connectivity_t     *conn = NULL;
  fclaw2d_domain_t	   *domain;

  /* ForestClaw options */
  amr_options_t            samr_options, *gparms = &samr_options;

  /* Constants */
  double pi = M_PI;

  /* Mapping  */
  fclaw2d_map_context_t    *cont = NULL, *brick = NULL;

  /* Clawpack options */
  fclaw2d_clawpack_parms_t  sclawpack_parms, *clawpack_parms = &sclawpack_parms;

  /* User variables */
  user_options_t suser_options, *user = &suser_options;

  double rotate[2];
  int mi, mj, a,b;
  int retval;

  /* -------------------------------------------------------------
     - Initialize application
     ------------------------------------------------------------- */
  app = fclaw_app_new (&argc, &argv, user);
  options = fclaw_app_get_options (app);
  mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);


  /* -------------------------------------------------------------
     - Add options to be read later, either from a .ini file or
     from the command line.
     - Order in which options are added is the order in which they
     will be printed out in help message.
     ------------------------------------------------------------- */

  /*  Register core options, including verbosity level. */
  fclaw_app_options_register_general (app, "fclaw_options.ini", gparms);

  /* [clawpack46] Add solver options */
  clawpack46_options_add(options,clawpack_parms);

  /* [main] User options in */
  sc_options_add_int (options, 0, "main:example", &user->example, 0,
                      "[main] 1 = cart; 2 = torus; 3 = lat-long; 4 = annulus [2]");

  sc_options_add_double (options, 0, "main:alpha", &user->alpha, 0.4,
                         "[main] Ratio r/R, r=outer radius, R=inner radius " \
                         "(used for torus) [0.4]");

  fclaw_options_add_double_array(options, 0, "main:latitude", &user->latitude_string,
                                 "-50 50", &user->latitude, 2,
                                 "[main] Latitude range (degrees) [-50 50]");

  fclaw_options_add_double_array(options, 0, "main:longitude", &user->longitude_string,
                                 "0 360", &user->longitude, 2,
                                 "[main] Longitude range (degrees) [0 360]");

  sc_options_add_double (options, 0, "main:beta", &user->beta, 0.4,
                         "[main] Inner radius of annulus [0.4]");

  /* -------------------------------------------------------------
     - Read options from fclaw_options.ini
     - Register core options (?)
     - Parse command line
     - postprocess array input
     - checkparms
     ------------------------------------------------------------- */

  /* Read fclaw_options.ini */
  retval = fclaw_options_read_from_file(options);

  /* Parse command line and post-process */
  vexit = fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

  clawpack46_postprocess_parms(clawpack_parms);

  if (user->example == 3)
  {
      fclaw_options_convert_double_array (user->latitude_string, &user->latitude,2);
      fclaw_options_convert_double_array (user->longitude_string, &user->longitude,2);
  }

  /* retval = retval || fclaw_options_check (options, gparms); */
  retval = retval || clawpack46_checkparms(options,clawpack_parms,gparms);
  retval = retval || torus_checkparms (user->example);  /* Nothing more to check here */
  /* -------------------------------------------------------------
     - Run program
     ------------------------------------------------------------- */
  if (!retval & !vexit)
  {
      /* Only print options if verbosity >= info */
      fclaw_options_print_summary(options);

      if (gparms->trapfpe == 1)
      {
          fclaw_global_infof("Enabling floating point traps\n");
          feenableexcept(FE_INVALID);
      }

      if (gparms->mpi_debug == 1)
      {
          fclaw2d_mpi_debug();
      }

      /* ---------------------------------------------------------------
         Mapping geometry
         --------------------------------------------------------------- */
      mi = gparms->mi;
      mj = gparms->mj;
      rotate[0] = pi*gparms->theta/180.0;
      rotate[1] = pi*gparms->phi/180.0;
      a = gparms->periodic_x;
      b = gparms->periodic_y;

      switch (user->example)
      {
      case 1:
          conn = p4est_connectivity_new_brick(mi,mj,a,b);
          brick = fclaw2d_map_new_brick(conn,mi,mj);
          cont = fclaw2d_map_new_cart(brick,gparms->scale,gparms->shift,rotate);
          break;
      case 2:
          conn = p4est_connectivity_new_brick(mi,mj,a,b);
          brick = fclaw2d_map_new_brick(conn,mi,mj);
          cont = fclaw2d_map_new_torus(brick,gparms->scale,gparms->shift,rotate,user->alpha);
          break;
      case 3:
          /* Lat-long example */
          conn = p4est_connectivity_new_brick(mi,mj,a,b);
          brick = fclaw2d_map_new_brick(conn,mi,mj);
          cont = fclaw2d_map_new_latlong(brick,gparms->scale,gparms->shift,
                                         rotate,user->latitude,user->longitude,a,b);
          break;
      case 4:
          /* Annulus */
          conn = p4est_connectivity_new_brick(mi,mj,a,b);
          brick = fclaw2d_map_new_brick(conn,mi,mj);
          cont = fclaw2d_map_new_annulus(brick,gparms->scale,gparms->shift,
                                         rotate,user->beta);
          break;

      default:
          SC_ABORT_NOT_REACHED (); /* must be checked in torus_checkparms */
      }

      domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);


      /* ---------------------------------------------------------- */
#if 0
      /* TODO : Replace this with an updated version? */
      /* How do I set lp?  Using macros? */
      if (gparms->verbosity > 0)
      {
          fclaw2d_domain_list_levels(domain, lp);
          fclaw2d_domain_list_neighbors(domain, lp);
      }
#endif

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
  if (user->example == 3)
  {
      fclaw_options_destroy_array((void*) user->latitude);
      fclaw_options_destroy_array((void*) user->longitude);
  }

  /* Which of these do I still need? */
  fclaw2d_clawpack_parms_delete(clawpack_parms);

  fclaw_app_destroy (app);

  return 0;
}
