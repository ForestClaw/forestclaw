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

#include "fclaw2d_capi.h"
#include <forestclaw2d.h>
#include <fclaw_options.h>
#include <fclaw_mpi.h>

#include "fp_exception_glibc_extension.h"


typedef struct fclaw_options_general
{
    /* This will eventually be replaced with particular members of amropt, but
       we leave thsi here for backwards compatibility */
    amr_options_t *amropt;

    /* this is just for ForestClaw debugging, no need to adopt elsewhere */
    int is_registered;     /**< Internal variable to double-check the flow */
}
fclaw_options_general_t;

static void *
options_register_general (fclaw_app_t * a, void *package, sc_options_t * opt)
{
    fclaw_options_general_t *general = (fclaw_options_general_t *) package;

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (opt != NULL);

    /* allocated storage for this package's option values */
    FCLAW_ASSERT (general != NULL);

    fclaw_options_add_general (opt, general->amropt);

    /* we do not need to work with the return value */
    general->is_registered = 1;
    return NULL;
}

fclaw_exit_type_t
fclaw_options_postprocess (fclaw_options_t * amropt)
{
  amropt->mthbc = NULL;
  fclaw_options_convert_int_array (amropt->mthbc_string, &amropt->mthbc,
                                   fclaw2d_NumFaces);

  amropt->scale = NULL;
  fclaw_options_convert_double_array (amropt->scale_string, &amropt->scale, 3);

  amropt->shift = NULL;
  fclaw_options_convert_double_array (amropt->shift_string, &amropt->shift, 3);

  return FCLAW_NOEXIT;
}

static fclaw_exit_type_t
options_postprocess_general (fclaw_app_t * a, void *package, void *registered)
{

    fclaw_options_general_t *general = (fclaw_options_general_t *) package;

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    /* errors from the key-value options would have showed up in parsing */

    /* postprocess this package */
    FCLAW_ASSERT (general != NULL);
    FCLAW_ASSERT (general->is_registered);

    /* Convert strings to arrays */
    return fclaw_options_postprocess (general->amropt);
}

fclaw_exit_type_t
fclaw_options_check (fclaw_options_t * gparms)
{
    /* Check outstyle. */
    if (gparms->outstyle == 1 && gparms->use_fixed_dt)
    {
        double dT_outer = gparms->tfinal / gparms->nout;
        double dT_inner = gparms->initial_dt;
        int nsteps = (int) floor (dT_outer / dT_inner + .5);
        if (fabs (nsteps * dT_inner - dT_outer) > 1e-8)
        {
            fclaw_global_essentialf
              ("For fixed dt, initial time step size must divide"
               " tfinal/nout exactly.\n");
            return FCLAW_EXIT_ERROR;
        }
    }

    /* Could also do basic sanity checks on mx,my,... */

    if (gparms->mpi_debug)
    {
        fclaw_global_infof("Entering mpi_debug session");
        fclaw_mpi_debug ();
    }

    if (gparms->trapfpe)
    {
        fclaw_global_infof("Enabling floating point traps\n");
        feenableexcept(FE_INVALID);
    }

    return FCLAW_NOEXIT;
}

fclaw_exit_type_t
options_check_general (fclaw_app_t * app, void *package, void *registered)
{
    fclaw_options_general_t *general = (fclaw_options_general_t *) package;

    return fclaw_options_check (general->amropt);
}

void
fclaw_options_reset (fclaw_options_t * amropt)
{
    FCLAW_FREE (amropt->mthbc);
    FCLAW_FREE (amropt->scale);
    FCLAW_FREE (amropt->shift);
}

static void
options_destroy_general (fclaw_app_t * a, void *package, void *registered)
{
    fclaw_options_general_t *general = (fclaw_options_general_t *) package;
    amr_options_t *amropt;

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    /* free this package */
    FCLAW_ASSERT (general != NULL);
    FCLAW_ASSERT (general->is_registered);

    /* Destroy option arrays created in postprocess */
    amropt = general->amropt;

    fclaw_options_reset (amropt);
    FCLAW_FREE(amropt);

    FCLAW_FREE (general);
}

static const fclaw_app_options_vtable_t options_vtable_general = {
    options_register_general,
    options_postprocess_general,
    options_check_general,
    options_destroy_general
};

amr_options_t* fclaw_options_register_general (fclaw_app_t * a, const char *configfile)
{
    fclaw_options_general_t* general;
    amr_options_t* gparms;
    FCLAW_ASSERT (a != NULL);

    /* Basic options for print out help message, current options, etc */
    fclaw_app_options_register_core (a, configfile);

    /* allocate storage for fclaw_options */
    /* we will free it in the options_destroy callback */
    gparms = FCLAW_ALLOC(amr_options_t,1);
    general = FCLAW_ALLOC_ZERO (fclaw_options_general_t, 1);

    general->amropt = gparms;

    /* Could also pass in a section header (set to NULL for now) */
    fclaw_app_options_register (a,NULL,
                                configfile,
                                &options_vtable_general,
                                general);
    return gparms;
}

void fclaw_options_add_general (sc_options_t * opt, amr_options_t* amropt)
{
    sc_options_add_int (opt, 0, "mx", &amropt->mx, 8,
                        "Number of grid cells per patch in x [8]");

    sc_options_add_int (opt, 0, "my", &amropt->my, 8,
                        "Number of grid cells per patch in y [8]");

    sc_options_add_double (opt, 0, "initial_dt", &amropt->initial_dt, 0.1,
                           "Initial time step size [0.1]");

    sc_options_add_int (opt, 0, "outstyle", &amropt->outstyle, 1,
                        "Output style (1,2,3) [1]");

    /* If outstyle == 1 */
    sc_options_add_double (opt, 0, "tfinal", &amropt->tfinal, 1.0,
                           "Final time [1.0]");

    sc_options_add_int (opt, 0, "nout", &amropt->nout, 10,
                        "Number of time steps, used with outstyle=3 [10]");

    /* Only needed if outstyle == 3 */
    sc_options_add_int (opt, 0, "nstep", &amropt->nstep, 1,
                        "Steps between output files, used with outstyle=3 [1]");


    /* This is a hack to control the VTK output while still in development.
     * The values are numbers which can be bitwise-or'd together.
     * 0 - no VTK output ever.
     * 1 - output for all stages of amrinit.
     * 2 - output whenever amr_output() is called.
     */
    sc_options_add_int (opt, 0, "vtkout", &amropt->vtkout, 0,
                        "VTK output method [F]");
    sc_options_add_double (opt, 0, "vtkspace", &amropt->vtkspace, 0.,
                           "VTK visual spacing [F]");
    sc_options_add_int (opt, 0, "vtkwrite", &amropt->vtkwrite, 0,
                        "VTK write variant [F]");

    /* output options */

#if 0
    /* This is handled now by core options (fclaw_base.c) */
    sc_options_add_int (opt, 0, "verbosity", &amropt->verbosity, 0,
                        "Verbosity mode [0]");
#endif

    sc_options_add_bool (opt, 0, "serialout", &amropt->serialout, 1,
                            "Enable serial output [F]");
    sc_options_add_string (opt, 0, "prefix", &amropt->prefix, "fort",
                           "Output file prefix [fort]");

    /* more clawpack options */
    sc_options_add_double (opt, 0, "max_cfl", &amropt->max_cfl, 1,
                           "Maximum CFL allowed [1]");

    sc_options_add_double (opt, 0, "desired_cfl", &amropt->desired_cfl, 0.9,
                           "Maximum CFL allowed [0.9]");

    sc_options_add_int (opt, 0, "meqn", &amropt->meqn, 1,
                        "Number of equations [1]");

    sc_options_add_int (opt, 0, "mbc", &amropt->mbc, 2,
                        "Number of ghost cells [2]");

    /* Array of NumFaces many values */
    fclaw_options_add_int_array (opt, 0, "mthbc", &amropt->mthbc_string, "1 1 1 1",
        &amropt->mthbc, fclaw2d_NumFaces,
        "Physical boundary condition type [1 1 1 1]");


    /* At this point amropt->mthbc is allocated. Set defaults if desired. */

    sc_options_add_int (opt, 0, "minlevel", &amropt->minlevel, 0,
                        "Minimum refinement level [0]");

    sc_options_add_int (opt, 0, "maxlevel", &amropt->maxlevel, 0,
                        "Maximum refinement level[0]");

    sc_options_add_int (opt, 0, "regrid_interval", &amropt->regrid_interval,
                        1, "Regridding frequency [1]");

    sc_options_add_int (opt, 0, "refratio", &amropt->refratio,
                        2, "Refinement ratio [2]");

    sc_options_add_bool (opt, 0, "smooth-refine", &amropt->smooth_refine,
                         0, "Refinement smoothing[F]");

    sc_options_add_int (opt, 0, "coarsen-delay", &amropt->coarsen_delay,
                        0, "Number skipped coarsenings[0]");

    sc_options_add_double (opt, 0, "refine_threshold", &amropt->refine_threshold,
                           0.5, "Refinement threshold [0.5]");

    sc_options_add_double (opt, 0, "coarsen_threshold", &amropt->coarsen_threshold,
                           0.1, "Coarsening threshold [0.1]");


    sc_options_add_double (opt, 0, "ax", &amropt->ax, 0, "xlower " \
                           "(used only with manifold=0) [0]");
    sc_options_add_double (opt, 0, "bx", &amropt->bx, 1, "xupper " \
                           "(used only with manifold=0) [1]");
    sc_options_add_double (opt, 0, "ay", &amropt->ay, 0, "ylower " \
                           "(used only with manifold=0) [0]");
    sc_options_add_double (opt, 0, "by", &amropt->by, 1, "yupper " \
                           "(used only with manifold=0) [1]");

    sc_options_add_bool (opt, 0, "manifold", &amropt->manifold, 0,
                         "Solution is on manifold [F]");

    sc_options_add_bool (opt, 0, "use_fixed_dt", &amropt->use_fixed_dt, 0,
                         "Use fixed coarse grid time step [F]");

    sc_options_add_bool (opt, 0, "run_diagnostics",
                         &amropt->run_diagnostics,0,
                         "Run diagnostics [F]");

    sc_options_add_bool (opt, 0, "conservation_check",
                         &amropt->conservation_check,0,
                         "Conservation check [F]");

    sc_options_add_bool (opt, 0, "subcycle", &amropt->subcycle, 1,
                         "Use subcycling in time [F]");

    sc_options_add_bool (opt, 0, "noweightedp", &amropt->noweightedp, 0,
                           "No weighting when subcycling [F]");

    /* ---------------------- Control execution -------------------------- */
    sc_options_add_bool (opt, 0, "trapfpe", &amropt->trapfpe,1,
                         "Trap floating point exceptions [1]");

    sc_options_add_bool (opt, 0, "mpi_debug", &amropt->mpi_debug,0,
                        "Start MPI debug session (for attaching processes in gdb) [0]");


    /* ---------------------- Mapping options -------------------------- */
    sc_options_add_int (opt, 0, "mi", &amropt->mi, 1,
                        "Number of blocks in x direction [1]");

    sc_options_add_int (opt, 0, "mj", &amropt->mj, 1,
                        "Number of blocks in y direction  [1]");

    sc_options_add_bool (opt, 0, "periodic_x", &amropt->periodic_x, 0,
                        "Periodic in x direction [F]");

    sc_options_add_bool (opt, 0, "periodic_y", &amropt->periodic_y, 0,
                        "Periodic in y direction  [F]");

    fclaw_options_add_double_array (opt,0, "scale",
                                    &amropt->scale_string, "1 1 1",
                                    &amropt->scale, 3,
                                    "Scale factor [1 1 1]");

    fclaw_options_add_double_array (opt,0, "shift",
                                    &amropt->shift_string, "0 0 0",
                                    &amropt->shift, 3,
                                    "Shift array [0 0 0]");

    sc_options_add_double (opt, 0, "phi", &amropt->phi, 0,
                           "Rotation angle about x axis (degrees) [0]");

    sc_options_add_double (opt, 0, "theta", &amropt->theta, 0,
                           "Rotation angle about z axis (degrees) [0]");

    /* -----------------------------------------------------------------------
       Options will be read from this file, if a '-F' flag is used at the command
       line.  Use this file for local modifications that are not tracked by Git.
       ----------------------------------------------------------------------- */
    sc_options_add_inifile (opt, 'F', "inifile",
                            "File used to override one or more options " \
                            "in fclaw_options.ini [empty]");

}

/* This is still needed as long as the config file option in app isn't yet implemented */
int fclaw_options_read_from_file(sc_options_t* opt)
{
    int retval;

    int fclaw_package_id;
    fclaw_package_id = fclaw_get_package_id ();
    retval = sc_options_load (fclaw_package_id, FCLAW_VERBOSITY_ESSENTIAL, opt,
                              "fclaw_options.ini");
    if (retval < 0)
    {
        fclaw_global_essentialf( \
                            "Problem reading fclaw_options.ini.\n");
    }
    else
    {
        fclaw_global_infof ("Reading file fclaw_options.ini.\n");
    }
    return retval;
}


void
fclaw_options_add_int_array (sc_options_t * opt,
                             int opt_char, const char *opt_name,
                             const char **array_string,
                             const char *default_string,
                             int **int_array, int initial_length,
                             const char *help_string)
{
    *int_array = NULL;
    sc_options_add_string (opt, opt_char, opt_name,
                           array_string, default_string, help_string);
}

void
fclaw_options_add_double_array (sc_options_t * opt,
                                int opt_char,
                                const char *opt_name,
                                const char **array_string,
                                const char *default_string,
                                double **double_array,
                                int initial_length,
                                const char *help_string)
{
    *double_array = NULL;
    sc_options_add_string (opt, opt_char, opt_name,
                           array_string, default_string, help_string);
}

void fclaw_options_convert_int_array (const char *array_string,
                                      int **int_array, int new_length)
{
    int i;
    const char *beginptr;
    char *endptr;

    new_length = SC_MAX (new_length, 0);
    *int_array = FCLAW_REALLOC (*int_array, int, new_length);

    beginptr = array_string;
    for (i = 0; i < new_length; ++i)
    {
        if (beginptr == NULL)
        {
            (*int_array)[i] = 0;
        }
        else
        {
            (*int_array)[i] = (int) strtol (beginptr, &endptr, 10);
            beginptr = endptr;
        }
    }
}

void fclaw_options_convert_double_array (const char *array_string,
                                        double **double_array, int new_length)
{
    int i;
    const char *beginptr;
    char *endptr;

    new_length = SC_MAX (new_length, 0);
    *double_array = FCLAW_REALLOC (*double_array, double, new_length);

    beginptr = array_string;
    for (i = 0; i < new_length; ++i)
    {
        if (beginptr == NULL)
        {
            (*double_array)[i] = 0;
        }
        else
        {
            (*double_array)[i] = (double) strtod (beginptr, &endptr);
            beginptr = endptr;
        }
    }
}

/* This is here to complement the array conversion routines above */
void fclaw_options_destroy_array(void* array)
{
    FCLAW_FREE (array);
}

/* ---------------------------------------------------------
   AMR access functions
   --------------------------------------------------------- */

amr_options_t* fclaw_forestclaw_get_options(fclaw_app_t* app)
{
    amr_options_t*  gparms;
    gparms = (amr_options_t*) fclaw_app_get_attribute(app,"fclaw-options",NULL);
    FCLAW_ASSERT(gparms != NULL);
    return gparms;
}
