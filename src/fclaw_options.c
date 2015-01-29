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

/** This is the internal state of an options structure for core variables. */

typedef struct fclaw_options_general
{
    amr_options_t *amropt;

    int print_help;        /**< Option variable to activate help message */
    int print_version;     /**< Option variable to print the version */
    int fclaw_verbosity;   /**< Option variable for ForestClaw verbosity */
    int lib_verbosity;     /**< Option variable for p4est, sc, and others */
    sc_keyvalue_t *kv_verbosity;      /**< Holds key-values for log levels */

    /* this is just for ForestClaw debugging, no need to adopt elsewhere */
    int is_registered;     /**< Internal variable to double-check the flow */
}
fclaw_options_general_t;

static void *
options_register_general (fclaw_app_t * a, void *package, sc_options_t * opt)
{
    sc_keyvalue_t *kv;
    fclaw_options_general_t *core = (fclaw_options_general_t *) package;

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (opt != NULL);

    /* allocated storage for this package's option values */
    FCLAW_ASSERT (core != NULL);
    FCLAW_ASSERT (!core->is_registered);

    /* this key-value pair understands the verbosity levels */
    kv = core->kv_verbosity = sc_keyvalue_new ();
    sc_keyvalue_set_int (kv, "default", FCLAW_VERBOSITY_DEFAULT);
    sc_keyvalue_set_int (kv, "debug", FCLAW_VERBOSITY_DEBUG);
    sc_keyvalue_set_int (kv, "info", FCLAW_VERBOSITY_INFO);
    sc_keyvalue_set_int (kv, "production", FCLAW_VERBOSITY_PRODUCTION);
    sc_keyvalue_set_int (kv, "essential", FCLAW_VERBOSITY_ESSENTIAL);
    sc_keyvalue_set_int (kv, "silent", FCLAW_VERBOSITY_SILENT);

    /* set the options for the core package */
    sc_options_add_switch (opt, 'h', "help", &core->print_help,
                           "Print usage information");
    sc_options_add_switch (opt, 'v', "version", &core->print_version,
                           "Print ForestClaw version");
    sc_options_add_keyvalue (opt, 'V', "verbosity", &core->fclaw_verbosity,
                             "default", kv, "Set ForestClaw verbosity");
    sc_options_add_keyvalue (opt, '\0', "lib-verbosity", &core->lib_verbosity,
                             "essential", kv, "Set verbosity for libraries");

    fclaw_options_add_general (opt, core->amropt);

    /* we do not need to work with the return value */
    core->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
options_postprocess_general (fclaw_app_t * a, void *package, void *registered)
{
    int fclaw_package_id;
    fclaw_options_general_t *core = (fclaw_options_general_t *) package;

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    /* errors from the key-value options would have showed up in parsing */

    /* postprocess this package */
    FCLAW_ASSERT (core != NULL);
    FCLAW_ASSERT (core->is_registered);

    /* go through this packages options */
    fclaw_package_id = fclaw_get_package_id();
    sc_package_set_verbosity (sc_package_id, core->lib_verbosity);
    sc_package_set_verbosity (p4est_package_id, core->lib_verbosity);
    sc_package_set_verbosity (fclaw_package_id, core->fclaw_verbosity);

    /* print help and/or version information and exit gracefully */
    if (core->print_help)
    {
        return FCLAW_EXIT_USAGE;
    }
    if (core->print_version)
    {
        fclaw_global_essentialf ("ForestClaw version %s\n",
                                 FCLAW_PACKAGE_VERSION);
        return FCLAW_EXIT_QUIET;
    }

    fclaw_options_postprocess(core->amropt);

    /* at this point there are no errors to report */
    return FCLAW_NOEXIT;
}

static void
options_destroy_general (fclaw_app_t * a, void *package, void *registered)
{
    fclaw_options_general_t *core = (fclaw_options_general_t *) package;

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    /* free this package */
    FCLAW_ASSERT (core != NULL);
    FCLAW_ASSERT (core->is_registered);
    FCLAW_ASSERT (core->kv_verbosity != NULL);
    sc_keyvalue_destroy (core->kv_verbosity);

    fclaw_options_destroy_arrays (core->amropt);
#if 0
    FCLAW_FREE(core->amropt);
#endif

    FCLAW_FREE (core);
}

static const fclaw_app_options_vtable_t options_vtable_general = {
    options_register_general,
    options_postprocess_general,
    NULL,
    options_destroy_general
};

void fclaw_app_options_register_general (fclaw_app_t * a, const char *configfile,
                                         amr_options_t* gparms)
{
    fclaw_options_general_t *core;

    FCLAW_ASSERT (a != NULL);

    /* allocate storage for core's option values */
    /* we will free it in the options_destroy callback */
    core = FCLAW_ALLOC_ZERO (fclaw_options_general_t, 1);
#if 0
    /* Allocate in calling program?  Or here? In any case, we allocate/
       deallocate option arrays here */
    *gparms = fclaw_options_new();  /* Requires amr_options_t** as input */
#endif
    core->amropt = gparms;

    /* sneaking the version string into the package pointer */
    /* when there are more parameters to pass, create a structure to pass */
    fclaw_app_options_register (a,NULL, configfile, &options_vtable_general,
                                core);
}


/* Use this with 'fclaw_options_destroy' */
amr_options_t* fclaw_options_new ()
{
    amr_options_t* amropt;
    amropt = FCLAW_ALLOC_ZERO (amr_options_t, 1);

    return amropt;
}

/* Use this with 'fclaw_options_new' */
void fclaw_options_destroy(amr_options_t* amropt)
{
    FCLAW_FREE (amropt);
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
    sc_options_add_int (opt, 0, "verbosity", &amropt->verbosity, 0,
                        "Verbosity mode [0]");
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

    sc_options_add_int (opt, 0, "refratio", &amropt->refratio, 2,
                        "Refinement ratio (fixed) [2]");

    sc_options_add_int (opt, 0, "minlevel", &amropt->minlevel, 0,
                        "Minimum refinement level [0]");

    sc_options_add_int (opt, 0, "maxlevel", &amropt->maxlevel, 0,
                        "Maximum refinement level[0]");

    sc_options_add_int (opt, 0, "regrid_interval", &amropt->regrid_interval,
                        1, "Regridding frequency [1]");


    sc_options_add_double (opt, 0, "ax", &amropt->ax, 0, "xlower " \
                           "(used only with manifold=0) [0]");
    sc_options_add_double (opt, 0, "bx", &amropt->bx, 1, "xupper " \
                           "(used only with manifold=0)[1]");
    sc_options_add_double (opt, 0, "ay", &amropt->ay, 0, "ylower " \
                           "(used only with manifold=0)[0]");
    sc_options_add_double (opt, 0, "by", &amropt->by, 1, "yupper " \
                           "(used only with manifold=0)[1]");

    sc_options_add_bool (opt, 0, "manifold", &amropt->manifold, 0,
                           "Solution is on manifold [F]");
    sc_options_add_bool (opt, 0, "use_fixed_dt", &amropt->use_fixed_dt, 0,
                           "Use fixed coarse grid time step [F]");
    sc_options_add_bool (opt, 0, "run_diagnostics",
                         &amropt->run_diagnostics,0,
                         "Run diagnostics [F]");
    sc_options_add_bool (opt, 0, "subcycle", &amropt->subcycle, 1,
                           "Use subcycling in time [F]");
    sc_options_add_bool (opt, 0, "noweightedp", &amropt->noweightedp, 0,
                           "No weighting when subcycling [F]");

    /* ---------------------- Usage information -------------------------- */
    sc_options_add_bool (opt, 0, "help", &amropt->help, 0,
                           "Print usage information (same as --usage) [F]");
    sc_options_add_bool (opt, 0, "usage", &amropt->help, 0,
                           "Print usage information (same as --help) [F]");

    sc_options_add_bool (opt, 0, "print_options", &amropt->print_options, 0,
                         "Print current option settings [F]");

    /* ---------------------- Control execution -------------------------- */
    sc_options_add_bool (opt, 0, "trapfpe", &amropt->trapfpe, 1,
                         "Trap floating point exceptions [1]");

    sc_options_add_bool (opt, 0, "mpi_debug", &amropt->mpi_debug, 0,
                        "Start MPI debug session (for attaching processes in gdb) [0]");


    sc_options_add_int (opt, 0, "mi", &amropt->mi, 1,
                        "Number of blocks in x direction [1]");

    sc_options_add_int (opt, 0, "mj", &amropt->mj, 1,
                        "Number of blocks in y direction  [1]");

    sc_options_add_bool (opt, 0, "periodic_x", &amropt->periodic_x, 0,
                        "Periodic in x direction [F]");

    sc_options_add_bool (opt, 0, "periodic_y", &amropt->periodic_y, 0,
                        "Periodic in y direction  [F]");

    /* --------------------------------------------------------------------
       Scale
       --------------------------------------------------------------------*/
    fclaw_options_add_double_array (opt,0, "scale",
                                    &amropt->scale_string, "1 1 1",
                                    &amropt->scale, 3,
                                    "Scale factor [1 1 1]");

    /* --------------------------------------------------------------------
       Shift
       --------------------------------------------------------------------*/
    fclaw_options_add_double_array (opt,0, "shift",
                                    &amropt->shift_string, "0 0 0",
                                    &amropt->shift, 3,
                                    "Shift array [0 0 0]");

    /* --------------------------------------------------------------------
       Rotate
       --------------------------------------------------------------------*/
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

void fclaw_options_destroy_arrays (amr_options_t * amropt)
{
    fclaw_options_destroy_array((void*) amropt->mthbc);
    fclaw_options_destroy_array((void*) amropt->scale);
    fclaw_options_destroy_array((void*) amropt->shift);
}


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


int fclaw_options_parse_command_line(sc_options_t * opt,
                                     int argc, char **argv)
{
    int retval;

    retval = sc_options_parse (sc_package_id, FCLAW_VERBOSITY_ESSENTIAL, opt, argc, argv);
    if (retval < 0)
    {
        fclaw_global_essentialf ("Command line option parsing failed.  " \
                                 "Use --help option to see valid option settings.\n");
        return retval;
    }
    else
    {
        fclaw_global_infof ("Reading command line.\n");
    }
    if (sc_is_root ())
    {
        retval = sc_options_save (sc_package_id, SC_LP_ERROR, opt,
                                  "fclaw_options.ini.used");
        if (retval < 0)
        {
            fclaw_global_essentialf("Cannot save fclaw_options.ini.used.  ");
            return retval;
        }
    }
    return retval;
}

void fclaw_options_print_summary(sc_options_t *opt)
{
    int fclaw_package_id;
    fclaw_package_id = fclaw_get_package_id();
    /* This is printed assuming vebosity level 'INFO' */
    sc_options_print_summary (fclaw_package_id, FCLAW_VERBOSITY_INFO,opt);
}

void fclaw_options_print_usage(sc_options_t *opt)
{
    int fclaw_package_id;
    fclaw_package_id = fclaw_get_package_id();
    /* This is printed assuming vebosity level 'INFO' */
    sc_options_print_usage (fclaw_package_id, FCLAW_VERBOSITY_INFO,opt,NULL);
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

void fclaw_options_destroy_array(void* array)
{
    FCLAW_FREE (array);
}

void fclaw_options_postprocess (amr_options_t * amropt)
{
    fclaw_options_convert_int_array (amropt->mthbc_string, &amropt->mthbc,
                                     fclaw2d_NumFaces);
    fclaw_options_convert_double_array (amropt->scale_string, &amropt->scale,3);
    fclaw_options_convert_double_array (amropt->shift_string, &amropt->shift,3);
}


/* -----------------------------------------------------------------
   Check input parms
   ----------------------------------------------------------------- */
int fclaw_options_check (sc_options_t * options, amr_options_t * gparms)
{
    /* Check for user help argument */
    if (gparms->help || gparms->print_options)
    {
        /* Both help and current settings can be printed */
        if (gparms->help)
        {
            /* This prints the help message for each options */
            fclaw_options_print_usage(options);
        }
        if (gparms->print_options)
        {
            /* This prints the current values of the parameters */
            fclaw_options_print_summary(options);
        }
        return -1;
    }

    /* Check outstyle. */
    if (gparms->outstyle == 1 && gparms->use_fixed_dt)
    {
        double dT_outer = gparms->tfinal / gparms->nout;
        double dT_inner = gparms->initial_dt;
        int nsteps = (int) floor (dT_outer / dT_inner + .5);
        if (fabs (nsteps * dT_inner - dT_outer) > 1e-8)
        {
            fclaw_global_essentialf( "For fixed dt, initial time step size must divide" \
                                     " tfinal/nout exactly.\n");
            return -1;
        }
    }

    /* Could also do basic sanity checks on mx,my,... */

    /* Everything is good */
    return 0;
}

#if 0
amr_options_t* fclaw_get_amr_options(fclaw_app_t *app)
{
    void* pkg = fclaw_app_get_options_package(app);
    fclaw_options_general_t* opt = (fclaw_options_general_t*) pkg;
    return opt->amropt;
}
#endif
