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

#include <fclaw2d_forestclaw.h>
#include <fclaw_mpi.h>

/* amr_options_t will eventually be replaced with fclaw_options_t. */
typedef struct fclaw_options
{
    amr_options_t *gparms;
} fclaw_options_t;


static void* 
fclaw_register (fclaw_options_t* fclaw_opt, sc_options_t * opt)
{
    amr_options_t *gparms = fclaw_opt->gparms;

    /* ---------------------- Time stepping control -------------------------- */

    sc_options_add_double (opt, 0, "initial_dt", &gparms->initial_dt, 0.1,
                           "Initial time step size [0.1]");

    sc_options_add_double (opt, 0, "max_cfl", &gparms->max_cfl, 1,
                           "Maximum CFL allowed [1]");

    sc_options_add_double (opt, 0, "desired_cfl", &gparms->desired_cfl, 0.9,
                           "Maximum CFL allowed [0.9]");

    sc_options_add_bool (opt, 0, "use_fixed_dt", &gparms->use_fixed_dt, 0,
                         "Use fixed coarse grid time step [F]");

    sc_options_add_int (opt, 0, "outstyle", &gparms->outstyle, 1,
                        "Output style (1,2,3) [1]");

    /* If outstyle == 1 */
    sc_options_add_double (opt, 0, "tfinal", &gparms->tfinal, 1.0,
                           "Final time [1.0]");

    sc_options_add_int (opt, 0, "nout", &gparms->nout, 10,
                        "Number of time steps, used with outstyle=3 [10]");

    /* Only needed if outstyle == 3 */
    sc_options_add_int (opt, 0, "nstep", &gparms->nstep, 1,
                        "Steps between output files, used with outstyle=3 [1]");

    sc_options_add_bool (opt, 0, "advance-one-step",
                         &gparms->advance_one_step, 1,
                         "Advance from t to t+dt in one global step (subcycle=F) [F]");

    sc_options_add_bool (opt, 0, "outstyle-uses-maxlevel",
                         &gparms->outstyle_uses_maxlevel, 1,
                         "Expand nout/nstep to global (subcycle=F) [F]");

    sc_options_add_bool (opt, 0, "subcycle", &gparms->subcycle, 1,
                         "Use subcycling in time [F]");

    sc_options_add_bool (opt, 0, "weighted_partition", &gparms->weighted_partition, 0,
                         "Weight grids when partitioning [F]");

    /* ---------------------- Output options -------------------------- */

    sc_options_add_bool (opt, 0, "output", &gparms->output, 0,
                            "Enable output [F]");


    /* tikz output */
    sc_options_add_bool (opt, 0, "tikz-out", &gparms->tikz_out, 0,
                         "Enable tikz output for gridlines [F]");

    fclaw_options_add_double_array (opt, 0, "tikz-figsize", &gparms->tikz_figsize_string,
                                    "8 6",&gparms->tikz_figsize,2,
                                    "Figure size used by tikz (inches) [8,6]");

    sc_options_add_string (opt, 0, "tikz-plot-prefix", &gparms->tikz_plot_prefix, "plot",
                           "Figure prefix for plotting [plot_]");    


    sc_options_add_string (opt, 0, "tikz-plot-suffix", &gparms->tikz_plot_suffix, "png",
                           "Figure suffix for plotting [png]");    

    /* This is a hack to control the VTK output while still in development.
     * The values are numbers which can be bitwise-or'd together.
     * 0 - no VTK output ever.
     * 1 - output for all stages of amrinit.
     * 2 - output whenever amr_output() is called.
     */

    sc_options_add_string (opt, 0, "prefix", &gparms->prefix, "fort",
                           "Output file prefix [fort]");

    sc_options_add_int (opt, 0, "vtkout_debug", &gparms->vtkout_debug, 0,
                        "VTK debug output method [F]");
 
    sc_options_add_double (opt, 0, "vtkspace", &gparms->vtkspace, 0.,
                           "VTK visual spacing [F]");

    sc_options_add_int (opt, 0, "vtkwrite", &gparms->vtkwrite, 0,
                        "VTK write variant [F]");

    /* ---------------------- Regridding options -------------------------- */

    sc_options_add_bool (opt, 0, "init_ghostcell", &gparms->init_ghostcell, 1,
                        "Initialize ghost cells [T]");

    sc_options_add_int (opt, 0, "minlevel", &gparms->minlevel, 0,
                        "Minimum refinement level [0]");

    sc_options_add_int (opt, 0, "maxlevel", &gparms->maxlevel, 0,
                        "Maximum refinement level[0]");

    sc_options_add_int (opt, 0, "regrid_interval", &gparms->regrid_interval,
                        1, "Regridding frequency [1]");

    sc_options_add_int (opt, 0, "refratio", &gparms->refratio,
                        2, "Refinement ratio [2]");

    sc_options_add_bool (opt, 0, "smooth-refine", &gparms->smooth_refine,
                         0, "Refinement smoothing[F]");

    sc_options_add_int (opt, 0, "smooth-level", &gparms->smooth_refine_level,
                        0, "Lowest level for smooth-refine[0]");

    sc_options_add_int (opt, 0, "coarsen-delay", &gparms->coarsen_delay,
                        0, "Number skipped coarsenings[0]");

    sc_options_add_double (opt, 0, "refine_threshold", &gparms->refine_threshold,
                           0.5, "Refinement threshold [0.5]");

    sc_options_add_double (opt, 0, "coarsen_threshold", &gparms->coarsen_threshold,
                           0.1, "Coarsening threshold [0.1]");

     /* ---------------------- Diagnostics -------------------------- */

   sc_options_add_bool (opt, 0, "run-user-diagnostics",
                         &gparms->run_user_diagnostics,0,
                         "Run user diagnostics [F]");

    sc_options_add_bool (opt, 0, "compute-error",
                         &gparms->compute_error,0,
                         "Compute error [F]");

    sc_options_add_bool (opt, 0, "conservation-check",
                         &gparms->conservation_check,0,
                         "Conservation check [F]");

    sc_options_add_bool (opt, 0, "report-timing",
                         &gparms->report_timing,1,
                         "Report timing results [T]");

    /* ---------------------- Ghost packing options -------------------------- */

    sc_options_add_bool (opt, 0, "ghost_patch_pack_area", &gparms->ghost_patch_pack_area,1,
                         "Pack area for parallel comm. of ghost patches [T]");

    sc_options_add_bool (opt, 0, "ghost_patch_pack_extra", &gparms->ghost_patch_pack_extra,
                        0, "Pack extra fields for parallel comm. of ghost patches [F]");

    sc_options_add_int (opt, 0, "ghost_patch_pack_numextrafields", 
                        &gparms->ghost_patch_pack_numextrafields,
                        0, "Number of extra fields to pack [0]");

    /* ---------------------- Debugging -------------------------- */

    sc_options_add_bool (opt, 0, "trapfpe", &gparms->trapfpe,1,
                         "Trap floating point exceptions [T]");

    sc_options_add_bool (opt, 0, "mpi_debug", &gparms->mpi_debug,0,
                        "Start MPI debug session (for attaching processes in gdb) [F]");


    /* ---------------------- Domain geometry options -------------------------- */

    sc_options_add_double (opt, 0, "ax", &gparms->ax, 0, "xlower " \
                           "(used only with manifold=0) [0]");

    sc_options_add_double (opt, 0, "bx", &gparms->bx, 1, "xupper " \
                           "(used only with manifold=0) [1]");

    sc_options_add_double (opt, 0, "ay", &gparms->ay, 0, "ylower " \
                           "(used only with manifold=0) [0]");

    sc_options_add_double (opt, 0, "by", &gparms->by, 1, "yupper " \
                           "(used only with manifold=0) [1]");

    sc_options_add_bool (opt, 0, "manifold", &gparms->manifold, 0,
                         "Solution is on manifold [F]");

    sc_options_add_int (opt, 0, "mi", &gparms->mi, 1,
                        "Number of blocks in x direction [1]");

    sc_options_add_int (opt, 0, "mj", &gparms->mj, 1,
                        "Number of blocks in y direction  [1]");

    sc_options_add_bool (opt, 0, "periodic_x", &gparms->periodic_x, 0,
                        "Periodic in x direction [F]");

    sc_options_add_bool (opt, 0, "periodic_y", &gparms->periodic_y, 0,
                        "Periodic in y direction  [F]");

    fclaw_options_add_double_array (opt,0, "scale",
                                    &gparms->scale_string, "1 1 1",
                                    &gparms->scale, 3,
                                    "Scale factor [1 1 1]");

    fclaw_options_add_double_array (opt,0, "shift",
                                    &gparms->shift_string, "0 0 0",
                                    &gparms->shift, 3,
                                    "Shift array [0 0 0]");

    sc_options_add_double (opt, 0, "phi", &gparms->phi, 0,
                           "Rotation angle about x axis (degrees) [0]");

    sc_options_add_double (opt, 0, "theta", &gparms->theta, 0,
                           "Rotation angle about z axis (degrees) [0]");


    /* ---------------------- advanced options -------------------------- */
    sc_options_add_int (opt, 0, "interp_stencil_width",
                        &gparms->interp_stencil_width,
                        3, "Interpolation stencil width [3]");


    /* -----------------------------------------------------------------------
       Options will be read from this file, if a '-F' flag is used at the command
       line.  Use this file for local modifications that are not tracked by Git.
       ----------------------------------------------------------------------- */
    sc_options_add_inifile (opt, 'F', "inifile",
                            "File used to override one or more options " \
                            "in fclaw_options.ini [empty]");

    gparms->is_registered = 1;

    return NULL;
}


static fclaw_exit_type_t 
fclaw_postprocess (fclaw_options_t * fclaw_opt)
{
    amr_options_t *gparms = fclaw_opt->gparms;

    fclaw_options_convert_double_array (gparms->scale_string, &gparms->scale, 3);
    fclaw_options_convert_double_array (gparms->shift_string, &gparms->shift, 3);

    fclaw_options_convert_double_array (gparms->tikz_figsize_string, &gparms->tikz_figsize, 2);

  return FCLAW_NOEXIT;
}

static fclaw_exit_type_t
fclaw_check (fclaw_options_t * fclaw_opt, 
             fclaw2d_clawpatch_options_t *clawpatch_opt)
{

    amr_options_t *gparms = fclaw_opt->gparms;

    if (gparms->interp_stencil_width/2 > clawpatch_opt->mbc)
    {
        fclaw_global_essentialf("Interpolation width is too large for number of " \
                                "ghost cells (mbc) specifed.  We should have " \
                                "(width)/2 <= mbc");
    }

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

void
fclaw_destroy(fclaw_options_t* fclaw_opt)
{
    amr_options_t *gparms = fclaw_opt->gparms;
    FCLAW_FREE (gparms->scale);
    FCLAW_FREE (gparms->shift);
    FCLAW_FREE (gparms->tikz_figsize);
    FCLAW_FREE (gparms);
}


/* ------------------------------------------------------------------------
  Generic functions - these call the functions above
  ------------------------------------------------------------------------ */
static void *
options_register(fclaw_app_t * a, void *package, sc_options_t * opt)
{
    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (opt != NULL);

    fclaw_options_t *fclaw_opt = (fclaw_options_t *) package;

    /* we do not need to work with the return value */
    return fclaw_register(fclaw_opt,opt);
}

static fclaw_exit_type_t
options_postprocess (fclaw_app_t * a, void *package, void *registered)
{
    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    /* errors from the key-value options would have showed up in parsing */
    fclaw_options_t *fclaw_opt = (fclaw_options_t *) package;

    /* post-process this package */
    FCLAW_ASSERT(fclaw_opt->gparms->is_registered);

    /* Convert strings to arrays */
    return fclaw_postprocess (fclaw_opt);
}


static fclaw_exit_type_t
options_check (fclaw_app_t * app, void *package, void *registered)
{
    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    fclaw_options_t *fclaw_opt = (fclaw_options_t *) package;
    FCLAW_ASSERT(fclaw_opt->gparms->is_registered);

    fclaw2d_clawpatch_options_t *clawpatch_opt;
    clawpatch_opt = fclaw_app_get_attribute(app,"clawpatch",NULL);
    FCLAW_ASSERT(clawpatch_opt->is_registered);

    return fclaw_check (fclaw_opt,clawpatch_opt);
}

static void
options_destroy (fclaw_app_t * a, void *package, void *registered)
{
    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    fclaw_options_t *fclaw_opt = (fclaw_options_t*) package;

    /* free this package */
    FCLAW_ASSERT (fclaw_opt->gparms->is_registered);

    /* Destroy option arrays created in post-process */
    fclaw_destroy (fclaw_opt);
    FCLAW_FREE(fclaw_opt->gparms);
}

static const fclaw_app_options_vtable_t options_vtable = {
    options_register,
    options_postprocess,
    options_check,
    options_destroy
};


/* ---------------------------------------------------------
   Public interface to ForestClaw options
   --------------------------------------------------------- */
amr_options_t* fclaw_options_register (fclaw_app_t * a, 
                                       const char *configfile)
{
    fclaw_options_t* fclaw_opt;

    FCLAW_ASSERT (a != NULL);

    /* Basic options for print out help message, current options, etc */
    fclaw_app_options_register_core (a, configfile);

    /* allocate storage for fclaw_options */
    /* we will free it in the options_destroy callback */
    fclaw_opt = FCLAW_ALLOC(fclaw_options_t,1);
    fclaw_opt->gparms = FCLAW_ALLOC(amr_options_t, 1);

    /* Could also pass in a section header (set to NULL for now) */
    fclaw_app_options_register (a,NULL,
                                configfile,
                                &options_vtable,
                                fclaw_opt);
    
    /* this is to retrieve the option key-value pair */
    fclaw_app_set_attribute(a,"Options",fclaw_opt->gparms);

    return fclaw_opt->gparms;
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
        fclaw_global_essentialf("Problem reading fclaw_options.ini.\n");
    }
    else
    {
        fclaw_global_infof ("Reading file fclaw_options.ini.\n");
    }
    return retval;
}

#if 0
void fclaw2d_options_store (fclaw2d_global_t *glob, amr_options_t* gparms)
{
    int id;

    FCLAW_ASSERT(s_forestclaw_package_id == -1);
    id = fclaw_package_container_add_pkg(glob,gparms);
    s_forestclaw_package_id = id;
}

amr_options_t* fclaw2d_get_options(fclaw2d_global_t* glob)
{
    amr_options_t *gp = (amr_options_t*) 
            fclaw_package_get_options(glob, s_forestclaw_package_id);
    return gp;
}
#endif

/* ------------------------------------------------------------------------
  Some utility functions
  ------------------------------------------------------------------------ */

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

