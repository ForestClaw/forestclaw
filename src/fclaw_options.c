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

#include <fclaw_options.h>
#include <fclaw_timer.h>
#include <fclaw_mpi.h>

/* Get whatever definitions exist already */
#ifdef FCLAW_HAVE_FENV_H
#ifndef __USE_GNU
#define __USE_GNU
#endif
#include <fenv.h>
#endif

/* Use as an alternate to GNU feenableexcept */
#ifndef FCLAW_HAVE_FEENABLEEXCEPT
#include <fp_exception_glibc_extension.h>
#endif

#ifdef FCLAW_HAVE_SIGNAL_H
#include <signal.h>
#endif

#ifdef FCLAW_HAVE_UNISTD_H
#include <unistd.h>    /* To get process ids */
#endif

static void* 
fclaw_register (fclaw_options_t* fclaw_opt, sc_options_t * opt)
{
    /* -------------------------- Time stepping control ------------------------------- */

    sc_options_add_double (opt, 0, "initial_dt", &fclaw_opt->initial_dt, 0.1,
                           "Initial time step size [0.1]");

    sc_options_add_double (opt, 0, "max_cfl", &fclaw_opt->max_cfl, 1,
                           "Maximum CFL allowed [1]");

    sc_options_add_double (opt, 0, "desired_cfl", &fclaw_opt->desired_cfl, 0.9,
                           "Maximum CFL allowed [0.9]");

    sc_options_add_bool (opt, 0, "reduce-cfl", &fclaw_opt->reduce_cfl, 1,
                           "Get maximum CFL over all processors [T]");

    sc_options_add_bool (opt, 0, "use_fixed_dt", &fclaw_opt->use_fixed_dt, 0,
                         "Use fixed coarse grid time step [F]");

    sc_options_add_int (opt, 0, "outstyle", &fclaw_opt->outstyle, 1,
                        "Output style (1,2,3) [1]");

    /* If outstyle == 1 */
    sc_options_add_double (opt, 0, "tfinal", &fclaw_opt->tfinal, 1.0,
                           "Final time [1.0]");

    sc_options_add_int (opt, 0, "nout", &fclaw_opt->nout, 10,
                        "Number of time steps, used with outstyle=3 [10]");

    /* Only needed if outstyle == 3 */
    sc_options_add_int (opt, 0, "nstep", &fclaw_opt->nstep, 1,
                        "Steps between output files, used with outstyle=3 [1]");

    sc_options_add_bool (opt, 0, "advance-one-step",
                         &fclaw_opt->advance_one_step, 0,
                         "Advance from t to t+dt in one global step (subcycle=F) [F]");

    sc_options_add_bool (opt, 0, "outstyle-uses-maxlevel",
                         &fclaw_opt->outstyle_uses_maxlevel, 0,
                         "Expand nout/nstep to global (subcycle=F) [F]");

    sc_options_add_bool (opt, 0, "subcycle", &fclaw_opt->subcycle, 1,
                         "Use subcycling in time [T]");

    sc_options_add_bool (opt, 0, "ghost-fill-uses-time-interp", &fclaw_opt->timeinterp2fillghost, 1,
                         "Use linear time interpolation when subcycling [T]");

    sc_options_add_bool (opt, 0, "weighted_partition", &fclaw_opt->weighted_partition, 1,
                         "Weight grids when partitioning [T]");

    /* ------------------------------ Conservation fix -------------------------------- */

    sc_options_add_bool (opt, 0, "time-sync", &fclaw_opt->time_sync, 0,
                         "Synchronize coarse/fine grids to include conservation fix [F]");

    sc_options_add_bool (opt, 0, "flux-correction", &fclaw_opt->flux_correction, 1,
                        "[user] Include flux correction [T]");

    sc_options_add_bool (opt, 0, "fluctuation-correction", 
                         &fclaw_opt->fluctuation_correction, 1,
                        "[user] Include fluctuation correction [T]");

    /* ------------------------------- Output options --------------------------------- */

    sc_options_add_bool (opt, 0, "output", &fclaw_opt->output, 0,
                            "Enable output [F]");


    /* -------------------------------------- Gauges  --------------------------------- */
    /* Gauge options */
    sc_options_add_bool (opt, 0, "output-gauges", &fclaw_opt->output_gauges, 0,
                            "Print gauge output [F]");

    sc_options_add_int(opt, 0, "gauge-buffer-length",
                       &fclaw_opt->gauge_buffer_length, 1,
                       "Number of lines of gauge output to buffer before printing [1]");

    /* -------------------------------- tikz output ----------------------------------- */
    sc_options_add_bool (opt, 0, "tikz-out", &fclaw_opt->tikz_out, 0,
                         "Enable tikz output for gridlines [F]");

    fclaw_options_add_double_array (opt, 0, "tikz-figsize", 
                                    &fclaw_opt->tikz_figsize_string,
                                    "8 6",&fclaw_opt->tikz_figsize,2,
                                    "Figure size used by tikz (inches) [8,6]");

    sc_options_add_string (opt, 0, "tikz-plot-prefix", 
                           &fclaw_opt->tikz_plot_prefix, 
                           "plot","Figure prefix for plotting [plot_]");    


    sc_options_add_string (opt, 0, "tikz-plot-suffix", 
                           &fclaw_opt->tikz_plot_suffix, 
                           "png","Figure suffix for plotting [png]");    

    sc_options_add_bool (opt, 0, "tikz-mesh-only", 
                         &fclaw_opt->tikz_mesh_only,0,
                         "Do not include .png file in tikz .tex output [0]");    

    /* Deprecated */
    sc_options_add_bool (opt, 0, "tikz-plot-fig", &fclaw_opt->tikz_plot_fig,1,
                           "Include figure (png,jpg,etc) in .tex file [1]");    



    /* --------------------------------- VTK output ----------------------------------- */
    /* This is a hack to control the VTK output while still in development.
     * The values are numbers which can be bitwise-or'd together.
     * 0 - no VTK output ever.
     * 1 - output for all stages of amrinit.
     * 2 - output whenever amr_output() is called.
     */

    /* Deprecated */
    sc_options_add_string (opt, 0, "prefix", &fclaw_opt->prefix, "fort",
                           "Output file prefix [fort]");

    /* Deprecated */
    sc_options_add_double (opt, 0, "vtkspace", &fclaw_opt->vtkspace, 0.,
                           "VTK visual spacing [F]");

    /* ----------------------------- Regridding options ------------------------------- */

    sc_options_add_bool (opt, 0, "init_ghostcell", 
                         &fclaw_opt->init_ghostcell, 0,
                        "Initialize ghost cells [F]");

    sc_options_add_int (opt, 0, "minlevel", &fclaw_opt->minlevel, 0,
                        "Minimum refinement level [0]");

    sc_options_add_int (opt, 0, "maxlevel", &fclaw_opt->maxlevel, 0,
                        "Maximum refinement level[0]");

    sc_options_add_int (opt, 0, "regrid_interval", 
                        &fclaw_opt->regrid_interval,
                        1, "Regridding frequency [1]");

    sc_options_add_int (opt, 0, "refratio", &fclaw_opt->refratio,
                        2, "Refinement ratio [2]");

    sc_options_add_bool (opt, 0, "smooth-refine",
                         &fclaw_opt->smooth_refine,
                         0, "Refinement smoothing[F]");

    sc_options_add_int (opt, 0, "smooth-level", 
                        &fclaw_opt->smooth_level,
                        0, "Lowest level for smooth-refine[0]");

    sc_options_add_int (opt, 0, "coarsen-delay", 
                        &fclaw_opt->coarsen_delay,
                        0, "Number skipped coarsenings[0]");

    sc_options_add_double (opt, 0, "refine_threshold", 
                           &fclaw_opt->refine_threshold,
                           0.5, "Refinement threshold [0.5]");

    sc_options_add_double (opt, 0, "coarsen_threshold", 
                           &fclaw_opt->coarsen_threshold,
                           0.1, "Coarsening threshold [0.1]");

    /* ---------------------------------- Diagnostics --------------------------------- */

    sc_options_add_bool (opt, 0, "run-user-diagnostics",
                         &fclaw_opt->run_user_diagnostics,0,
                         "Run user diagnostics [F]");

    sc_options_add_bool (opt, 0, "compute-error",
                         &fclaw_opt->compute_error,0,
                         "Compute error [F]");

    sc_options_add_bool (opt, 0, "conservation-check",
                         &fclaw_opt->conservation_check,0,
                         "Conservation check [F]");

    sc_options_add_bool (opt, 0, "report-timing",
                         &fclaw_opt->report_timing,1,
                         "Report timing results [T]");


    /* Set verbosity level for reporting timing */
    sc_keyvalue_t *kv = fclaw_opt->kv_timing_verbosity = sc_keyvalue_new ();
    sc_keyvalue_set_int (kv, "wall",      FCLAW_TIMER_PRIORITY_WALL);
    sc_keyvalue_set_int (kv, "summary",   FCLAW_TIMER_PRIORITY_SUMMARY);
    sc_keyvalue_set_int (kv, "exclusive", FCLAW_TIMER_PRIORITY_EXCLUSIVE);
    sc_keyvalue_set_int (kv, "counters",  FCLAW_TIMER_PRIORITY_COUNTERS);
    sc_keyvalue_set_int (kv, "details",   FCLAW_TIMER_PRIORITY_DETAILS);
    sc_keyvalue_set_int (kv, "extra",     FCLAW_TIMER_PRIORITY_EXTRA);
    sc_keyvalue_set_int (kv, "all",       FCLAW_TIMER_PRIORITY_EXTRA);
    sc_options_add_keyvalue (opt, 0, "report-timing-verbosity", 
                             &fclaw_opt->report_timing_verbosity,
                             "summary", kv, "Set verbosity for timing output [summary]");


    /* ---------------------------- Ghost packing options ----------------------------- */

    sc_options_add_bool (opt, 0, "ghost_patch_pack_area", 
                         &fclaw_opt->ghost_patch_pack_area,1,
                         "Pack area for parallel comm. of ghost patches [T]");

    sc_options_add_bool (opt, 0, "ghost_patch_pack_extra", 
                         &fclaw_opt->ghost_patch_pack_extra,
                        0, "Pack extra fields for parallel comm. of ghost patches [F]");

    sc_options_add_int (opt, 0, "ghost_patch_pack_numextrafields", 
                        &fclaw_opt->ghost_patch_pack_numextrafields,
                        0, "Number of extra fields to pack [0]");

    /* ---------------------------------- Debugging ----------------------------------- */

    sc_options_add_bool (opt, 0, "trapfpe", &fclaw_opt->trapfpe,1,
                         "Trap floating point exceptions [T]");

    sc_options_add_bool (opt, 0, "mpi_debug", &fclaw_opt->mpi_debug,0,
                        "Start MPI debug session (for attaching processes in gdb) [F]");


    /* ---------------------------- Domain geometry options --------------------------- */

    sc_options_add_double (opt, 0, "ax", &fclaw_opt->ax, 0, "xlower " \
                           "(used only with manifold=0) [0]");

    sc_options_add_double (opt, 0, "bx", &fclaw_opt->bx, 1, "xupper " \
                           "(used only with manifold=0) [1]");

    sc_options_add_double (opt, 0, "ay", &fclaw_opt->ay, 0, "ylower " \
                           "(used only with manifold=0) [0]");

    sc_options_add_double (opt, 0, "by", &fclaw_opt->by, 1, "yupper " \
                           "(used only with manifold=0) [1]");

    sc_options_add_double (opt, 0, "az", &fclaw_opt->az, 0, "zlower " \
                           "(used only with manifold=0) [0]");

    sc_options_add_double (opt, 0, "bz", &fclaw_opt->bz, 1, "zupper " \
                           "(used only with manifold=0) [1]");

    sc_options_add_bool (opt, 0, "manifold", &fclaw_opt->manifold, 0,
                         "Solution is on manifold [F]");

    sc_options_add_int (opt, 0, "mi", &fclaw_opt->mi, 1,
                        "Number of blocks in x direction [1]");

    sc_options_add_int (opt, 0, "mj", &fclaw_opt->mj, 1,
                        "Number of blocks in y direction  [1]");

    sc_options_add_bool (opt, 0, "periodic_x", &fclaw_opt->periodic_x, 0,
                        "Periodic in x direction [F]");

    sc_options_add_bool (opt, 0, "periodic_y", &fclaw_opt->periodic_y, 0,
                        "Periodic in y direction  [F]");

    fclaw_options_add_double_array (opt,0, "scale",
                                    &fclaw_opt->scale_string, "1 1 1",
                                    &fclaw_opt->scale, 3,
                                    "Scale factor [1 1 1]");

    fclaw_options_add_double_array (opt,0, "shift",
                                    &fclaw_opt->shift_string, "0 0 0",
                                    &fclaw_opt->shift, 3,
                                    "Shift array [0 0 0]");

    sc_options_add_double (opt, 0, "phi", &fclaw_opt->phi, 0,
                           "Rotation angle about x axis (degrees) [0]");

    sc_options_add_double (opt, 0, "theta", &fclaw_opt->theta, 0,
                           "Rotation angle about z axis (degrees) [0]");

    /* -----------------------------------------------------------------------
       Options will be read from this file, if a '-F' flag is used at the command
       line.  Use this file for local modifications that are not tracked by Git.
       ----------------------------------------------------------------------- */
    sc_options_add_inifile (opt, 'F', "inifile",
                            "File used to override one or more options " \
                            "in fclaw_options.ini [empty]");

    fclaw_opt->is_registered = 1;

    return NULL;
}


fclaw_exit_type_t 
fclaw_options_postprocess (fclaw_options_t * fclaw_opt)
{
    fclaw_options_convert_double_array (fclaw_opt->scale_string, &fclaw_opt->scale, 3);
    fclaw_options_convert_double_array (fclaw_opt->shift_string, &fclaw_opt->shift, 3);

    fclaw_options_convert_double_array (fclaw_opt->tikz_figsize_string, 
                                        &fclaw_opt->tikz_figsize, 2);

  return FCLAW_NOEXIT;
}

fclaw_exit_type_t
fclaw_options_check (fclaw_options_t * fclaw_opt)
{
    /* Check outstyle. */
    if (fclaw_opt->outstyle == 1 && fclaw_opt->use_fixed_dt)
    {
        double dT_outer = fclaw_opt->tfinal / fclaw_opt->nout;
        double dT_inner = fclaw_opt->initial_dt;
        int nsteps = (int) floor (dT_outer / dT_inner + .5);
        if (fabs (nsteps * dT_inner - dT_outer) > 1e-8)
        {
            fclaw_global_essentialf
              ("For fixed dt, initial time step size must divide"
               " tfinal/nout exactly.\n");
            return FCLAW_EXIT_ERROR;
        }
    }
    if (!fclaw_opt->reduce_cfl & !fclaw_opt->use_fixed_dt)
    {
        fclaw_global_essentialf("The CFL reduce can only be skipped if"
                                " use_fixed_dt = True\n");
        return FCLAW_EXIT_ERROR;
    }

    /* TODO: move these blocks to the beginning of forestclaw's control flow */
    if (fclaw_opt->mpi_debug)
    {
        fclaw_global_infof("Entering mpi_debug session");
        fclaw_mpi_debug ();
    }
#ifdef FCLAW_HAVE_FEENABLEEXCEPT
    if (fclaw_opt->trapfpe)
    {
        fclaw_global_infof("Enabling floating point traps\n");
        // feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
        feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    }
#endif

#if FCLAW2D_PATCHDIM == 3
    if (fclaw_opt->manifold)
    {
        fclaw_global_essentialf("Options : Manifold must be false for 3d patches\n");
        return FCLAW_EXIT_ERROR;        
    }
#endif

    return FCLAW_NOEXIT;
}

void
fclaw_options_destroy(fclaw_options_t* fclaw_opt)
{
    FCLAW_FREE (fclaw_opt->scale);
    FCLAW_FREE (fclaw_opt->shift);
    FCLAW_FREE (fclaw_opt->tikz_figsize);

    FCLAW_ASSERT (fclaw_opt->kv_timing_verbosity != NULL);
    sc_keyvalue_destroy (fclaw_opt->kv_timing_verbosity);
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
    FCLAW_ASSERT(fclaw_opt->is_registered);

    /* Convert strings to arrays */
    return fclaw_options_postprocess (fclaw_opt);
}


static fclaw_exit_type_t
options_check (fclaw_app_t * app, void *package, void *registered)
{
    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    fclaw_options_t *fclaw_opt = (fclaw_options_t *) package;
    FCLAW_ASSERT(fclaw_opt->is_registered);

    return fclaw_options_check (fclaw_opt);
}

static void
options_destroy (fclaw_app_t * a, void *package, void *registered)
{
    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    fclaw_options_t *fclaw_opt = (fclaw_options_t*) package;

    /* free this package */
    FCLAW_ASSERT (fclaw_opt->is_registered);

    /* Destroy option arrays created in post-process */
    fclaw_options_destroy (fclaw_opt);
    FCLAW_FREE(fclaw_opt);
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
fclaw_options_t* fclaw_options_register (fclaw_app_t * a, 
                                         const char *section,
                                         const char *configfile)
{
    fclaw_options_t* fclaw_opt;

    FCLAW_ASSERT (a != NULL);

    /* Basic options for print out help message, current options, etc */
    if(!fclaw_app_options_core_registered(a)){
        fclaw_app_options_register_core (a, configfile);
    }

    /* allocate storage for fclaw_options */
    /* we will free it in the options_destroy callback */
    fclaw_opt = FCLAW_ALLOC(fclaw_options_t,1);

    /* Could also pass in a section header (set to NULL for now) */
    fclaw_app_options_register (a,
                                section,
                                configfile,
                                &options_vtable,
                                fclaw_opt);
    
    return fclaw_opt;
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

