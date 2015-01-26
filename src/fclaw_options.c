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
#include <fclaw_base.h>


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


void fclaw_options_register (sc_options_t * opt, amr_options_t* amropt)
{

    sc_options_add_int (opt, 0, "mx", &amropt->mx, 8,
                        "[Options] Number of grid cells per patch in x [8]");

    sc_options_add_int (opt, 0, "my", &amropt->my, 8,
                        "[Options] Number of grid cells per patch in y [8]");

    sc_options_add_double (opt, 0, "initial_dt", &amropt->initial_dt, 0.1,
                           "[Options] Initial time step size [0.1]");

    sc_options_add_int (opt, 0, "outstyle", &amropt->outstyle, 1,
                        "[Options] Output style (1,2,3) [1]");

    /* If outstyle == 1 */
    sc_options_add_double (opt, 0, "tfinal", &amropt->tfinal, 1.0,
                           "[Options] Final time [1.0]");

    sc_options_add_int (opt, 0, "nout", &amropt->nout, 10,
                        "[Options] Number of time steps, used with outstyle=3 [10]");

    /* Only needed if outstyle == 3 */
    sc_options_add_int (opt, 0, "nstep", &amropt->nstep, 1,
                        "[Options] Steps between output files, used with outstyle=3 [1]");


    /* This is a hack to control the VTK output while still in development.
     * The values are numbers which can be bitwise-or'd together.
     * 0 - no VTK output ever.
     * 1 - output for all stages of amrinit.
     * 2 - output whenever amr_output() is called.
     */
    sc_options_add_int (opt, 0, "vtkout", &amropt->vtkout, 0,
                        "[Options] VTK output method [F]");
    sc_options_add_double (opt, 0, "vtkspace", &amropt->vtkspace, 0.,
                           "[Options] VTK visual spacing [F]");
    sc_options_add_int (opt, 0, "vtkwrite", &amropt->vtkwrite, 0,
                        "[Options] VTK write variant [F]");

    /* output options */
    sc_options_add_int (opt, 0, "verbosity", &amropt->verbosity, 0,
                        "[Options] Verbosity mode [0]");
    sc_options_add_bool (opt, 0, "serialout", &amropt->serialout, 1,
                            "[Options] Enable serial output [F]");
    sc_options_add_string (opt, 0, "prefix", &amropt->prefix, "fort",
                           "[Options] Output file prefix [fort]");

    /* more clawpack options */
    sc_options_add_double (opt, 0, "max_cfl", &amropt->max_cfl, 1,
                           "[Options] Maximum CFL allowed [1]");

    sc_options_add_double (opt, 0, "desired_cfl", &amropt->desired_cfl, 0.9,
                           "[Options] Maximum CFL allowed [0.9]");


    sc_options_add_int (opt, 0, "meqn", &amropt->meqn, 1,
                        "[Options] Number of equations [1]");

    sc_options_add_int (opt, 0, "mbc", &amropt->mbc, 2,
                        "[Options] Number of ghost cells [2]");

    /* Array of NumFaces many values */
    fclaw_options_add_int_array (opt, 0, "mthbc", &amropt->mthbc_string, "1 1 1 1",
        &amropt->mthbc, fclaw2d_NumFaces,
        "[Options] Physical boundary condition type [1 1 1 1]");


    /* At this point amropt->mthbc is allocated. Set defaults if desired. */

    sc_options_add_int (opt, 0, "refratio", &amropt->refratio, 2,
                        "[Options] Refinement ratio (fixed) [2]");

    sc_options_add_int (opt, 0, "minlevel", &amropt->minlevel, 0,
                        "[Options] Minimum refinement level [0]");

    sc_options_add_int (opt, 0, "maxlevel", &amropt->maxlevel, 0,
                        "[Options] Maximum refinement level[0]");

    sc_options_add_int (opt, 0, "regrid_interval", &amropt->regrid_interval,
                        1, "[Options] Regridding frequency [1]");


    sc_options_add_double (opt, 0, "ax", &amropt->ax, 0, "[Options] xlower " \
                           "(used only with manifold=0) [0]");
    sc_options_add_double (opt, 0, "bx", &amropt->bx, 1, "[Options] xupper " \
                           "(used only with manifold=0)[1]");
    sc_options_add_double (opt, 0, "ay", &amropt->ay, 0, "[Options] ylower " \
                           "(used only with manifold=0)[0]");
    sc_options_add_double (opt, 0, "by", &amropt->by, 1, "[Options] yupper " \
                           "(used only with manifold=0)[1]");

    /* -------------------------------------------------------------------
       CB: sc now has a new type of option, the bool.  While switch
       increments the value of the variable on each call, the bool
       can be initialized to either true or false and changed both ways.

       DC : Cool - thanks!
     */
    sc_options_add_bool (opt, 0, "manifold", &amropt->manifold, 0,
                           "[Options] Solution is on manifold [F]");
    sc_options_add_bool (opt, 0, "use_fixed_dt", &amropt->use_fixed_dt, 0,
                           "[Options] Use fixed coarse grid time step [F]");
    sc_options_add_bool (opt, 0, "run_diagnostics",
                         &amropt->run_diagnostics,0,
                         "[Options] Run diagnostics [F]");
    sc_options_add_bool (opt, 0, "subcycle", &amropt->subcycle, 1,
                           "[Options] Use subcycling in time [F]");
    sc_options_add_bool (opt, 0, "noweightedp", &amropt->noweightedp, 0,
                           "[Options] No weighting when subcycling [F]");

    /* ---------------------- Usage information -------------------------- */
    sc_options_add_bool (opt, 0, "help", &amropt->help, 0,
                           "[Options] Print usage information (same as --usage) [F]");
    sc_options_add_bool (opt, 0, "usage", &amropt->help, 0,
                           "[Options] Print usage information (same as --help) [F]");

    sc_options_add_bool (opt, 0, "print_options", &amropt->print_options, 0,
                         "[Options] Print current option settings [F]");

    /* ---------------------- Control execution -------------------------- */
    sc_options_add_bool (opt, 0, "trapfpe", &amropt->trapfpe, 1,
                         "[Options] Trap floating point exceptions [1]");

    sc_options_add_bool (opt, 0, "mpi_debug", &amropt->mpi_debug, 0,
                        "[Options] Start MPI debug session (for attaching processes in gdb) [0]");


    sc_options_add_int (opt, 0, "mi", &amropt->mi, 1,
                        "[Options] Number of blocks in x direction [1]");

    sc_options_add_int (opt, 0, "mj", &amropt->mj, 1,
                        "[Options] Number of blocks in y direction  [1]");

    sc_options_add_int (opt, 0, "periodic_x", &amropt->periodic_x, 0,
                        "[Options] Periodic in x direction [0]");

    sc_options_add_int (opt, 0, "periodic_y", &amropt->periodic_y, 0,
                        "[Options] Periodic in y direction  [0]");

    /* --------------------------------------------------------------------
       Scale
       --------------------------------------------------------------------*/
    fclaw_options_add_double_array (opt,0, "scale",
                                    &amropt->scale_string, "1 1 1",
                                    &amropt->scale, 3,
                                    "[Options] Scale factor [1 1 1]");

    /* --------------------------------------------------------------------
       Shift
       --------------------------------------------------------------------*/
    fclaw_options_add_double_array (opt,0, "shift",
                                    &amropt->shift_string, "0 0 0",
                                    &amropt->shift, 3,
                                    "[Options] Shift array [0 0 0]");

    /* --------------------------------------------------------------------
       Rotate
       --------------------------------------------------------------------*/
    sc_options_add_double (opt, 0, "phi", &amropt->phi, 0,
                           "[Options] Rotation angle about x axis (degrees) [0]");

    sc_options_add_double (opt, 0, "theta", &amropt->theta, 0,
                           "[Options] Rotation angle about z axis (degrees) [0]");


    /* -----------------------------------------------------------------------
       Options will be read from this file, if a '-F' flag is used at the command
       line.  Use this file for local modifications that are not tracked by Git.
       ----------------------------------------------------------------------- */
    sc_options_add_inifile (opt, 'F', "inifile",
                            "[Options] File used to override one or more options " \
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
    retval = sc_options_load (sc_package_id, SC_LP_ALWAYS, opt,
                              "fclaw_options.ini");
    if (retval < 0)
    {
        fclaw_global_essentialf( \
                            "Cannot find (or cannot open) fclaw_options.ini.\n");
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
    if (gparms->help)
    {
        fclaw_options_print_summary(options);
        return -1;
    }
    if (gparms->print_options)
    {
        fclaw_options_print_summary(options);
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

void fclaw_set_verbosity(sc_options_t* options,int *fclaw_verbosity, int p4est_verbosity)
{
    sc_keyvalue_t *kv_verbosity;
    kv_verbosity = sc_keyvalue_new ();

    sc_keyvalue_set_int (kv_verbosity, "default", FCLAW_VERBOSITY_DEFAULT);
    sc_keyvalue_set_int (kv_verbosity, "debug", FCLAW_VERBOSITY_DEBUG);
    sc_keyvalue_set_int (kv_verbosity, "info", FCLAW_VERBOSITY_INFO);
    sc_keyvalue_set_int (kv_verbosity, "production",FCLAW_VERBOSITY_PRODUCTION);
    sc_keyvalue_set_int (kv_verbosity, "essential",FCLAW_VERBOSITY_ESSENTIAL);
    sc_keyvalue_set_int (kv_verbosity, "silent", FCLAW_VERBOSITY_SILENT);

    sc_options_add_keyvalue (options, 'V', "fclaw-verbosity", fclaw_verbosity,
                             "default", kv_verbosity, "Set verbosity level");

#if 0
    /* This probably doesn't go here... */
    sc_package_set_verbosity (p4est_package_id,p4est_verbosity);
#endif
}
