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

#include <fclaw2d_options.h>
#include <fclaw_options.h>
#include <fclaw2d_base.h>
#include <amr_options.h>
#include <forestclaw2d.h>

/* Proposed naming convention:
 * Parameter names in config file (= long option names) identical to the
 * C variable members of fclaw2d_options_t, except "-" in parameter name
 * corresponds to "_" in C variable.
 * For example the short option would be -F <filename> and the long option
 * --new-datafile=<Filename>.
 */

/* Use this with 'fclaw2d_options_destroy' */
amr_options_t* fclaw2d_options_new ()
{
    amr_options_t* amropt;
    amropt = FCLAW_ALLOC_ZERO (amr_options_t, 1);

    return amropt;
}

/* Use this with 'fclaw2d_options_new' */
void fclaw2d_options_destroy(amr_options_t* amropt)
{
    FCLAW_FREE (amropt);
}


void
fclaw2d_options_destroy_arrays (amr_options_t * amropt)
{
    FCLAW_FREE (amropt->mthbc);
}

#if 0
int fclaw2d_options_read_from_file(sc_options_t* opt, int log_priority)
{
    int retval;
    retval = sc_options_load (sc_package_id, SC_LP_ALWAYS, opt,
                              "fclaw2d_defaults.ini");
    if (retval < 0)
    {
        fclaw2d_global_log (log_priority,      \
                            "Cannot find (or cannot open) fclaw2d_defaults.ini.\n");
    }
    else
    {
        fclaw2d_global_log (log_priority,      \
                            "Reading file fclaw2d_defaults.ini.\n");
    }
    return retval;
}
#endif


void fclaw2d_options_register (sc_options_t * opt, amr_options_t* amropt)
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

    /* -----------------------------------------------------------------------
       Options will be read from this file, if a '-F' flag is used at the command
       line.  Use this file for local modifications that are not tracked by Git.
       ----------------------------------------------------------------------- */
    sc_options_add_inifile (opt, 'F', "inifile",
                            "[Options] File used to override one or more options " \
                            "in fclaw_options.ini [empty]");

}

void
fclaw2d_postprocess_parms (amr_options_t * amropt)
{
      fclaw_options_convert_int_array (amropt->mthbc_string, &amropt->mthbc,
                                         fclaw2d_NumFaces);
}


/* -----------------------------------------------------------------
   Check input parms
   ----------------------------------------------------------------- */
int
fclaw2d_checkparms (sc_options_t * options, amr_options_t * gparms, int lp)
{
    /* Check for user help argument */
    if (gparms->help)
    {
        sc_options_print_usage (sc_package_id, lp, options, NULL);
        return -1;
    }
    if (gparms->print_options)
    {
        fclaw_options_print_summary(options,lp);
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
            SC_GEN_LOG (sc_package_id, SC_LC_GLOBAL, lp,
                        "For fixed dt, initial time step size must divide"
                        " tfinal/nout exactly.\n");
            return -1;
        }
    }

    /* Could also do basic sanity checks on mx,my,... */

    /* Everything is good */
    return 0;
}
