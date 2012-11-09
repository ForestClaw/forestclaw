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

#include "amr_options.h"

#include "fclaw_defs.H"

/* Proposed naming convention:
 * Parameter names in config file (= long option names) identical to the
 * C variable members of amr_options_t, except "-" in parameter name
 * corresponds to "_" in C variable.
 * For example the short option would be -F <filename> and the long option
 * --new-datafile=<Filename>.
 */

void
amr_options_add_int_array (sc_options_t * opt,
                           int opt_char, const char *opt_name,
                           const char **array_string,
                           const char *default_string,
                           int **int_array, int initial_length,
                           const char *help_string)
{
    sc_options_add_string (opt, opt_char, opt_name,
                           array_string, default_string, help_string);
    amr_options_convert_int_array (*array_string, int_array, initial_length);
}

void
amr_options_convert_int_array (const char *array_string,
                               int **int_array, int new_length)
{
    int i;
    const char *beginptr;
    char *endptr;

    new_length = SC_MAX (new_length, 0);
    *int_array = SC_REALLOC (*int_array, int, new_length);

    beginptr = array_string;
    for (i = 0; i < new_length; ++i)
    {
        if (beginptr == NULL)
        {
            (*int_array)[i] = 0;
        }
        else
        {
            (*int_array)[i] = strtod (beginptr, &endptr);
            beginptr = endptr;
        }
    }
}

static void
amr_options_convert_arrays (amr_options_t * amropt)
{
    amr_options_convert_int_array (amropt->order_string, &amropt->order,
                                   SpaceDim);
    amr_options_convert_int_array (amropt->mthlim_string, &amropt->mthlim,
                                   amropt->mwaves);
    amr_options_convert_int_array (amropt->mthbc_string, &amropt->mthbc,
                                   CubeFaces);
}

amr_options_t *
amr_options_new (sc_options_t * opt)
{
    amr_options_t *amropt;

    amropt = SC_ALLOC_ZERO (amr_options_t, 1);

    sc_options_add_int (opt, 0, "mx", &amropt->mx, 0,
                        "Subdivision of each patch in x");

    sc_options_add_int (opt, 0, "my", &amropt->my, 0,
                        "Subdivision of each patch in y");

    sc_options_add_double (opt, 0, "initial_dt", &amropt->initial_dt, 0.0,
                           "Initial time step size");

    sc_options_add_switch (opt, 0, "use_fixed_dt", &amropt->use_fixed_dt,
                           "Use fixed coarse grid time step [F]");

    /*
       sc_options_add_double (opt, 0, "tfinal", &amropt->tfinal, 0.0,
       "Final time");
     */

    sc_options_add_int (opt, 0, "outstyle", &amropt->outstyle, 0.0,
                        "Output style (1,2,3)");

    sc_options_add_double (opt, 0, "max_cfl", &amropt->max_cfl, 1.0,
                           "Maximum CFL number [1.0]");

    sc_options_add_double (opt, 0, "desired_cfl", &amropt->desired_cfl, 0.9,
                           "Desired CFL number [0.9]");


    /*
       sc_options_add_int (opt, 0, "nout", &amropt->nout, 0,
       "Number of time steps");
     */

    /* Array of SpaceDim many values, with no defaults is set to all 0's */
    amr_options_add_int_array (opt, 0, "order", &amropt->order_string, NULL,
                               &amropt->order, SpaceDim,
                               "Normal and transverse orders");
    /* At this point amropt->order is allocated. Set defaults if desired. */

    sc_options_add_int (opt, 0, "verbosity", &amropt->verbosity, 0,
                        "Verbosity mode [0]");
    sc_options_add_int (opt, 0, "src_term", &amropt->src_term, 0,
                        "Source term option [0]");
    sc_options_add_int (opt, 0, "mcapa", &amropt->mcapa, -1,
                        "Location of capacity function in aux array [-1]");
    sc_options_add_int (opt, 0, "maux", &amropt->maux, 0,
                        "Number of auxiliary variables [0]");

    sc_options_add_int (opt, 0, "meqn", &amropt->meqn, 1,
                        "Number of equations");

    sc_options_add_int (opt, 0, "mwaves", &amropt->mwaves, 1,
                        "Number of waves");

    /* Array of mwaves many values */
    amr_options_add_int_array (opt, 0, "mthlim", &amropt->mthlim_string, NULL,
                               &amropt->mthlim, amropt->mwaves,
                               "Waves limiters (one for each wave)");
    /* At this point amropt->mthlim is allocated. Set defaults if desired. */

    sc_options_add_int (opt, 0, "mbc", &amropt->mbc, 2,
                        "Number of ghost cells [2]");

    /* Array of CubeFaces many values */
    amr_options_add_int_array (opt, 0, "mthbc", &amropt->mthbc_string, NULL,
                               &amropt->mthbc, CubeFaces,
                               "Physical boundary condition type");
    /* At this point amropt->mthbc is allocated. Set defaults if desired. */

    sc_options_add_int (opt, 0, "refratio", &amropt->refratio, 2,
                        "Refinement ratio (fixed) [2]");

    sc_options_add_int (opt, 0, "minlevel", &amropt->minlevel, 0,
                        "Minimum refinement level [0]");

    sc_options_add_int (opt, 0, "maxlevel", &amropt->maxlevel, 0,
                        "Maximum refinement level");

    sc_options_add_int (opt, 0, "regrid_interval", &amropt->regrid_interval,
                        0, "Regrid every ''regrid_interval'' steps");

#if 0                           /* bool is not allocated, disable for now */
    /* Does bool get allocated somewhere? */
    sc_options_add_string (opt, 0, "manifold", &bool, "F", "Manifold [F]");

    amropt->manifold = bool[0] == 'T' ? 1 : 0;

    sc_options_add_string (opt, 0, "mapped", &bool, "F", "Mapped grid [F]");
    amropt->mapped = bool[0] == 'T' ? 1 : 0;
#endif

#if 0
    /* Right now all switch options default to false, need to change that */
    sc_options_add_switch (opt, 0, "subcycle", &amropt->subcycle,
                           "T", "Use subcycling in time [T]");
#endif

    /* There is no default for this option.  A default will be read later. */
    sc_options_add_inifile (opt, 'F', "inifile",
                            "Read options from this file");

    /* It would be nice to have a default file that gets read,
       in case none is specified at the command line. */
    sc_options_load (sc_package_id, SC_LP_ALWAYS, opt, "fclaw_defaults.ini");
    /* Can check return value for -1 to see if file was not found */

    amr_options_convert_arrays (amropt);

    return amropt;
}

void
amr_options_parse (sc_options_t * opt, amr_options_t * amropt,
                   int argc, char **argv, int log_priority)
{
    int retval;

    retval = sc_options_parse (sc_package_id, SC_LP_ERROR, opt, argc, argv);
    if (retval < 0)
    {
        sc_options_print_usage (sc_package_id, log_priority, opt, NULL);
        sc_abort_collective ("Option parsing failed");
    }
    sc_options_print_summary (sc_package_id, log_priority, opt);
    if (sc_is_root ())
    {
        retval = sc_options_save (sc_package_id, SC_LP_ERROR, opt,
                                  "claw2ez.data.used");
        SC_CHECK_ABORT (!retval, "Option save failed");
    }

    amr_options_convert_arrays (amropt);
}

void
amr_options_destroy (amr_options_t * amropt)
{
    SC_FREE (amropt->order);
    SC_FREE (amropt->mthlim);
    SC_FREE (amropt->mthbc);
    SC_FREE (amropt);
}
