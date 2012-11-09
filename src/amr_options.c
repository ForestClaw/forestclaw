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

static void
amr_options_convert_arrays (amr_options_t * amropt)
{
    int i;
    int retval;
    int *o = amropt->order;

    if (amropt->order_string == NULL)
    {
        for (i = 0; i < SpaceDim; ++i)
        {
            o[i] = 0;           /* what would be a good value here? */
        }
    }
    else
    {
        retval = 0;
        switch (SpaceDim)
        {
        case 2:
            retval = sscanf (amropt->order_string, "%d %d", &o[0], &o[1]);
            break;
        case 3:
            retval = sscanf (amropt->order_string, "%d %d %d",
                             &o[0], &o[1], &o[2]);
            break;
        }
        if (retval != SpaceDim)
        {
            sc_abort_collective ("Option \"order\" needs DIM many values");
        }
    }
}

void
amr_options_register (sc_options_t * opt, amr_options_t * amropt)
{
    const char *bool;
    int vector_length;

    sc_options_add_int (opt, 0, "mx", &amropt->mx, 0,
                        "Subdivision of each patch in x");

    sc_options_add_int (opt, 0, "my", &amropt->my, 0,
                        "Subdivision of each patch in y");

    sc_options_add_double (opt, 0, "initial_dt", &amropt->initial_dt, 0.0,
                           "Initial time step size");

    sc_options_add_string (opt, 0, "use_fixed_dt", &bool, "F",
                           "Use fixed coarse grid time step [F]");
    amropt->use_fixed_dt = bool[0] == 'T' ? 1 : 0;

    /*
    sc_options_add_double (opt, 0, "tfinal", &amropt->tfinal, 0.0,
    "Final time");
    */

    sc_options_add_int(opt,0,"outstyle",&amropt->outstyle,0.0,"Output style (1,2,3)");

    sc_options_add_double (opt, 0, "max_cfl", &amropt->max_cfl, 1.0,
                           "Maximum CFL number [1.0]");

    sc_options_add_double (opt, 0, "desired_cfl", &amropt->desired_cfl, 0.9,
                           "Desired CFL number [0.9]");


    /*
    sc_options_add_int (opt, 0, "nout", &amropt->nout, 0,
    "Number of time steps");
    */

    /* Array of SpaceDim many values */
    sc_options_add_string (opt, 0, "order", &amropt->order_string, NULL,
                           "Normal and transverse orders");

    sc_options_add_int (opt, 0, "verbosity", &amropt->verbosity, 0,
                        "Verbosity mode [0]");
    sc_options_add_int (opt, 0, "src_term", &amropt->src_term, 0,
                        "Source term option [0]");
    sc_options_add_int (opt, 0, "mcapa", &amropt->mcapa, -1,
                        "Location of capacity function in aux array [-1]");
    sc_options_add_int (opt, 0, "maux", &amropt->maux, 0,
                        "Number of auxiliary variables [0]");

    sc_options_add_int (opt, 0, "meqn", &amropt->meqn,1,
                        "Number of equations");

    sc_options_add_int (opt, 0, "mwaves", &amropt->mwaves, 1,
                        "Number of waves");

    /*Length should be mwaves */
    vector_length = amropt->mwaves;
    /*
       sc_options_add_int_array (opt, 0, "mthlim", &amropt->mthlim,vector_length,
       "Waves limiters (one for each wave)");
     */

    sc_options_add_int (opt, 0, "mbc", &amropt->mbc, 2,
                        "Number of ghost cells [2]");

    /*
       vector_length = 2*SpaceDim;
       sc_options_add_int_array (opt, 0, "mthbc", &amropt->mthbc,vector_length,
       "Physical boundary condition type");
     */

    sc_options_add_int (opt, 0, "refratio", &amropt->refratio, 2,
                        "Refinement ratio (fixed) [2]");

    sc_options_add_int (opt, 0, "minlevel", &amropt->minlevel, 0,
                        "Minimum refinement level [0]");

    sc_options_add_int (opt, 0, "maxlevel", &amropt->maxlevel,0,
                        "Maximum refinement level");

    sc_options_add_int (opt, 0, "regrid_interval", &amropt->regrid_interval,0,
                        "Regrid every ''regrid_interval'' steps");

#if 0 /* bool is not allocated, disable for now */
    /* Does bool get allocated somewhere? */
    sc_options_add_string (opt, 0, "manifold", &bool, "F", "Manifold [F]");

    amropt->manifold = bool[0] == 'T' ? 1 : 0;

    sc_options_add_string (opt, 0, "mapped", &bool, "F", "Mapped grid [F]");
    amropt->mapped = bool[0] == 'T' ? 1 : 0;

    sc_options_add_string (opt, 0, "subcycle", &bool,
                           "T", "Use subcycling in time [T]");
    amropt->subcycle = bool[0] == 'T' ? 1 : 0;
#endif


    sc_options_add_inifile (opt, 'F', "fclaw_defaults.ini",
                            "Read options from this file");

    /* It would be nice to have a default file that gets read,
       in case none is specified at the command line. */

    amr_options_convert_arrays (amropt);
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
amr_options_delete (amr_options_t * amropt)
{
    // Need to delete this memory, but how?
    // delete [] amropts->mthlim;
}
