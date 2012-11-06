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
amr_options_register (sc_options_t * opt, amr_options_t * amropt)
{
    const char *bool;
    int vector_length;

    sc_options_add_int (opt, 0, "mx", &amropt->mx, 0,
                        "Subdivision of each patch in x [0]");

    sc_options_add_int (opt, 0, "my", &amropt->my, 0,
                        "Subdivision of each patch in y [0]");

    sc_options_add_double (opt, 0, "initial_dt", &amropt->initial_dt, 0.0,
                           "Initial time step size [0.0]");

    sc_options_add_double (opt, 0, "tfinal", &amropt->tfinal, 0.0,
                           "Final time [0.0]");

    sc_options_add_double (opt, 0, "max_cfl", &amropt->max_cfl, 1.0,
                           "Maximum CFL number [1.0]");

    sc_options_add_double (opt, 0, "desired_cfl", &amropt->desired_cfl, 0.9,
                           "Desired CFL number [0.9]");

    sc_options_add_int (opt, 0, "nout", &amropt->nout, 0,
                        "Number of time steps [0]");

    /*
      // Specify order of accuracy in each space dimension.
      vector_length = SpaceDim;
      sc_options_add_int_array(opt,0,"order",&amropt->order,vector_length,
                 "Normal and transverse order");
    */

    sc_options_add_int (opt, 0, "verbosity", &amropt->verbosity, 0,
                        "Verbosity mode [0]");
    sc_options_add_int (opt, 0, "src_term",&amropt->src_term,0,
                        "Source term option [0]");
    sc_options_add_int (opt, 0, "mcapa", &amropt->mcapa,-1,
                        "Location of capacity function in aux array [-1]");
    sc_options_add_int (opt, 0, "maux", &amropt->maux,0,
                        "Number of auxiliary variables [0]");

    sc_options_add_int (opt, 0, "meqn", &amropt->meqn,1,
                        "Number of equations [1]");
    sc_options_add_int (opt, 0, "mwaves", &amropt->mwaves, 1,
                        "Number of waves [1]");

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

    sc_options_add_int (opt, 0, "refratio", &amropt->refratio,2,
                        "Refinement ratio (fixed) [2]");

    sc_options_add_int (opt, 0, "minlevel", &amropt->minlevel,0,
                        "Minimum refinement level [0]");

    sc_options_add_int (opt, 0, "maxlevel", &amropt->maxlevel,0,
                        "Maximum refinement level [0]");

    /* Does bool get allocated somewhere? */
    sc_options_add_string (opt, 0, "manifold", &bool,
                           "F", "Manifold [F]");

    amropt->manifold = bool[0] == 'T' ? 1 : 0;

    sc_options_add_string (opt, 0, "mapped", &bool, "F", "Mapped grid [F]");
    amropt->mapped = bool[0] == 'T' ? 1 : 0;

    sc_options_add_string (opt, 0, "subcycle", &bool,
                           "T", "Use subcycling in time [T]");
    amropt->subcycle = bool[0] == 'T' ? 1 : 0;

    /* Is there another file that gets read if this file is not specified at the
       command line ? */
    sc_options_add_inifile (opt, 'F', "fclaw_defaults.ini",
                            "Read options from this file");
}

void
amr_options_parse (sc_options_t * opt, int argc, char **argv,
                   int log_priority)
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
}
