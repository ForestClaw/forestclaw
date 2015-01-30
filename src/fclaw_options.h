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

/** \file
 *
 * Routines for handling general ForestClaw input options.
 *
 */


#ifndef FCLAW_OPTIONS_H
#define FCLAW_OPTIONS_H

#include <fclaw_base.h>
#include <fclaw_options.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

/* Plan is to replace amr_options_t with fclaw_options_t */
typedef struct amr_options amr_options_t;


/** Create storage for option values specific to ForestClaw.
 * \return     New options container (values undefined).
 */
amr_options_t* fclaw_options_new ();

/** Clean up option storage.  Use fclaw_options_destroy_arrays to deallocate
 * storage for array options.
 * \param [in,out] amropt  Option storage will be deallocated.
 */
void fclaw_options_destroy(amr_options_t* amropt);


/** Free option arrays defined in ForestClaw container amr_options_t
 * \param [out] fclawopt        Option container for ForestClaw
 */

void fclaw_app_options_register_general (fclaw_app_t * a,
                                         const char *configfile,
                                         amr_options_t *gparms);


void fclaw_options_destroy_arrays (amr_options_t * fclawopt);

void fclaw_options_add_general (sc_options_t * opt, amr_options_t* amropt);

int  fclaw_options_read_from_file(sc_options_t* opt);

int  fclaw_options_parse_command_line (sc_options_t * opt, int argc,
                                       char **argv);

void fclaw_options_postprocess (amr_options_t * amropt);

int  fclaw_options_check (sc_options_t * options, amr_options_t * gparms);

void fclaw_options_print_summary(sc_options_t *opt);

void fclaw_options_print_usage(sc_options_t *opt);

/** Add a string option and prepare using it for an integer array.
 * \param [in,out] opt          Option container (see sc/sc_options.h).
 * \param [in] opt_char         Option character for command line (or 0).
 * \param [in] opt_name         Long option name for command line (or NULL).
 * \param [in,out] array_string Address that will point to the option string.
 * \param [in] default_string   Default string to be used or NULL.
 * \param [in,out] int_array    Pointer to an int array that gets resized
 *                              and populated with values from the string.
 * \param [in] initial_length   Initial length of int_array.
 * \param [in] help_string      Help message (used with --help).
 */
void fclaw_options_add_int_array (sc_options_t * opt,
                                int opt_char, const char *opt_name,
                                const char **array_string,
                                const char *default_string,
                                int **int_array, int initial_length,
                                const char *help_string);

void
fclaw_options_add_double_array (sc_options_t * opt,
                                int opt_char, const char *opt_name,
                                const char **array_string,
                                const char *default_string,
                                double **double_array, int initial_length,
                                const char *help_string);


/** Convert a string with multiple integers into an integer array.
 * \param [in] array_string     A string of space-separated integers.
 * \param [in,out] int_array    Pointer to an int array that gets resized
 *                              and populated with values from the string.
 *                              If string too short or NULL, set to 0.
 * \param [in] new_length       Length of int_array.
 */
void fclaw_options_convert_int_array (const char *array_string,
                                    int **int_array, int new_length);

void fclaw_options_convert_double_array (const char *array_string,
                                         double **double_array, int new_length);


void fclaw_options_destroy_array(void* array);

#if 0
amr_options_t* fclaw_get_amr_options(fclaw_app_t* app);
#endif



/* Plan is to replace amr_options_t with fclaw_options_t.
   Maybe use a macro as an intermediate step? */

struct amr_options
{
    /* Fixed grid size for each grid */
    int mx;      /**< Number of cells in x direction (fixed for all grids) */
    int my;      /**< Number of cells in y direction (fixed for all grids) */
    int mbc;     /**< Number of ghost cells in each grid */

    /* Time stepping */
    double initial_dt;  /**< Initial time step size */
    double tfinal;      /**< Final time */
    int outstyle;       /**< Outstyle */
    int nout;
    int nstep;
    int subcycle;               /**< Only relevant when subcycling. */
    int use_fixed_dt;
    double max_cfl;
    double desired_cfl;
    double *tout;

    /* Number of equations in the system of PDEs */
    int meqn;

    const char *mthbc_string;
    int *mthbc;

    /* Refinement parameters */
    int refratio;
    int minlevel;
    int maxlevel;
    int regrid_interval;

    /* Mapping functions */
    int manifold;
    int mi;
    int mj;
    int periodic_x;
    int periodic_y;

    const char *scale_string;
    double *scale;

    const char *shift_string;
    double *shift;

    double theta;
    double phi;

    double ax;   /**< Only for the single block, unmapped case */
    double bx;   /**< Only for the single block, unmapped case */
    double ay;   /**< Only for the single block, unmapped case */
    double by;   /**< Only for the single block, unmapped case */

    /* Diagnostics */
    int run_diagnostics;
    int conservation_check;
    int trapfpe;
    int mpi_debug;
#if 0
    int print_options;
#endif
    int help;

    /* Output and console IO */
    int verbosity;              /**< TODO: Do we have guidelines here? */
    int serialout;              /**< Allow for serial output.  WARNING:
                                     Will kill all parallel performance. */

    const char *prefix;         /**< This is prepended to output files */

    /* VTK output control */
    int vtkout;      /**< 0 for no output, 1 for output during amrinit,
                          2 for output when in amr_output.  Can be or'd. */
    double vtkspace; /**< between 0. and 1. to separate patches visually */
    int vtkwrite;    /**< 0 for MPI_File_write_all, 1 for MPI_File_write */

    int noweightedp;            /**< Don't use weighted partition. */

    int is_registered;     /**< Internal variable to double-check the flow */
};


#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW2D_OPTIONS_H */
