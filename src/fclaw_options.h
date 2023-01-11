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

/** \file
 *
 * Routines for handling general ForestClaw input options.
 *
 */

#ifndef FCLAW_OPTIONS_H
#define FCLAW_OPTIONS_H

#include <fclaw_base.h>

#ifdef __cplusplus
extern "C"
{
#endif

/* Plan is to replace fclaw_options_t with fclaw_options_t */
typedef struct fclaw_options fclaw_options_t;

/**
 * @brief Register options in SC
 * 
 * @param a the app context
 * @param section the section name
 * @param configfile the config file
 * @return fclaw2d_options_t* a newly allocated options struct
 */
fclaw_options_t* fclaw_options_register (fclaw_app_t * a,
                                         const char *section,
                                         const char *configfile);

/**
 * @brief Get the packing vtable for the options
 * 
 * @return const fclaw_userdata_vtable_t* the vtable
 */
const fclaw_userdata_vtable_t* fclaw_options_get_packing_vtable();

/* These can be called from external routines (in torthem, for example?) */
fclaw_exit_type_t 
fclaw_options_postprocess (fclaw_options_t * fclaw_opt);

fclaw_exit_type_t
fclaw_options_check (fclaw_options_t * fclaw_opt);

void fclaw_options_destroy(fclaw_options_t* fclaw_opt);



int fclaw_options_read_from_file(sc_options_t* opt);

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


/* Plan is to replace fclaw_options_t with fclaw_options_t.
   Maybe use a macro as an intermediate step? */

struct fclaw_options
{
    int dim;

    /* Run Options */
    const char* run_directory; /**< Directory where solver should be run */

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
    int reduce_cfl;   /* Do an all-reduce to get max. cfl */
    double *tout;

    /* Refinement parameters */
    int refratio;
    int minlevel;
    int maxlevel;
    int regrid_interval;
    int smooth_refine;
    int smooth_level;
    double refine_threshold;

    /* Conservation */
    int time_sync;

    /* Gauges */
    int output_gauges;
    int gauge_buffer_length;       

    int output_rays;

    /* Mapping functions */
    int manifold;
    int mi;
    int mj;
    int periodic_x;
    int periodic_y;

    /* Advanced options */
    int flux_correction;
    int fluctuation_correction;

    int coarsen_delay;
    double coarsen_threshold;

    /* Initialization of ghost cell */
    int init_ghostcell;

    /* Return after each time step  */
    int advance_one_step;

    /* nout, when used with outstyle option 3 refers to number of fine grid steps */
    int outstyle_uses_maxlevel;

    /* Do not use time interpolation.  This will be used for higher order schemes that 
       may use extra ghost cells, for example. */
    int timeinterp2fillghost;  


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
    // TODO 
    double az;
    double bz;

    /* Diagnostics */
    int run_user_diagnostics;
    int compute_error;
    int conservation_check;
    int trapfpe;
    int report_timing;
    int report_timing_verbosity;
    sc_keyvalue_t *kv_timing_verbosity;

    /* Parallel options */
    int mpi_debug;

    /* Parallel ghost patch comm. */
    int ghost_patch_pack_area;
    int ghost_patch_pack_extra;
    int ghost_patch_pack_numextrafields;

    /* Output and console IO */
    int verbosity;              /**< TODO: Do we have guidelines here? */

    int output;                    
    int tikz_out;      /* Boolean */

    const char *tikz_figsize_string;
    double *tikz_figsize;  /* In inches, e.g. [8,2] */

    int tikz_plot_fig;
    const char *tikz_plot_prefix;  /* For plotting */
    const char *tikz_plot_suffix;  /* For plotting */
    int tikz_mesh_only;

    const char *prefix;         /**< This is prepended to output files */

    double vtkspace; /**< between 0. and 1. to separate patches visually */

    int weighted_partition;            /**< Use weighted partition. */

    int is_registered;

    const char * logging_prefix; /**< prefix presented in logging ie. [prefix] */
};

#ifdef __cplusplus
}
#endif

#endif /* !FCLAW2D_OPTIONS_H */
