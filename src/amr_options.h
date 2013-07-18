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

#ifndef AMR_OPTIONS_H
#define AMR_OPTIONS_H

#include <sc_options.h>
#include "fclaw2d_defs.H"

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

typedef struct amr_options amr_options_t;


typedef void *(*fclaw2d_user_parms_new_t) (sc_options_t * opt);

typedef void (*fclaw2d_user_parms_destroy_t) (amr_options_t * gparms);


struct amr_options
{
    /* Fixed grid size for each grid */
    int mx, my;

    /* Time stepping */
    double initial_dt;
    double tfinal;
    int outstyle;
    int nout;
    double *tout;
    int nstep;

    /* VTK output control */
    int vtkout;

    int verbosity;

    /* wave prop parameters */
    double max_cfl;
    double desired_cfl;

    /* Number of equations in the system of PDEs */
    int meqn;

    /* Boundary condition information */
    int mbc;

    const char *mthbc_string;
    int *mthbc;

    /* Refinement parameters */
    int refratio;
    int minlevel;
    int maxlevel;
    int regrid_interval;

    /* Boolean values (switches) */
    int manifold;
    int mapped;
    int subcycle;
    int check_conservation;
    int use_fixed_dt;

    /* xlower,xupper,ylower,yupper for block.
       This is probably only useful in the single block case  */
    double ax;
    double bx;
    double ay;
    double by;

};

/** Convert a string with multiple integers into an integer array.
 * \param [in] array_string     A string of space-separated integers.
 * \param [in,out] int_array    Pointer to an int array that gets resized
 *                              and populated with values from the string.
 *                              If string too short or NULL, set to 0.
 * \param [in] new_length       Length of int_array.
 */
void amr_options_convert_int_array (const char *array_string,
                                    int **int_array, int new_length);

/** Add a string option and prepare using it for an integer array.
 * \param [in,out] opt          Option container (see sc/sc_options.h).
 * \param [in] opt_char         Option character for command line (or 0).
 * \param [in] opt_name         Long option name for command line (or NULL).
 * \param [in,out] array_string Address that will point to the option string.
 * \param [in] default_string   Default string to be used or NULL.
 * \param [in,out] int_array    Pointer to an int array that gets resized
 *                              and populated with values from the string.
 * \param [in] initial_length   Initial length of int_array.
 */
void amr_options_add_int_array (sc_options_t * opt,
                                int opt_char, const char *opt_name,
                                const char **array_string,
                                const char *default_string,
                                int **int_array, int initial_length,
                                const char *help_string);

/** Create storage for option values specific to forestclaw.
 * \param [in,out] opt          Used for command line parsing.
 * \return                      Options with preset default values.
 */
amr_options_t *amr_options_new (sc_options_t * opt);

/* Parse options and populate values in registered amr_options structure.
 */
#if 0
void amr_options_parse (sc_options_t * opt, amr_options_t * amropt,
                        int argc, char **argv, int log_priority);
#endif
void amr_options_parse (sc_options_t * opt, int argc, char **argv,
                        int log_priority);


/** Clean up option storage.
 * \param [in,out]              Option storage will be deallocated.
 */
void amr_options_destroy (amr_options_t * amropt);


void amr_postprocess_parms (amr_options_t * amropt);
void amr_checkparms (amr_options_t * amropt);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !AMR_OPTIONS_H */
