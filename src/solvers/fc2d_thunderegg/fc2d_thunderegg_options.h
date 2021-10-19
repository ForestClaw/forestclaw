/*
Copyright (c) 2019 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#ifndef FCLAW2D_THUNDEREGG_OPTIONS_H
#define FCLAW2D_THUNDEREGG_OPTIONS_H

#include <fclaw_base.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

struct fclaw2d_global;

typedef struct fc2d_thunderegg_options fc2d_thunderegg_options_t;


typedef enum {
    STARPATCH = 0,  /* ThunderEgg operator */
    FIVEPOINT,      /* Laplacian (no beta) */
    VARPOISSON,      /* Variable Poisson operator */
    HEAT,      /* Variable Poisson operator */
    USER_OPERATOR
} fc2d_thunderegg_operator_types;

typedef enum {
    FFT = 0,    /* Must use starpatch or fivepoint */
    BICG,       /* Can be used with any operator */
    USER_SOLVER
} fc2d_thunderegg_solver_types;


struct fc2d_thunderegg_options
{
    /* Boundary conditions */
    int *boundary_conditions;
    const char *bc_cond_string;

    /* Output */
    int ascii_out;
    int vtk_out;

    /* bicgstab settings */
    int mg_prec;
    int max_it;
    double tol;

    /* thunderegg cyle settings */
    int pre_sweeps;
    int post_sweeps;
    int mid_sweeps;
    int coarse_sweeps;
    const char *cycle_type;

    int verbosity_level;

    /* bicgstab patch solver settings*/
    int patch_bcgs_max_it;
    double patch_bcgs_tol;

    /* Used by starpatch only */
    const char *patch_solver_type;

    int patch_operator;
    sc_keyvalue_t *kv_patch_operator;

    int patch_solver;
    sc_keyvalue_t *kv_patch_solver;


    int is_registered;
};

fclaw_exit_type_t fc2d_thunderegg_postprocess (fc2d_thunderegg_options_t *
                                               mg_opt);

fclaw_exit_type_t fc2d_thunderegg_check (fc2d_thunderegg_options_t * mg_opt);

void fc2d_thunderegg_reset (fc2d_thunderegg_options_t * mg_opt);

fc2d_thunderegg_options_t*  fc2d_thunderegg_options_register (fclaw_app_t * app,
                                                              const char *configfile);

void fc2d_thunderegg_package_register(fclaw_app_t* app,
                                     fc2d_thunderegg_options_t* mg_opt);

fc2d_thunderegg_options_t* fc2d_thunderegg_get_options(struct fclaw2d_global *glob);

void fc2d_thunderegg_options_store (struct fclaw2d_global* glob, 
                                   fc2d_thunderegg_options_t* mg_opt);

void fc2d_thunderegg_output(struct fclaw2d_global *glob, int iframe);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
