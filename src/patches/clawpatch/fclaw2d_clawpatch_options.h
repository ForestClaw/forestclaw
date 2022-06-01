/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun
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

#ifndef FCLAW2D_CLAWPATCH_OPTIONS_H
#define FCLAW2D_CLAWPATCH_OPTIONS_H

/** 
 * @file
 * Routines for handling clawpatch input options.
 */

#include <fclaw_base.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif



/* Criteria for tagging patches */
/** Refine based on value */
#define  FCLAW_REFINE_CRITERIA_VALUE          0
/** Refine based on derivative */
#define  FCLAW_REFINE_CRITERIA_DIFFERENCE     1
/** Refine based on difference between min and max */
#define  FCLAW_REFINE_CRITERIA_MINMAX         2
/** Refine based on gradient */
#define  FCLAW_REFINE_CRITERIA_GRADIENT       3
/** Refine based on user provided function */
#define  FCLAW_REFINE_CRITERIA_USER           4


struct fclaw2d_global;

/** Typedef for fclaw2d_clwaptch_opitons structure */
typedef struct fclaw2d_clawpatch_options fclaw2d_clawpatch_options_t;


/**
 * @brief Clawpatch options
 */
struct fclaw2d_clawpatch_options
{
    /* These are constant for all clawpatch's */
    int mx; /**< number of cells in the x direction */
    int my; /**< number of cells in the y direction */
    int maux; /**< number of aux equations */
    int mbc; /**< the number of ghost cells */

    int meqn; /**< number fields in solution */
    int rhs_fields; /**< number of rhs fields for elliptic problems */

    int refinement_criteria; /**< The refinement criteria */
    sc_keyvalue_t *kv_refinement_criteria; /**< The refinement criteria */

    /* Advanced options */
    int interp_stencil_width; /**< The width of the interpolation stencil */
    int ghost_patch_pack_aux; /**< True if aux equations should be packed */
    int save_aux;             /**< Save the aux array when retaking a time step */


    int is_registered; /**< true if options have been registered */

};

/**
 * @brief Register options from SC
 * 
 * @param app the app context
 * @param configfile the config file
 * @return fclaw2d_clawpatch_options_t* a newly allocated options struct
 */
fclaw2d_clawpatch_options_t *
fclaw2d_clawpatch_options_register(fclaw_app_t* app, const char* configfile);

/**
 * @brief Store the options in the global context
 * 
 * @param glob the global context
 * @param clawpatch_options the options
 */
void fclaw2d_clawpatch_options_store (struct fclaw2d_global *glob, 
                                      fclaw2d_clawpatch_options_t* clawpatch_options);

/**
 * @brief Get the options from the global context
 * 
 * @param glob the global context
 * @return fclaw2d_clawpatch_options_t* the options
 */
fclaw2d_clawpatch_options_t* fclaw2d_clawpatch_get_options(struct fclaw2d_global* glob);


#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* FCLAW2D_CLAWPATCH_OPTIONS_H */
