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

#ifndef FCLAW_CLAWPATCH_OPTIONS_H
#define FCLAW_CLAWPATCH_OPTIONS_H

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


struct fclaw_global;

/** Typedef for fclaw_clawpatch_options structure */
typedef struct fclaw_clawpatch_options fclaw_clawpatch_options_t;


/**
 * @brief Clawpatch options
 */
struct fclaw_clawpatch_options
{
    /* These are constant for all clawpatch's */
    int patch_dim; /**< dimension of clawpatch */
    int mx; /**< number of cells in the x direction */
    int my; /**< number of cells in the y direction */
    int mz; /**< number of cells in the z direction */
    int maux; /**< number of aux equations */
    int mbc; /**< the number of ghost cells */

    int meqn; /**< number fields in solution */
    int rhs_fields; /**< number of rhs fields for elliptic problems */

    int refinement_criteria; /**< The refinement criteria */
    sc_keyvalue_t *kv_refinement_criteria; /**< The refinement criteria */
    int threshold_variable;

    /* Advanced options */
    int interp_stencil_width; /**< The width of the interpolation stencil */
    int ghost_patch_pack_aux; /**< True if aux equations should be packed */
    int save_aux;             /**< Save the aux array when retaking a time step */


    int vtk_patch_threshold; /**< The buffer threshold for vtk output */
    int hdf5_patch_threshold; /**< The buffer threshold for hdf5 output */
    int hdf5_compression_level; /**< The compression level for hdf5 output */

    int is_registered; /**< true if options have been registered */

};

/**
 * @brief Allocate a new fclaw_clawpatch_options_t structure with the given dimension.
 *
 * This shouldn't have be called by the user.
 * Used mainly for testing.
 * 
 * @param dim The dimension of the clawpatch.
 * @return fclaw_clawpatch_options_t* A newly allocated fclaw_clawpatch_options_t structure.
 */
fclaw_clawpatch_options_t* fclaw_clawpatch_options_new (int dim);

/**
 * @brief Deallocate the given fclaw_clawpatch_options_t structure.
 * 
 * This shouldn't have be called by the user.
 * Used mainly for testing.
 *
 * @param clawpatch_opt The fclaw_clawpatch_options_t structure to deallocate.
 */
void fclaw_clawpatch_options_destroy (fclaw_clawpatch_options_t *clawpatch_opt);

/**
 * @brief Register 2d options in SC
 * 
 * @param app the app context
 * @param name the name of the options group
 * @param configfile the config file
 * @return fclaw_clawpatch_options_t* a newly allocated options struct
 */
fclaw_clawpatch_options_t *
fclaw_clawpatch_2d_options_register(fclaw_app_t* app, const char* name, const char* configfile);

/**
 * @brief Register 3d options in SC
 * 
 * @param app the app context
 * @param name the name of the options group
 * @param configfile the config file
 * @return fclaw_clawpatch_options_t* a newly allocated options struct
 */
fclaw_clawpatch_options_t *
fclaw_clawpatch_3d_options_register(fclaw_app_t* app, const char* name, const char* configfile);

/**
 * @brief Register dimension independent options. 
 *        This will expose a dim option in the clawpatch options.
 * 
 * @param app the app context
 * @param name the name of the options group
 * @param configfile the config file
 * @return fclaw_clawpatch_options_t* a newly allocated options struct
 */
fclaw_clawpatch_options_t *
fclaw_clawpatch_dim_ind_options_register(fclaw_app_t* app, const char* name, const char* configfile);



/**
 * @brief Store the options in the global context
 * 
 * @param glob the global context
 * @param clawpatch_options the options
 */
void fclaw_clawpatch_options_store (struct fclaw_global *glob, 
                                    fclaw_clawpatch_options_t* clawpatch_options);

/**
 * @brief Get the options from the global context
 * 
 * @param glob the global context
 * @return fclaw_clawpatch_options_t* the options
 */
fclaw_clawpatch_options_t* fclaw_clawpatch_get_options(struct fclaw_global* glob);

/**
 * @brief Get the packing vtable for fclaw_clawpatch_options_t
 * 
 * @return const fclaw_packing_vtable_t* the vtable
 */
const fclaw_packing_vtable_t* fclaw_clawpatch_options_get_packing_vtable();


#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* FCLAW2D_CLAWPATCH_OPTIONS_H */
