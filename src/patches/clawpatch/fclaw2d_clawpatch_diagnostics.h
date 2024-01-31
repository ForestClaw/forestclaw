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

#ifndef FCLAW2D_CLAWPATCH_DIAGNOSTICS_H
#define FCLAW2D_CLAWPATCH_DIAGNOSTICS_H

/** 
 *  @file
 *  Clawpatch diagnostics related structures and routines
 */

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

struct fclaw_global;
struct fclaw_patch;

/**
 * @brief Data structure for default diagnostic routines
 */
typedef struct {
    double* local_error;  /**< meqn x 3 array of errors on a patch */
    double* global_error; /**< meqn x 3 array of errors after gather */
    double area; /**< the area */
    double *mass; /**< the mass per eqn */
    double *mass0;  /**< Mass at initial time per eqn */
    double *c_kahan; /**< c value in Kahan summation algorithm per eqn */
} error_info_t;

/*------------------------------------ typedefs -----------------------------*/

/*-------------------------------- virtualized functions ----------------------------*/

/**
 * @brief Initialize a user defined data structer
 * 
 * @param glob the global context
 * @param patch_acc the data structure
 */
void fclaw2d_clawpatch_diagnostics_initialize(struct fclaw_global *glob,
                                              void** patch_acc);

/**
 * @brief Iterate over patches and perform computations
 * 
 * @param glob the global context
 * @param patch_acc the user defined data structure
 */
void fclaw2d_clawpatch_diagnostics_compute(struct fclaw_global *glob,
                                           void* patch_acc);

/**
 * @brief Perform final computations after fclaw2d_clawpatch_diagnostics_computer() has been called
 * 
 * @param glob the glboal context
 * @param patch_acc the user defined data structure
 * @param init_flag true if in init stage
 */
void fclaw2d_clawpatch_diagnostics_gather(struct fclaw_global *glob,
                                          void* patch_acc, int init_flag);

/**
 * @brief reset the user defined data structure
 * 
 * @param glob the global context
 * @param patch_acc the user defined data structure
 */
void fclaw2d_clawpatch_diagnostics_reset(struct fclaw_global *glob,
                                         void* patch_acc);

/**
 * @brief Deallocate the user defined data structure
 * 
 * @param glob the global context
 * @param patch_acc the user defined data structure
 */
void fclaw2d_clawpatch_diagnostics_finalize(struct fclaw_global *glob,
                                            void** patch_acc);

/**
 * @brief Initialize a global vtable
 * 
 * fclaw2d_diagnotistcs_vtable has to be initialized first
 * 
 * @param glob the global context
 */
void fclaw2d_clawpatch_diagnostics_vtable_initialize(struct fclaw_global* glob);

/**
 * @brief Calls the function in fclaw2d_clawpatch_vtable.fort_conservation_check
 * 
 * @param glob the global context
 * @param patch the patch context
 * @param blockno the block number
 * @param patchno the patch number
 * @param error_data user data pointer
 */
void fclaw2d_clawpatch_diagnostics_cons_default(struct fclaw_global *glob,
                                                struct fclaw_patch *patch,
                                                int blockno,
                                                int patchno,
                                                void *error_data);

/**
 * @brief Calls the function in fclaw2d_clawpatch_vtable.fort_compute_patch_error
 * 
 * @param glob the global context
 * @param patch the patch context
 * @param blockno the block number
 * @param patchno the patch number
 * @param error_data user data pointer
 */
void fclaw2d_clawpatch_diagnostics_error_default(struct fclaw_global *glob,
                                                 struct fclaw_patch *patch,
                                                 int blockno,
                                                 int patchno,
                                                 void *error_data);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
