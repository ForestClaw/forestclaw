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
/**
 * @file
 * Single step routines
 */

#ifndef FCLAW2D_UPDATE_SINGLE_STEP_H
#define FCLAW2D_UPDATE_SINGLE_STEP_H


#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

struct fclaw_global;

/**
 * @brief Buffer data for cudaclaw
 */
typedef struct fclaw2d_single_step_buffer_data
{
    /** Number of patches in buffer */
    int total_count;
    /** Current patch index in buffer */
    int iter;
    /**  Buffer pointer */
    void* user;    
} fclaw2d_single_step_buffer_data_t;


/**
 * @brief Struct for single step iteration over patches
 */
typedef struct fclaw2d_single_step_data
{
    /** The time */
    double t;
    /** The time step */
    double dt;
    /** The maxcfl */
    double maxcfl;
    /** The buffer data */
    fclaw2d_single_step_buffer_data_t buffer_data;
} fclaw2d_single_step_data_t;

/**
 * @brief Advance the level using a single explicit time step.
 * 
 * Assumption is that the call back function 'cb_single_step'
 * can update each patch independently of all the other
 * patches at the level.  A second assumption is that only
 * boundary conditions at the the current patch time are
 * needed.  These are set before entering the patch update
 * routine.
 *
 * This function is analogous to the MOL step solver
 * fclaw_mol_step.cpp in that upon return, all the patches at
 * the given level have been updated at the new time.
 * 
 * @param glob the global context
 * @param level the level to advance
 * @param t the current time
 * @param dt the time step
 * @return double the maxcfl
 */
double fclaw2d_update_single_step(struct fclaw_global *glob,
                                  int level,
                                  double t, double dt);


#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif
