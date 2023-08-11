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

#ifndef FCLAW_VTABLE_H
#define FCLAW_VTABLE_H

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

/** 
 *  @file 
 *  @brief ForestClaw vtable
 */

struct fclaw_global;
struct fclaw_patch;

/* ---------------------------------- Typedefs ---------------------------------------- */  
/**
 * @brief Sets up common blocks in Fortran
 * 
 * @param glob the global context
 */
typedef void (*fclaw_problem_setup_t)(struct fclaw_global *glob);

/**
 * @brief Sets up common blocks in Fortran
 * 
 * @param glob the global context
 * @param iframe frame number
 */
typedef void (*fclaw_output_frame_t)(struct fclaw_global * glob, int iframe);

/**
 * @brief Called after each regridding
 *  
 * @param glob the global context
 */
typedef void (*fclaw_after_regrid_t)(struct fclaw_global *glob);

/* ------------------------------------ vtable ---------------------------------------- */  
/**
 * @brief vtable for general ForestClaw functions
 */
typedef struct fclaw_vtable
{
	/** @brief sets up common blocks in Fortran */
	fclaw_problem_setup_t              problem_setup;

	/** @brief called after each regridding */
	fclaw_after_regrid_t               after_regrid;

	/** @brief called for output */
	fclaw_output_frame_t               output_frame;

	/** @brief true if vtable is set */
	int is_set;

} fclaw_vtable_t;


/**
 * @brief get the fclaw2d vtable
 * 
 * @param glob the global context
 */
fclaw_vtable_t* fclaw_vt(struct fclaw_global* glob);

/**
 * @brief Initialize fclaw2d vtable
 * 
 * @param glob the global context
 */
void fclaw_vtable_initialize(struct fclaw_global* glob);

/**
 * @brief Called after each regridding
 * 
 * @param glob the global context
 */
void fclaw_after_regrid(struct fclaw_global *glob);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
