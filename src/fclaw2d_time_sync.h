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


#ifndef FCLAW2D_TIME_SYNC_H
#define FCLAW2D_TIME_SYNC_H

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


struct fclaw_global;

/** 
 *  @file
 *  Routines that add corrections from fine grids to coarse grids when timesteps are synced.
 */

/**
 * @brief The type of reset to perform
 */
typedef enum fclaw2d_time_sync_type
{   /** Reset registers at coarse-fine boundaries */
    FCLAW2D_TIME_SYNC_RESET_F2C = 1,   
    /** Reset registers between same size grids */
    FCLAW2D_TIME_SYNC_RESET_SAMESIZE,  
    /** Reset registers at physical boundary */
    FCLAW2D_TIME_SYNC_RESET_PHYS             
} fclaw_time_sync_type_t;


/**
 * @brief Add corrections from  fine grids to coarse grid.  
 * This is is done when two or more levels are time synchronized.
 *	   - All parallel patches are valid 
 *	   - Iterate over boundary patches only, since correction occurs
 *	     only at coarse/fine boundaries.
 * 
 * @param glob the global context
 * @param minlevel the minimum level
 * @param maxlevel the maximmum level
 */
void fclaw2d_time_sync(struct fclaw_global *glob, int minlevel, int maxlevel);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif

