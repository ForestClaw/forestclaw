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

#ifndef FCLAW2D_CAPI_H
#define FCLAW2D_CAPI_H

#include "amr_options.h"
#include "forestclaw2d.h"

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

/* -----------------------------------------------------------
   Data needed for time stepping
   ----------------------------------------------------------- */
typedef struct fclaw2d_level_time_data
{
    /* Single step data. This always has to be set. */
    double dt;
    double t_initial;
    double t_level;
    double t_coarse;

    /* Needed for explicit CFL limited schemes */
    double maxcfl;

    /* Extra data that might be needed for more complicated time stepping
     * Not always set. */
    double alpha;         /* Fraction of coarser dt completed. */
    double dt_coarse;
    int is_coarsest;
}
fclaw2d_level_time_data_t;

typedef void (*fclaw2d_level_advance_t)(fclaw2d_domain_t *domain,
                                        int level,
                                        fclaw2d_level_time_data_t *time_data);

typedef double (*fclaw2d_single_step_patch_t)(fclaw2d_domain_t *domain,
                                              fclaw2d_patch_t *this_patch,
                                              int this_block_idx,
                                              int this_patch_idx,
                                              double t,
                                              double dt);

/* -----------------------------------------------------------------
 * Some lazy helper functions that really do make things easier...
 * Defined in amr_utils.cpp
 * Need to be prefixed to clean up namespace
 * ---------------------------------------------------------------*/
void allocate_user_data(fclaw2d_domain_t *domain);
const amr_options_t* get_domain_parms(fclaw2d_domain_t *domain);
void set_domain_time(fclaw2d_domain_t *domain, double time);
double get_domain_time(fclaw2d_domain_t *domain);

/* int corners_per_patch = FCLAW_CORNERS_PER_PATCH; */
const int get_corners_per_patch(fclaw2d_domain_t *domain);
const int get_faces_per_patch(fclaw2d_domain_t *domain);
const int get_siblings_per_patch(fclaw2d_domain_t *domain);
const int get_p4est_refineFactor(fclaw2d_domain_t *domain);

/* Misc. routines */
int num_patches(fclaw2d_domain_t *domain, int level);
int pow_int(int a, int n);

/* Functions with C prototypes to use forestclaw from C code */

void fclaw2d_allocate_domain_data (fclaw2d_domain_t * domain,
                                   amr_options_t * gparms,
                                   fclaw2d_level_advance_t level_advance_cb,
                                   fclaw2d_single_step_patch_t
                                   single_step_patch_cb);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif
