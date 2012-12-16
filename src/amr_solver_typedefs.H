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

#ifndef AMR_SOLVER_TYPEDEFS_H
#define AMR_SOLVER_TYPEDEFS_H

#include "forestclaw2d.h"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

struct fclaw2d_level_time_data;
struct amr_options;

typedef void (*fclaw2d_readparms_t)(sc_options_t *opt,
                                    struct amr_options *gparms);

typedef void (*fclaw2d_checkparms_t)(struct amr_options *gparms);


typedef void (*fclaw2d_patch_setup_t)(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int this_block_idx,
                                      int this_patch_idx);

typedef void (*fclaw2d_patch_initialize_t)(fclaw2d_domain_t *domain,
                                           fclaw2d_patch_t *this_patch,
                                           int this_block_idx,
                                           int this_patch_idx);

typedef void (*fclaw2d_patch_physbc_t)(fclaw2d_domain *domain,
                                       fclaw2d_patch_t *this_patch,
                                       int this_block_idx,
                                       int this_patch_idx,
                                       double t,
                                       double dt,
                                       fclaw_bool *intersects_bc);


/* Interface to MOL solver */
typedef void (*fclaw2d_ode_solver_level_t)(int neqn, double q[], double t, double dt);


/* These return a max cfl number, if needed */
typedef double (*fclaw2d_ode_solver_rhs_patch_t)(fclaw2d_domain_t *domain,
                                                 fclaw2d_patch_t *this_patch,
                                                 int this_block_idx,
                                                 int this_patch_idx,
                                                 double t,
                                                 double *rhs);



typedef double (*fclaw2d_single_step_level_t)(fclaw2d_domain_t *domain,
                                            int level,
                                            double t, double dt);



typedef double (*fclaw2d_single_step_update_patch_t)(fclaw2d_domain_t *domain,
                                                     fclaw2d_patch_t *this_patch,
                                                     int this_block_idx,
                                                     int this_patch_idx,
                                                     double t,
                                                     double dt);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif
