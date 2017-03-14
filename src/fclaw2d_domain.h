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

#ifndef FCLAW2D_DOMAIN_H
#define FCLAW2D_DOMAIN_H

#include <fclaw_timer.h>
#include <fclaw2d_partition.h>
#include <fclaw2d_map.h>
#include <fclaw2d_patch.h>
#include <fclaw2d_global.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif



typedef struct fclaw2d_domain_data
{
    /* Debug counters and timers */
    int count_set_clawpatch, count_delete_clawpatch;
    int count_amr_advance, count_ghost_exchange, count_amr_regrid;
    int count_amr_new_domain;
    int count_single_step;
    int count_multiproc_corner;
    int count_grids_per_proc;
    int count_grids_remote_boundary;
    int count_grids_local_boundary;
    int is_latest_domain;
    fclaw2d_timer_t timers[FCLAW2D_TIMER_COUNT];

    /* Time at start of each subcycled time step */
    double curr_time;

    /* This should not be copied, but needs to be redone for every new domain */
    fclaw2d_domain_exchange_t *domain_exchange;
    fclaw2d_domain_indirect_t *domain_indirect;
} fclaw2d_domain_data_t;

void fclaw2d_domain_data_new(fclaw2d_domain_t *domain);
void fclaw2d_domain_data_delete(fclaw2d_domain_t* domain);
void fclaw2d_domain_data_copy(fclaw2d_domain_t *old_domain,
                              fclaw2d_domain_t *new_domain);

void fclaw2d_domain_setup(fclaw2d_global_t* glob,
                          fclaw2d_domain_t* new_domain);

void fclaw2d_domain_reset(fclaw2d_global_t* glob);


/* ----------------------------------------------------
   Access functions for domain member data, stored in
   domain->user.  These include :

      -- fclaw_app_t
      -- time
      -- mapping context
   --------------------------------------------------- */
fclaw2d_domain_data_t*
fclaw2d_domain_get_data(fclaw2d_domain_t *domain);

/* fclaw_app_t */
fclaw_app_t*
fclaw2d_domain_get_app(fclaw2d_domain_t* domain);

void
fclaw2d_domain_set_app(fclaw2d_domain_t* domain,
                       fclaw_app_t* app);


/* time */
void
fclaw2d_domain_set_time(fclaw2d_global_t *glob, double time);

double
fclaw2d_domain_get_time(fclaw2d_global_t *glob);


/* Mapping context */
fclaw2d_map_context_t*
fclaw2d_domain_get_map_context(fclaw2d_global_t* glob);

int
fclaw2d_domain_get_num_patches(fclaw2d_domain_t* domain);


/* Options */
const amr_options_t*
fclaw2d_forestclaw_get_options(fclaw2d_domain_t *domain);

const amr_options_t*
get_domain_parms(fclaw2d_domain_t *domain);

void*
fclaw2d_domain_get_user_options(fclaw2d_domain_t* domain);


#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif



#endif
