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

#ifndef AMR_FORESTCLAW_H
#define AMR_FORESTCLAW_H


#include "fclaw2d_typedefs.h"

#include "forestclaw2d.h"

#include "subcycle_manager.H"
#include "amr_regrid.H"
#include "amr_utils.H"
#include "ClawPatch.H"
#include "fclaw2d_transform.h"

/* Note: either we make this a C .h file, or we remove the extern "C". */

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif


typedef enum fclaw2d_exchange_type
{
    FCLAW2D_COPY = 1,
    FCLAW2D_AVERAGE,
    FCLAW2D_INTERPOLATE,
} fclaw2d_exchange_type_t;

typedef enum fclaw2d_grid_type
{
    FCLAW2D_IS_COARSE = 1,
    FCLAW2D_IS_FINE,
} fclaw2d_grid_type_t;

typedef struct fclaw2d_exchange_info
{
    fclaw_bool time_interp;
    int level;
    fclaw2d_exchange_type_t exchange_type;
    fclaw2d_grid_type_t grid_type;
} fclaw2d_exchange_info_t;

// ----------------------------------------
// Routines needed to time step and produce
// output
// ----------------------------------------

// amr_advance.cpp
double advance_all_levels (fclaw2d_domain_t * domain,
                           subcycle_manager * a_time_stepper);

// amr_regrid.cpp
void regrid (fclaw2d_domain_t ** domain);

// amr_output.cpp
void amrout (fclaw2d_domain_t * domain, int iframe);
void patch_output_default (fclaw2d_domain_t * domain, fclaw2d_patch_t * patch,
                           int this_block_idx, int this_patch_idx,
                           int iframe, int num, int matlab_level);

// -------------------------------------------
// Routines needed to fill in ghost cells
// -------------------------------------------
// amr_levelexchange.cpp
// Pass in the enum value for a running timer or FCLAW2D_TIMER_NONE.
void level_exchange (fclaw2d_domain_t * domain, int a_level);

// amr_exchangecoarse.cpp
// Pass in the enum value for a running timer or FCLAW2D_TIMER_NONE.
void exchange_with_coarse (fclaw2d_domain_t * domain, int a_level,
                           double t_level, fclaw_bool time_interp);

void timeinterp(fclaw2d_domain_t *domain,
                int level, double alpha);

void update_ghost_all_levels(fclaw2d_domain_t* domain,
                             fclaw2d_timer_names_t running);


// amr_physbc.cpp
void set_phys_bc (fclaw2d_domain_t * domain, int a_level, double a_level_time,
                  fclaw_bool time_interp);


// -------------------------------------------
// Routines needed to communicate with p4est
// -------------------------------------------

#if 0
void get_block_boundary (fclaw2d_domain_t * domain,
                         int this_block_idx,
                         int this_patch_idx, bool * intersects_block);
#endif

void get_phys_boundary (fclaw2d_domain_t * domain,
                        int this_block_idx,
                        int this_patch_idx, fclaw_bool * intersects_bc);

void cb_corner_fill(fclaw2d_domain_t *domain,
                    fclaw2d_patch_t *this_patch,
                    int this_block_idx,
                    int this_patch_idx,
                    void *user);

void cb_face_fill(fclaw2d_domain_t *domain,
                  fclaw2d_patch_t *this_patch,
                  int this_block_idx,
                  int this_patch_idx,
                  void *user);

// ----------------------------
// Routines used by diagnostics
// ----------------------------

// amr_diagnostics.cpp
void check_conservation (fclaw2d_domain_t * domain);
void dump_patch (fclaw2d_domain_t * domain, int dump_patch);
void dump_last_patch (fclaw2d_domain_t * domain, int dump_patch);
void dump_time_interp_patch (fclaw2d_domain_t * domain, int dump_patch);
void dump_auxarray (fclaw2d_domain_t * domain, int dump_patchno);



#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif
