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

#ifndef AMR_UTILS_H
#define AMR_UTILS_H

#include <iostream>
#include <cstdlib>
#include <vector>

#include "fclaw2d_convenience.h"
#include "forestclaw2d.h"
#include "amr_options.h"

#include "fclaw_typedefs.h"
#include "fclaw_defs.H"

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

// -----------------------------------------------------------------
// Check parms
// -----------------------------------------------------------------

void check_amr_parms(amr_options_t *gparms);

// -----------------------------------------------------------------
// Initialize data
// -----------------------------------------------------------------
void init_domain_data(fclaw2d_domain_t *domain);
void init_block_data(fclaw2d_block_t *block);
void init_patch_data(fclaw2d_patch_t *patch);

// -----------------------------------------------------------------
// Return pointer to user data
// -----------------------------------------------------------------
fclaw2d_domain_data_t *get_domain_data(fclaw2d_domain_t *domain);
fclaw2d_block_data_t *get_block_data(fclaw2d_block_t *block);
fclaw2d_patch_data_t *get_patch_data(fclaw2d_patch_t *patch);

// -----------------------------------------------------------------
// Set user data with user defined variables, etc.
// -----------------------------------------------------------------
void set_domain_data(fclaw2d_domain_t *domain, const amr_options_t *gparms);
void set_block_data(fclaw2d_block_t *block, const int mthbc[]);
void set_patch_data(fclaw2d_patch_t *patch, ClawPatch* cp);

// -----------------------------------------------------------------
// Some lazy helper functions that really do make things easier..
// -----------------------------------------------------------------
void allocate_user_data(fclaw2d_domain_t *domain);
const amr_options_t* get_domain_parms(fclaw2d_domain_t *domain);
void set_domain_time(fclaw2d_domain_t *domain, double time);
double get_domain_time(fclaw2d_domain_t *domain);

// Will change the name of this to 'get_clawpatch' eventually
ClawPatch* get_clawpatch(fclaw2d_patch_t *patch);
const int get_refratio(fclaw2d_domain_t *domain);

// int corners_per_patch = FCLAW_CORNERS_PER_PATCH;
const int get_corners_per_patch(fclaw2d_domain_t *domain);
const int get_faces_per_patch(fclaw2d_domain_t *domain);
const int get_siblings_per_patch(fclaw2d_domain_t *domain);
const int get_p4est_refineFactor(fclaw2d_domain_t *domain);

// Misc. routines
int num_patches(fclaw2d_domain_t *domain, int level);
int pow_int(int a, int n);


#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif
