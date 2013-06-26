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

#include "amr_options.h"
#include "forestclaw2d.h"

#include "amr_forestclaw.H"
#include "fclaw2d_typedefs.h"

#include "fclaw2d_capi.h"

/* Note: either we make this a C .h file, or we remove the extern "C". */

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

/* -----------------------------------------------------------------
   Initialize data
   ----------------------------------------------------------------- */
void init_domain_data(fclaw2d_domain_t *domain);
void init_block_data(fclaw2d_block_t *block);
void init_patch_data(fclaw2d_patch_t *patch);
void init_block_and_patch_data(fclaw2d_domain_t *domain);

/* -----------------------------------------------------------------
   Return pointer to user data
   ---------------------------------------------------------------- */
fclaw2d_domain_data_t *get_domain_data(fclaw2d_domain_t *domain);
fclaw2d_block_data_t  *get_block_data(fclaw2d_block_t *block);
fclaw2d_patch_data_t  *get_patch_data(fclaw2d_patch_t *patch);

/* -----------------------------------------------------------------
   Set user data with user defined variables, etc.
   ----------------------------------------------------------------- */
void set_block_data(fclaw2d_block_t *block, const int mthbc[]);
void set_patch_data(fclaw2d_patch_t *patch, ClawPatch* cp);


void copy_domain_data(fclaw2d_domain_t *old_domain, fclaw2d_domain_t *new_domain);

/* -----------------------------------------------------------------
   Link routines that don't have any other place to go...
   ----------------------------------------------------------------- */
void link_problem_setup(fclaw2d_domain_t* domain, fclaw2d_problem_setup_t f_problem_setup);
void link_patch_output(fclaw2d_domain_t* domain, fclaw2d_patch_output_t f_patch_output);

#if 0 /* moved into fclaw2d_capi.h */

/* -----------------------------------------------------------------
   Some helper functions that really do make things easier..
   ----------------------------------------------------------------- */
const amr_options_t* get_domain_parms(fclaw2d_domain_t *domain);
void set_domain_time(fclaw2d_domain_t *domain, double time);
double get_domain_time(fclaw2d_domain_t *domain);

// Misc. routines
int num_patches(fclaw2d_domain_t *domain, int level,int include_shadow);
int pow_int(int a, int n);

#endif

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif
