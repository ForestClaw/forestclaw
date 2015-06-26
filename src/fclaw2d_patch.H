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

#ifndef FCLAW2D_PATCH_H
#define FCLAW2D_PATCH_H

#include <forestclaw2d.h>
#include <p4est_base.h>
#include <fclaw2d_block.h>
#include <fclaw2d_domain.h>

class ClawPatch;    /* Incomplete type */


#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/* Opaque pointer to patch data */
typedef struct fclaw2d_patch_data fclaw2d_patch_data_t;

void init_patch_data(fclaw2d_patch_t *patch);

void fclaw2d_patch_delete_data(fclaw2d_patch_t *patch);

fclaw2d_patch_data_t *get_patch_data(fclaw2d_patch_t *patch);

ClawPatch* fclaw2d_patch_get_cp(fclaw2d_patch_t* this_patch);
ClawPatch* fclaw2d_patch_new_cp(fclaw2d_patch_t* this_patch);
void fclaw2d_patch_delete_cp(fclaw2d_patch_t* this_patch);

void fclaw2d_patch_user_data_new(fclaw2d_domain_t* domain,
                                  fclaw2d_patch_t* this_patch);

void fclaw2d_patch_user_data_delete(fclaw2d_domain_t* domain,
                                    fclaw2d_patch_t *patch);
void fclaw2d_patch_user_data_cleanup(fclaw2d_domain_t* domain);




void init_block_and_patch_data(fclaw2d_domain_t *domain);

#if 0
static void cb_num_patches(fclaw2d_domain_t *domain,
	fclaw2d_patch_t *patch, int block_no, int patch_no, void *user)
{
  (*(int *) user)++;
}

int num_patches(fclaw2d_domain_t *domain, int level, int include_shadow)
{
    int count = 0;
    if (include_shadow == 0)
    {
        fclaw2d_domain_iterate_level(domain, level,
                                     cb_num_patches,
                                     &count);
    }
    else
    {
        // Include shadow patches
    }
    return count;
}
#endif


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
