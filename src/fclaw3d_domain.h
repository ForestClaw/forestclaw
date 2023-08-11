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
 * Domain structures and routines
 */

#ifndef FCLAW3D_DOMAIN_H
#define FCLAW3D_DOMAIN_H

#include <forestclaw3d.h>  /* Needed for domain_exchange/domain_indirect info */

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

struct fclaw3d_global;

typedef struct fclaw3d_domain_data
{
    /* Debug counters and timers */
    int count_set_patch;
    int count_delete_patch;

    fclaw3d_domain_exchange_t *domain_exchange;
    fclaw3d_domain_indirect_t *domain_indirect;

} fclaw3d_domain_data_t;

void fclaw3d_domain_data_new(struct fclaw3d_domain *domain);

void fclaw3d_domain_data_delete(struct fclaw3d_domain* domain);

void fclaw3d_domain_setup(struct fclaw3d_global* glob,
                          struct fclaw3d_domain* new_domain);

void fclaw3d_domain_reset(struct fclaw3d_global* glob);

fclaw3d_domain_data_t* fclaw3d_domain_get_data(struct fclaw3d_domain *domain);

/* OpenMP iterator (not part of forestclaw3d.h */
void fclaw3d_domain_iterate_level_mthread (struct fclaw3d_domain * domain, int level,
                                           fclaw3d_patch_callback_t pcb, void *user);

/* below are the functions needed for dimension independence */

/** safeguard value for dimension-independent domain */
#define FCLAW3D_DOMAIN_MAGIC 0x67890303

void fclaw3d_domain_iterate_cb
  (fclaw3d_domain_t * d3, fclaw_patch_t * patch,
   int blockno, int patchno, void *user);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* FCLAW3D_DOMAIN_H */
