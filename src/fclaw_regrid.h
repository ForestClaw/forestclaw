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

#ifndef FCLAW_REGRID_H
#define FCLAW_REGRID_H

#include <forestclaw.h>    /* Needed to define fclaw_patch_relation_t */

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

struct fclaw_global;
struct fclaw_domain;
struct fclaw_patch;

/* Called from both fclaw2d_initialize and fclaw2d_regrid */
void cb_fclaw_regrid_tag4refinement(struct fclaw_domain *domain,
                                    struct fclaw_patch *this_patch,
                                    int this_block_idx,
                                    int this_patch_idx,
                                    void *user);

void cb_fclaw_regrid_tag4coarsening(fclaw_domain_t *domain,
                                   fclaw_patch_t *fine_patches,
                                   int blockno, int fine0_patchno,
                                   void *user);


void cb_fclaw_regrid_repopulate(struct fclaw_domain * old_domain,
                                struct fclaw_patch * old_patch,
                                struct fclaw_domain * new_domain,
                                struct fclaw_patch * new_patch,
                                fclaw_patch_relation_t newsize,
                                int blockno,
                                int old_patchno,
                                int new_patchno,
                                void *user);

void fclaw_regrid_set_neighbor_types(struct fclaw_global *glob);

void fclaw_regrid(struct fclaw_global *glob);

void fclaw_after_regrid(struct fclaw_global *glob);

void fclaw_regrid_set_neighbor_types(struct fclaw_global *glob);




#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
