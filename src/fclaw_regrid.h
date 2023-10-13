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

#ifndef FCLAW2D_REGRID_H
#define FCLAW2D_REGRID_H

#include <forestclaw2d.h>    /* Needed to define fclaw2d_patch_relation_t */

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

struct fclaw2d_global;
struct fclaw2d_domain;
struct fclaw2d_patch;

/* Called from both fclaw2d_initialize and fclaw2d_regrid */
void cb_fclaw2d_regrid_tag4refinement(struct fclaw2d_domain *domain,
                                      struct fclaw2d_patch *this_patch,
                                      int this_block_idx,
                                      int this_patch_idx,
                                      void *user);

void cb_regrid_tag4coarsening(fclaw2d_domain_t *domain,
                              fclaw2d_patch_t *fine_patches,
                              int blockno, int fine0_patchno,
                              void *user);


void cb_fclaw2d_regrid_repopulate(struct fclaw2d_domain * old_domain,
                                  struct fclaw2d_patch * old_patch,
                                  struct fclaw2d_domain * new_domain,
                                  struct fclaw2d_patch * new_patch,
                                  fclaw2d_patch_relation_t newsize,
                                  int blockno,
                                  int old_patchno,
                                  int new_patchno,
                                  void *user);

void fclaw2d_regrid_set_neighbor_types(struct fclaw2d_global *glob);

void fclaw2d_regrid(struct fclaw2d_global *glob);

void fclaw2d_after_regrid(struct fclaw2d_global *glob);

void fclaw2d_regrid_set_neighbor_types(struct fclaw2d_global *glob);




#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
