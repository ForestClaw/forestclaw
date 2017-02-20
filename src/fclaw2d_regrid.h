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

#include "forestclaw2d.h"
#include <fclaw2d_regrid_default_fort.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


void cb_fclaw2d_regrid_tag4refinement(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int this_block_idx,
                                      int this_patch_idx,
                                      void *user);

void cb_fclaw2d_regrid_repopulate(fclaw2d_domain_t * old_domain,
                                  fclaw2d_patch_t * old_patch,
                                  fclaw2d_domain_t * new_domain,
                                  fclaw2d_patch_t * new_patch,
                                  fclaw2d_patch_relation_t newsize,
                                  int blockno,
                                  int old_patchno,
                                  int new_patchno,
                                  void *user);

void fclaw2d_regrid_set_neighbor_types(fclaw2d_domain_t *domain);



void fclaw2d_regrid(fclaw2d_domain_t **domain);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
