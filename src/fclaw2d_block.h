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

#ifndef FCLAW2D_BLOCK_H
#define FCLAW2D_BLOCK_H

#include <fclaw2d_forestclaw.h>
#include <fclaw2d_defs.h>
#include <fclaw2d_global.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef struct fclaw2d_block_data
{
    int mthbc[FCLAW2D_NUMFACES];  /* >=0 for physical bc types */
}
fclaw2d_block_data_t;


void
fclaw2d_block_data_new(fclaw2d_domain_t *domain);

fclaw2d_block_data_t*
fclaw2d_block_get_data(fclaw2d_block_t* block);

void
    fclaw2d_block_set_data(fclaw2d_block_t* block,const int mthbc[]);

void fclaw2d_block_get_block_boundary(fclaw2d_global_t * glob,
                                      fclaw2d_patch_t * patch,
                                      fclaw_bool *intersects_block);

#if 0
fclaw2d_block_data_t *get_block_data(fclaw2d_block_t *block);

void set_block_data(fclaw2d_block_t *block, const int mthbc[]);

void init_block_and_patch_data(fclaw2d_domain_t *domain);
#endif

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
