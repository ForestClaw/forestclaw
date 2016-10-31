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

#include <fclaw2d_forestclaw.h>
#include <fclaw2d_transform.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef void (*fclaw2d_patch_iterator_t) (fclaw2d_domain_t * domain, int level,
                                          fclaw2d_patch_callback_t pcb, void *user);


/* Opaque pointer */
typedef struct fclaw2d_patch_data fclaw2d_patch_data_t;

void
fclaw2d_patch_data_new(fclaw2d_domain_t* domain,
                       fclaw2d_patch_t* this_patch);
void
fclaw2d_patch_data_delete(fclaw2d_domain_t* domain,
                          fclaw2d_patch_t *patch);

struct fclaw2d_patch_data*
fclaw2d_patch_get_data(fclaw2d_patch_t* patch);


void
fclaw2d_domain_iterate_level_mthread (fclaw2d_domain_t * domain, int level,
                                      fclaw2d_patch_callback_t pcb, void *user);

void*
fclaw2d_patch_get_user_patch(fclaw2d_patch_t* patch);

/* --------------------------------------------------------------
   Routines that return information about connectivity.
   This information is obtained with each new regridding.
   ------------------------------------------------------------ */

int fclaw2d_patch_on_parallel_boundary (const fclaw2d_patch_t * patch);


void fclaw2d_patch_set_face_type(fclaw2d_patch_t *patch,int iface,
                                 fclaw2d_patch_relation_t face_type);

void fclaw2d_patch_set_corner_type(fclaw2d_patch_t *patch,int icorner,
                                   fclaw2d_patch_relation_t corner_type);

void fclaw2d_patch_set_missing_corner(fclaw2d_patch_t *patch,int icorner);

fclaw2d_patch_relation_t fclaw2d_patch_get_face_type(fclaw2d_patch_t* patch,
                                                        int iface);
fclaw2d_patch_relation_t fclaw2d_patch_get_corner_type(fclaw2d_patch_t* patch,
                                                          int icorner);

int fclaw2d_patch_corner_is_missing(fclaw2d_patch_t* patch,
                                    int icorner);

void fclaw2d_patch_neighbors_set(fclaw2d_patch_t* patch);

void fclaw2d_patch_neighbors_reset(fclaw2d_patch_t* patch);

int fclaw2d_patch_has_finegrid_neighbors(fclaw2d_patch_t *patch);

int fclaw2d_patch_on_coarsefine_interface(fclaw2d_patch_t *patch);

int* fclaw2d_patch_block_corner_count(fclaw2d_domain_t* domain,
                                      fclaw2d_patch_t* this_patch);

void fclaw2d_patch_set_block_corner_count(fclaw2d_domain_t* domain,
                                          fclaw2d_patch_t* this_patch,
                                          int icorner, int block_corner_count);


/* -----------------------------------------------------
   Ghost exchange
   ---------------------------------------------------- */
#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
