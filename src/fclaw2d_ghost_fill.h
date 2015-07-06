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

#ifndef FCLAW2D_GHOST_H
#define FCLAW2D_GHOST_H

#include <fclaw2d_vtable.h>
#include <fclaw2d_physical_bc.h>
#include <fclaw2d_forestclaw.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/* -------------------------------------------
   Routines needed to fill in ghost cells
   ------------------------------------------- */

typedef enum fclaw2d_ghost_fill_parallel_mode
{
    FCLAW2D_BOUNDARY_INTERIOR_ONLY = 0,  /* Don't read parallel patches */
    FCLAW2D_BOUNDARY_GHOST_ONLY,   /* read parallel patches */
    FCLAW2D_BOUNDARY_ALL   /* read parallel patches */
} fclaw2d_ghost_fill_parallel_mode_t;


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
    fclaw_bool read_parallel_patches;   /* before we have done a parallel exchange */
    fclaw2d_exchange_type_t exchange_type;
    fclaw2d_grid_type_t grid_type;
} fclaw2d_exchange_info_t;

void level_exchange (fclaw2d_domain_t * domain, int a_level);

void exchange_with_coarse (fclaw2d_domain_t * domain, int a_level,
                           double t_level, fclaw_bool time_interp);

void update_ghost_all_levels(fclaw2d_domain_t* domain,
                             fclaw2d_timer_names_t running);

void get_phys_boundary (fclaw2d_domain_t * domain,
                        int this_block_idx,
                        int this_patch_idx, fclaw_bool * intersects_bc);

/* In fclaw2d_corner_neighbors.cpp */
void fclaw2d_ghost_get_block_boundary(fclaw2d_domain_t * domain,
                                      fclaw2d_patch_t * patch,
                                      fclaw_bool *intersects_block);

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

void fclaw2d_ghost_update(fclaw2d_domain_t* domain,
                          int fine_level,
                          int coarse_level,
                          int time_interp,
                          fclaw2d_timer_names_t running);

void fclaw2d_ghost_copy4timeinterp(fclaw2d_domain_t* domain,
                                   int level);

#define FCLAW2D_GHOST_PACK FCLAW_F77_FUNC(fclaw2d_ghost_pack,FCLAW2D_GHOST_PACK)
void  FCLAW2D_GHOST_PACK(int *mx, int *my, int *mbc,
                         int *meqn, int *mint,
                         double qdata[], double area[],
                         double qpack[], int *psize,
                         int *packmode, int *pack_layers,
                         int *ierror);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
