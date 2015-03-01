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

#ifndef AMR_REGRID_H
#define AMR_REGRID_H

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

#if 0
typedef struct fclaw2d_regrid_functions fclaw2d_regrid_functions_t;

/* ------------------------------------------------------------------------- */
// Function pointers for regridding functions

void patch_average2coarse(fclaw2d_domain_t *domain, fclaw2d_patch_t *fine_siblings,
                          fclaw2d_patch_t *coarse_patch, int blockno, int fine_patchno,
                          int coarse_patchno);

void patch_interpolate2fine(fclaw2d_domain_t* domain, fclaw2d_patch_t *coarse_patch,
                            fclaw2d_patch_t* fine_patch,
                            int this_blockno, int coarse_patchno,
                            int fine_patchno, int igrid);

void patch_copy2samesize(fclaw2d_domain_t* domain, fclaw2d_patch_t *old_patch,
                         fclaw2d_patch_t* new_patch, int blockno, int old_patchno,
                         int new_patchno);
#endif


void build_initial_domain(fclaw2d_domain_t* domain);
void rebuild_domain(fclaw2d_domain_t* old_domain, fclaw2d_domain_t* new_domain);

/** Repartition all patches in parallel.
 * \param [in] mode             A level for amrinit, -1 for running simulation.
 */
void repartition_domain(fclaw2d_domain_t** domain, int mode);

// void exchange_ghost_patch_data(fclaw2d_domain_t* domain,fclaw_bool time_interp);

void exchange_ghost_patch_data_levels(fclaw2d_domain_t* domain,
                                      int exchange_minlevel, int exchange_maxlevel);

void set_boundary_patch_ptrs(fclaw2d_domain_t* domain,int exchange_minlevel,
                             int exchange_maxlevel);

/* No time interpolation exchange done here */
void exchange_ghost_patch_data_all(fclaw2d_domain_t* domain);

/* From amr_advance.cpp */
void update_ghost_cells(fclaw2d_domain_t *domain);


void delete_ghost_patches(fclaw2d_domain_t *domain);


void amr_print_patches_and_procs(fclaw2d_domain_t *domain);



#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
