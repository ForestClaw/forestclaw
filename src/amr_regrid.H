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

typedef struct fclaw2d_regrid_functions fclaw2d_regrid_functions_t;

/* ------------------------------------------------------------------------- */
// Function pointers for regridding functions
typedef fclaw_bool (*fclaw2d_patch_tag4refinement_t)(fclaw2d_domain_t *domain,
                                                     fclaw2d_patch_t *this_patch,
                                                     int this_block_idx, int this_patch_idx,
                                                     int initflag);

typedef fclaw_bool (*fclaw2d_patch_tag4coarsening_t)(fclaw2d_domain_t *domain,
                                                     fclaw2d_patch_t *this_patch,
                                                     int this_blockno,
                                                     int this_patchno);

typedef void (*fclaw2d_patch_interpolate2fine_t)(fclaw2d_domain_t* domain,
                                                 fclaw2d_patch_t *coarse_patch,
                                                 fclaw2d_patch_t* fine_patch,
                                                 int this_blockno, int coarse_patchno,
                                                 int fine_patchno, int igrid);

typedef void (*fclaw2d_patch_average2coarse_t)(fclaw2d_domain_t *domain,
                                               fclaw2d_patch_t *fine_siblings,
                                               fclaw2d_patch_t *coarse_patch,
                                               int blockno, int fine_patchno,
                                               int coarse_patchno);

typedef void (*fclaw2d_patch_copy2samesize_t)(fclaw2d_domain_t* domain,
                                              fclaw2d_patch_t *old_patch,
                                              fclaw2d_patch_t* new_patch,
                                              int blockno, int old_patchno,
                                              int new_patchno);

/* ------------------------------------------------------------------------- */

fclaw_bool patch_tag4refinement_dummy(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int this_block_idx, int this_patch_idx,
                                      int initflag);

fclaw_bool patch_tag4coarsening_dummy(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int this_blockno,
                                      int this_patchno);

/* ------------------------------------------------------------------------- */

void link_regrid_functions(fclaw2d_domain_t* domain,
                           fclaw2d_patch_tag4refinement_t patch_tag4refinement,
                           fclaw2d_patch_tag4coarsening_t patch_tag4coarsening);


void initialize_regrid_functions(fclaw2d_regrid_functions_t *regrid_functions);

void copy_regrid_functions(fclaw2d_regrid_functions_t *old_functions,
                           fclaw2d_regrid_functions_t *new_functions);

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

fclaw2d_regrid_functions_t* get_regrid_functions(fclaw2d_domain_t *domain);

void build_initial_domain(fclaw2d_domain_t* domain);
void rebuild_domain(fclaw2d_domain_t* old_domain, fclaw2d_domain_t* new_domain);

/** Repartition all patches in parallel.
 * \param [in] mode             A level for amrinit, -1 for running simulation.
 */
void repartition_domain(fclaw2d_domain_t** domain, int mode);
void exchange_ghost_patch_data(fclaw2d_domain_t* domain);
void delete_ghost_patches(fclaw2d_domain_t *domain);


void amr_print_patches_and_procs(fclaw2d_domain_t *domain);


/* ------------------------------------------------------------------------- */
struct fclaw2d_regrid_functions
{
    fclaw2d_patch_tag4refinement_t   f_patch_tag4refinement;
    fclaw2d_patch_tag4coarsening_t   f_patch_tag4coarsening;
    fclaw2d_patch_interpolate2fine_t f_patch_interpolate2fine;
    fclaw2d_patch_average2coarse_t   f_patch_average2coarse;
    fclaw2d_patch_copy2samesize_t    f_patch_copy2samesize;
};



#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
