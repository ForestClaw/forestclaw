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

#ifndef SWIRL_USER_H
#define SWIRL_USER_H

#include <fclaw2d_forestclaw.h>
#include <fclaw2d_clawpatch.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef struct user_options
{
    int example;
    double alpha;
    double tperiod;
    int claw_version;

    amr_options_t* gparms;

    int is_registered;

} user_options_t;

void swirl_problem_setup(fclaw2d_domain_t* domain);

#define SWIRL_SETPROB FCLAW_F77_FUNC(swirl_setprob, SWIRL_SETPROB)
void SWIRL_SETPROB(double* tperiod);


#define SETAUX_MANIFOLD FCLAW_F77_FUNC(setaux_manifold, SETAUX_MANIFOLD)
void SETAUX_MANIFOLD(const int* mbc,const int* mx, const int* my,
                     const double* xlower, const double* ylower,
                     const double* dx, const double* dy,
                     const int* maux, double aux[],
                     const int* blockno,
                     double xd[], double yd[], double zd[]);

#define B4STEP2_MANIFOLD FCLAW_F77_FUNC(b4step2_manifold,B4STEP2_MANIFOLD)
void B4STEP2_MANIFOLD(const int* mx,const int* my,const int* mbc,
                      const double* dx,const double* dy,
                      const int* blockno,
                      double xd[], double yd[], double zd[],
                      const double* time,
                      const double* dt, const int* maux, double aux[]);

void swirl_patch_b4step2(fclaw2d_domain_t *domain,
                         fclaw2d_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx,
                         double t,
                         double dt);

void swirl_link_solvers(fclaw2d_domain_t *domain);

void swirl_patch_setup(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx);

double swirl_patch_update(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt);

#if 0
fclaw_bool swirl_patch_tag4refinement(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int this_block_idx, int this_patch_idx,
                                      int initflag);

fclaw_bool swirl_patch_tag4coarsening(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int blockno,
                                      int patchno);
#endif

/* Mappings */
fclaw2d_map_context_t* fclaw2d_map_new_nomap();

fclaw2d_map_context_t* fclaw2d_map_new_identity();

fclaw2d_map_context_t* fclaw2d_map_new_cart(fclaw2d_map_context_t* brick,
                                            const double scale[],
                                            const double shift[],
                                            const double rotate[]);

fclaw2d_map_context_t* fclaw2d_map_new_fivepatch(const double scale[],
                                                 const double shift[],
                                                 const double rotate[],
                                                 const double alpha);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
