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

#ifndef TORUS_USER_H
#define TORUS_USER_H

#include <fclaw2d_domain.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


#define SETAUX_MANIFOLD FCLAW_F77_FUNC(setaux_manifold,SETAUX_MANIFOLD)

void SETAUX_MANIFOLD(const int* maxmx, const int* maxmy, const int* mbc,
                     const int* mx, const int* my,
                     const double* xlower, const double* ylower,
                     const double* dx, const double* dy,
                     const int* maux, double aux[]);

#define B4STEP2_MANIFOLD FCLAW_F77_FUNC(b4step2_manifold,B4STEP2_MANIFOLD)
void B4STEP2_MANIFOLD(const int* maxmx, const int* maxmy, const int* mbc,
                      const int* mx, const int* my, const int* meqn,
                      double q[], const double* xlower, const double* ylower,
                      const double* dx, const double* dy,
                      const double* t, const double* dt,
                      const int* maux, double aux[]);


void torus_output_write_header(fclaw2d_domain_t* domain,
                               int iframe);

void torus_output_write_file(fclaw2d_domain_t *domain,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx, int this_patch_idx,
                              int iframe, int patch_num,int level);


void torus_link_solvers(fclaw2d_domain_t *domain);

void torus_patch_manifold_setup(fclaw2d_domain_t *domain,
                                fclaw2d_patch_t *this_patch,
                                int this_block_idx,
                                int this_patch_idx);

#define SETPROB_TORUS FCLAW_F77_FUNC(setprob_torus,SETPROB_TORUS)
void SETPROB_TORUS(const int* example);

void torus_patch_setup(fclaw2d_domain_t *domain);


fclaw2d_map_context_t* fclaw2d_map_new_nomap();

fclaw2d_map_context_t* fclaw2d_map_new_cart(fclaw2d_map_context_t* brick,
                                            const double scale[],
                                            const double shift[],
                                            const double rotate[]);


#define TORUS_COMPUTE_ERROR FCLAW_F77_FUNC(torus_compute_error,TORUS_COMPUTE_ERROR)

void TORUS_COMPUTE_ERROR(int* blockno, int *mx, int *my, int* mbc, int* meqn,
                         double *dx, double *dy, double *xlower,
                         double *ylower, double *t, double q[],
                         double error[]);


#define  TORUS_FORT_WRITE_HEADER FCLAW_F77_FUNC(torus_fort_write_header, \
                                                TORUS_FORT_WRITE_HEADER)

void     TORUS_FORT_WRITE_HEADER(char* matname1, char* matname2,
                                 double* time, int* meqn, int* ngrids);


#define  TORUS_FORT_WRITE_FILE FCLAW_F77_FUNC(torus_fort_write_file,  \
                                              TORUS_FORT_WRITE_FILE)
void     TORUS_FORT_WRITE_FILE(char* matname1,
                               int* mx,        int* my,
                               int* meqn,      int* mbc,
                               double* xlower, double* ylower,
                               double* dx,     double* dy,
                               double q[],     double *error,
                               double *time,
                               int* patch_num, int* level,
                               int* blockno,   int* mpirank);


fclaw2d_map_context_t *
    fclaw2d_map_new_torus (fclaw2d_map_context_t* brick,
                           const double scale[],
                           const double shift[],
                           const double rotate[],
                           const double alpha);

fclaw2d_map_context_t *
    fclaw2d_map_new_twisted_torus (fclaw2d_map_context_t* brick,
                                   const double scale[],
                                   const double shift[],
                                   const double rotate[],
                                   const double alpha);

fclaw2d_map_context_t *
    fclaw2d_map_new_annulus (fclaw2d_map_context_t* brick,
                             const double scale[],
                             const double shift[],
                             const double rotate[],
                             const double alpha);

fclaw2d_map_context_t *
    fclaw2d_map_new_latlong (fclaw2d_map_context_t* brick,
                             const double scale[],
                             const double shift[],
                             const double rotate[],
                             const double lat[],
                             const double longitude[],
                             const int a, const int b);



#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !TORUS_USER_H */
