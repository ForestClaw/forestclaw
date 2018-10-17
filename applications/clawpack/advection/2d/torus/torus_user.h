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

#include <fclaw2d_include_all.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/* --------------------------
   Headers for both versions
   -------------------------- */

typedef struct user_options
{
    int example;
    double alpha;  /* Ratio of inner radius to outer radius */
    int claw_version;

    int is_registered;
}
user_options_t;

#define TORUS_SETPROB FCLAW_F77_FUNC(torus_setprob,TORUS_SETPROB)
void TORUS_SETPROB(const int* example, const double* alpha);

void torus_link_solvers(fclaw2d_global_t *glob);

void torus_problem_setup(fclaw2d_global_t *glob);

void torus_patch_setup(fclaw2d_global_t *glob,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx);


user_options_t* torus_options_register (fclaw_app_t * app,
                                       const char *configfile);


void torus_options_store (fclaw2d_global_t* glob, user_options_t* user);

const user_options_t* torus_get_options(fclaw2d_global_t* glob);

fclaw2d_map_context_t *
    fclaw2d_map_new_torus (fclaw2d_map_context_t* brick,
                           const double scale[],
                           const double shift[],
                           const double rotate[],
                           const double alpha,
                           const int example);

/* ----------------------
   Clawpack 4.6 headers
   ---------------------- */
#define TORUS46_COMPUTE_ERROR FCLAW_F77_FUNC(torus46_compute_error,TORUS46_COMPUTE_ERROR)

void TORUS46_COMPUTE_ERROR(int* blockno, int *mx, int *my, int* mbc, int* meqn,
                           double *dx, double *dy, double *xlower,
                           double *ylower, double *t, double q[],
                           double error[]);


#define TORUS46_SETAUX  FCLAW_F77_FUNC(torus46_setaux, TORUS46_SETAUX)
void TORUS46_SETAUX(const int* maxmx, const int* maxmy, const int* mbc,
                    const int* mx, const int* my,
                    const double* xlower, const double* ylower,
                    const double* dx, const double* dy,
                    const int* maux, double aux[]);

#define  TORUS46_FORT_WRITE_FILE FCLAW_F77_FUNC(torus46_fort_write_file,  \
                                                TORUS46_FORT_WRITE_FILE)
void     TORUS46_FORT_WRITE_FILE(char* matname1,
                                 int* mx,        int* my,
                                 int* meqn,      int* mbc,
                                 double* xlower, double* ylower,
                                 double* dx,     double* dy,
                                 double q[],     double *error,
                                 double *time,
                                 int* patch_num, int* level,
                                 int* blockno,   int* mpirank);

#define TORUS46_TAG4REFINEMENT FCLAW_F77_FUNC(torus46_tag4refinement, \
                                              TORUS46_TAG4REFINEMENT)
void  TORUS46_TAG4REFINEMENT(const int* mx,const int* my,
                             const int* mbc,const int* meqn,
                             const double* xlower, const double* ylower,
                             const double* dx, const double* dy,
                             const int* blockno,
                             double q[],
                             const double* tag_threshold,
                             const int* init_flag,
                             int* tag_patch);

#define  TORUS46_TAG4COARSENING FCLAW_F77_FUNC(torus46_tag4coarsening, \
                                              TORUS46_TAG4COARSENING)
void  TORUS46_TAG4COARSENING(const int* mx, const int* my,
                             const int* mbc, const int* meqn,
                             const double* xlower, const double* ylower,
                             const double* dx, const double* dy,
                             const int* blockno,
                             double q0[],double q1[],
                             double q2[],double q3[],
                             const double* tag_threshold,
                             int* tag_patch);



/* ----------------------
   Clawpack 5.x headers
   ---------------------- */

#define TORUS5_COMPUTE_ERROR FCLAW_F77_FUNC(torus5_compute_error,TORUS5_COMPUTE_ERROR)

void TORUS5_COMPUTE_ERROR(int* blockno, int *mx, int *my, int* mbc, int* meqn,
                          double *dx, double *dy, double *xlower,
                          double *ylower, double *t, double q[],
                          double error[]);

#define TORUS5_SETAUX  FCLAW_F77_FUNC(torus5_setaux,  TORUS5_SETAUX)
void TORUS5_SETAUX(const int* mbc,
                   const int* mx, const int* my,
                   const double* xlower, const double* ylower,
                   const double* dx, const double* dy,
                   const int* maux, double aux[]);


#define  TORUS5_FORT_WRITE_FILE FCLAW_F77_FUNC(torus5_fort_write_file,  \
                                                TORUS5_FORT_WRITE_FILE)
void     TORUS5_FORT_WRITE_FILE(char* matname1,
                                int* mx,        int* my,
                                int* meqn,      int* mbc,
                                double* xlower, double* ylower,
                                double* dx,     double* dy,
                                double q[],     double *error,
                                double *time,
                                int* patch_num, int* level,
                                int* blockno,   int* mpirank);

#define TORUS5_TAG4REFINEMENT FCLAW_F77_FUNC(torus5_tag4refinement, \
                                              TORUS5_TAG4REFINEMENT)
void  TORUS5_TAG4REFINEMENT(const int* mx,const int* my,
                             const int* mbc,const int* meqn,
                             const double* xlower, const double* ylower,
                             const double* dx, const double* dy,
                             const int* blockno,
                             double q[],
                             const double* tag_threshold,
                             const int* init_flag,
                             int* tag_patch);

#define  TORUS5_TAG4COARSENING FCLAW_F77_FUNC(torus5_tag4coarsening, \
                                              TORUS5_TAG4COARSENING)
void  TORUS5_TAG4COARSENING(const int* mx, const int* my,
                             const int* mbc, const int* meqn,
                             const double* xlower, const double* ylower,
                             const double* dx, const double* dy,
                             const int* blockno,
                             double q0[],double q1[],
                             double q2[],double q3[],
                             const double* tag_threshold,
                             int* tag_patch);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !TORUS_USER_H */
