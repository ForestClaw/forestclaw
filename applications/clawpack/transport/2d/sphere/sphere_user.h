/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun
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

#ifndef SPHERE_USER_H
#define SPHERE_USER_H

#include "../all/transport_user.h"

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
    int mapping;
    int initial_condition;
    int refine_pattern;

    double *omega;
    const char* omega_string;

    int claw_version;

    int is_registered;

} user_options_t;

struct fclaw2d_global;
struct fclaw2d_patch;

#if 0
/* So syntax highlighting works */
#endif

void sphere_link_solvers(struct fclaw2d_global *glob);

/* ---------------------------------- Options ----------------------------------------- */

const user_options_t* sphere_get_options(struct fclaw2d_global* glob);

void sphere_options_store (fclaw2d_global_t* glob, user_options_t* user);

user_options_t* sphere_options_register (fclaw_app_t * app, const char *configfile);

/* ---------------------------------- Compute Error ----------------------------------- */

#define SPHERE46_COMPUTE_ERROR FCLAW_F77_FUNC(sphere46_compute_error, \
                                              SPHERE46_COMPUTE_ERROR)

void SPHERE46_COMPUTE_ERROR(int* blockno, int *mx, int *my, int* mbc, int* meqn,
                           double *dx, double *dy, double *xlower,
                           double *ylower, double *t, double q[],
                           double error[], double soln[]);

#define SPHERE5_COMPUTE_ERROR FCLAW_F77_FUNC(sphere5_compute_error, \
                                              SPHERE5_COMPUTE_ERROR)

void SPHERE5_COMPUTE_ERROR(int* blockno, int *mx, int *my, int* mbc, int* meqn,
                           double *dx, double *dy, double *xlower,
                           double *ylower, double *t, double q[],
                           double error[], double soln[]);


/* -------------------------------------- Output -------------------------------------- */

#define  SPHERE46_FORT_WRITE_FILE FCLAW_F77_FUNC(sphere46_fort_write_file,  \
                                                SPHERE46_FORT_WRITE_FILE)
void  SPHERE46_FORT_WRITE_FILE(char* matname1,
                             int* mx,        int* my,
                             int* meqn,      int* mbc,
                             double* xlower, double* ylower,
                             double* dx,     double* dy,
                             double q[],     double error[], double soln[],
                             double *time,
                             int* patch_num, int* level,
                             int* blockno,   int* mpirank);

#define SPHERE46_FORT_HEADER_ASCII \
         FCLAW_F77_FUNC(sphere46_fort_header_ascii, \
                        SPHERE46_FORT_HEADER_ASCII)
void SPHERE46_FORT_HEADER_ASCII(char* matname1, char* matname2,
                               double* time, int* meqn, int* maux, 
                               int* ngrids);

#define  SPHERE5_FORT_WRITE_FILE FCLAW_F77_FUNC(sphere5_fort_write_file,  \
                                                SPHERE5_FORT_WRITE_FILE)
void  SPHERE5_FORT_WRITE_FILE(char* matname1,
                             int* mx,        int* my,
                             int* meqn,      int* mbc,
                             double* xlower, double* ylower,
                             double* dx,     double* dy,
                             double q[],     double error[], double soln[],
                             double *time,
                             int* patch_num, int* level,
                             int* blockno,   int* mpirank);

#define SPHERE5_FORT_HEADER_ASCII \
         FCLAW_F77_FUNC(sphere5_fort_header_ascii, \
                        SPHERE5_FORT_HEADER_ASCII)
void SPHERE5_FORT_HEADER_ASCII(char* matname1, char* matname2,
                               double* time, int* meqn, int* maux, 
                               int* ngrids);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
