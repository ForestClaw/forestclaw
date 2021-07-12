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

#ifndef SQUARE_USER_H
#define SQUARE_USER_H

#include "../all/transport_user.h"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

#if 0
/* Needed so syntax highlighing works */
#endif


typedef struct user_options
{
    int example;
    int mapping;

    int initial_condition;  /* Smooth or non-smooth */

    double alpha;    /* Used by five-patch mapping */

    double *velocity;
    const char *velocity_string;

    double *center;
    const char *center_string;

    int claw_version;

    int is_registered;

} user_options_t;

struct fclaw2d_global;
struct fclaw2d_patch;

void square_link_solvers(struct fclaw2d_global *glob);

/* ---------------------------------- Options ----------------------------------------- */

const user_options_t* square_get_options(struct fclaw2d_global* glob);

void square_options_store (fclaw2d_global_t* glob, user_options_t* user);


user_options_t* square_options_register (fclaw_app_t * app,
                                       const char *configfile);

/* ---------------------------------- Compute Error ------------------------------------- */

#define SQUARE46_COMPUTE_ERROR FCLAW_F77_FUNC(square46_compute_error,SQUARE46_COMPUTE_ERROR)

void SQUARE46_COMPUTE_ERROR(int* blockno, int *mx, int *my, int* mbc, int* meqn,
                           double *dx, double *dy, double *xlower,
                           double *ylower, double *t, double q[],
                           double error[], double soln[]);


#define SQUARE5_COMPUTE_ERROR FCLAW_F77_FUNC(square5_compute_error,SQUARE5_COMPUTE_ERROR)

void SQUARE5_COMPUTE_ERROR(int* blockno, int *mx, int *my, int* mbc, int* meqn,
                           double *dx, double *dy, double *xlower,
                           double *ylower, double *t, double q[],
                           double error[], double soln[]);

/* ---------------------------------------- OUTPUT --------------------------------------- */

#define  SQUARE46_FORT_WRITE_FILE FCLAW_F77_FUNC(square46_fort_write_file,  \
                                                SQUARE46_FORT_WRITE_FILE)
void     SQUARE46_FORT_WRITE_FILE(char* matname1,
                                 int* mx,        int* my,
                                 int* meqn,      int* mbc,
                                 double* xlower, double* ylower,
                                 double* dx,     double* dy,
                                 double q[],     double error[], double soln[],
                                 double *time,
                                 int* patch_num, int* level,
                                 int* blockno,   int* mpirank);

#define SQUARE46_FORT_HEADER_ASCII \
         FCLAW_F77_FUNC(square46_fort_header_ascii, \
                        SQUARE46_FORT_HEADER_ASCII)
void SQUARE46_FORT_HEADER_ASCII(char* matname1, char* matname2,
                               double* time, int* meqn, int* maux, 
                               int* ngrids);

#define  SQUARE5_FORT_WRITE_FILE FCLAW_F77_FUNC(square5_fort_write_file,  \
                                                SQUARE5_FORT_WRITE_FILE)
void     SQUARE5_FORT_WRITE_FILE(char* matname1,
                                 int* mx,        int* my,
                                 int* meqn,      int* mbc,
                                 double* xlower, double* ylower,
                                 double* dx,     double* dy,
                                 double q[],     double error[], double soln[],
                                 double *time,
                                 int* patch_num, int* level,
                                 int* blockno,   int* mpirank);

#define SQUARE5_FORT_HEADER_ASCII \
         FCLAW_F77_FUNC(square5_fort_header_ascii, \
                        SQUARE5_FORT_HEADER_ASCII)
void SQUARE5_FORT_HEADER_ASCII(char* matname1, char* matname2,
                               double* time, int* meqn, int* maux, 
                               int* ngrids);
                         
#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
