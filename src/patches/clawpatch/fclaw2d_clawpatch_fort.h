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

#ifndef FCLAW2D_CLAWPATCH_FORT_H
#define FCLAW2D_CLAWPATCH_FORT_H

/* this header file must come first */
#include <fclaw2d_defs.h>

#include <forestclaw2d.h>
#include <fclaw2d_convenience.h>

#include <fclaw_options.h>
#include <fclaw2d_transform.h>

#ifdef __cplusplus
extern "C"
{          /* Beginning of extern "C" */
#if 0
}
#endif
#endif

/* ------------------------------------------------------------------------ 
  Output functions
  ------------------------------------------------------------------------ */

#if 0 
    /* output functions in clawpatch */
    fclaw2d_fort_write_ascii_header_t  fort_write_ascii_header;
    fclaw2d_fort_write_ascii_file_t    fort_write_ascii_file;
    fclaw2d_patch_callback_t           cb_write_ascii_file;
#endif    


typedef void  (*fclaw2d_fort_write_ascii_header_t)(char* matname1,char* matname2,
                                                   double* time, int* meqn,
                                                   int* ngrids);

#define  FCLAW2D_FORT_WRITE_HEADER FCLAW_F77_FUNC(fclaw2d_fort_write_header, \
                                                  FCLAW2D_FORT_WRITE_HEADER)

void     FCLAW2D_FORT_WRITE_HEADER(char* matname1, char* matname2,
                                   double* time, int* meqn, int* ngrids);


/* Write out data */
typedef void (*fclaw2d_fort_write_ascii_file_t)(char* matname1,
                                                int* mx,        int* my,
                                                int* meqn,      int* mbc,
                                                double* xlower, double* ylower,
                                                double* dx,     double* dy,
                                                double q[],
                                                int* patch_num, int* level,
                                                int* blockno,   int* mpirank);

#define  FCLAW2D_FORT_WRITE_FILE FCLAW_F77_FUNC(fclaw2d_fort_write_file, \
                                                FCLAW2D_FORT_WRITE_FILE)
void  FCLAW2D_FORT_WRITE_FILE(char* matname1,
                              int* mx,        int* my,
                              int* meqn,      int* mbc,
                              double* xlower, double* ylower,
                              double* dx,     double* dy,
                              double q[],
                              int* patch_num, int* level,
                              int* blockno,   int* mpirank);


/* ----------------------------------------------------------------------------------
   Internal boundary conditions
   ---------------------------------------------------------------------------------- */

#if 0
#define FCLAW2D_CLAWPATCH_FORT_ERROR FCLAW_F77_FUNC(fclaw2d_clawpatch_fort_error, \
                                                    FCLAW2D_CLAWPATCH_FORT_ERROR)

void FCLAW2D_CLAWPATCH_FORT_ERROR(int* blockno, int *mx, int *my, int* mbc, int* meqn,
                                  double *dx, double *dy, double *xlower,
                                  double *ylower, double *t, double q[],
                                  double error[]);


/* Not obvious that the user would want to generalize the conservation routine,
   so they are not in a 'default' file.*/

#define FCLAW2D_CLAWPATCH_FORT_CONSCHECK FCLAW_F77_FUNC(fclaw2d_clawpatch_fort_conscheck, \
                                                        FCLAW2D_CLAWPATCH_FORT_CONSCHECK)

void FCLAW2D_CLAWPATCH_FORT_CONSCHECK(int *mx, int *my, int* mbc, int* meqn,
                                      double *dx, double *dy,
                                      double* area, double *q, double* sum);

/* These are only needed if the user has not supplied their own routine to compute
   the area of the entire domain. Not obvious taht the user would want to
   generalize these routines. */
#define FCLAW2D_CLAWPATCH_FORT_AREA FCLAW_F77_FUNC(fclaw2d_clawpatch_fort_area, \
                                                   FCLAW2D_CLAWPATCH_FORT_AREA)

double FCLAW2D_CLAWPATCH_FORT_AREA(int *mx, int* my, int*mbc, double* dx,
                                   double* dy, double area[]);


#define FCLAW2D_CLAWPATCH_FORT_NORM FCLAW_F77_FUNC(fclaw2d_clawpatch_fort_norm, \
                                                   FCLAW2D_CLAWPATCH_FORT_NORM)

void FCLAW2D_CLAWPATCH_FORT_NORM(int *mx, int *my, int *mbc, int *meqn,
                                 double *dx, double *dy, double area[],
                                 double error[], double error_norm[]);

#endif


#ifdef __cplusplus
#if 0
{
#endif
}           /* end of extern "C" */
#endif

#endif
