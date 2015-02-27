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

#ifndef FCLAW2D_OUTPUT_FORT_H
#define FCLAW2D_OUTPUT_FORT_H

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/* Keep theses out of files requiring C++ linkage */

/* Header files for serial output */
typedef void (*fclaw2d_output_write_tfile_t)(const int* iframe,
                                             const double* t, const int* meqn,
                                             const int *ngrids);

#define FCLAW2D_OUTPUT_WRITE_TFILE FCLAW_F77_FUNC(fclaw2d_output_write_tfile, \
                                                  FCLAW2D_OUTPUT_WRITE_TFILE)

void FCLAW2D_OUTPUT_WRITE_TFILE(const int *iframe, const double *time,
                                const int *meqn, const int *ngrids);

/* Actual file that does the writing  (for serial output) */
typedef void (*fclaw2d_output_write_qfile_t)(const int* mx,        const int* my,
                                             const int* meqn,      const int* mbc,
                                             const double* xlower, const double* ylower,
                                             const double* dx,     const double* dy,
                                             double q[],           const int* iframe,
                                             const int* patch_idx, const int* level,
                                             const int* blockno,   const int* mpirank);

/* File that does the writing */
#define FCLAW2D_OUTPUT_WRITE_QFILE FCLAW_F77_FUNC(fclaw2d_output_write_qfile, \
                                                  FCLAW2D_OUTPUT_WRITE_QFILE)
void FCLAW2D_OUTPUT_WRITE_QFILE(const int* mx,        const int* my,
                                const int* meqn,      const int* mbc,
                                const double* xlower, const double* ylower,
                                const double* dx,     const double* dy,
                                double q[],           const int* iframe,
                                const int* patch_idx, const int* level,
                                const int* blockno,   const int* mpirank);

/* Create new frame file */
#define FCLAW2D_OUTPUT_NEW_QFILE FCLAW_F77_FUNC(fclaw2d_output_new_qfile, \
                                                FCLAW2D_OUTPUT_NEW_QFILE)
void FCLAW2D_OUTPUT_NEW_QFILE(const int *iframe);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
