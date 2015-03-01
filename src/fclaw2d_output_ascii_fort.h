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

#ifndef FCLAW2D_OUTPUT_ASCII_fort_H
#define FCLAW2D_OUTPUT_ASCII_fort_H

#include <fclaw2d_vtable.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


typedef void  (*fclaw2d_fort_write_header_t)(const int* iframe,const double* time,
                                             const int* meqn,const int* ngrids);

typedef void (*fclaw2d_fort_write_file_t)(const int* mx,        const int* my,
                                          const int* meqn,      const int* mbc,
                                          const double* xlower, const double* ylower,
                                          const double* dx,     const double* dy,
                                          double q[],           const int* iframe,
                                          const int* patch_num, const int* level,
                                          const int* blockno,   const int* mpirank);


#define  FCLAW2D_FORT_WRITE_HEADER FCLAW_F77_FUNC(fclaw2d_fort_write_header, \
                                                  FCLAW2D_FORT_WRITE_HEADER)
void  FCLAW2D_FORT_WRITE_HEADER(const int* iframe,const double* time,
                                const int* meqn,const int* ngrids);

#define FCLAW2D_FORT_WRITE_FILE FCLAW_F77_FUNC(fclaw2d_fort_write_file, \
                                               FCLAW2D_FORT_WRITE_FILE)
void FCLAW2D_FORT_WRITE_FILE(const int* mx,        const int* my,
                             const int* meqn,      const int* mbc,
                             const double* xlower, const double* ylower,
                             const double* dx,     const double* dy,
                             double q[],           const int* iframe,
                             const int* patch_num, const int* level,
                             const int* blockno,   const int* mpirank);






#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
