
/*
Copyright (c) 2012-2020 Carsten Burstedde, Donna Calhoun, Erik Chudzik, Christiane Helzel
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

#ifndef FC2D_ACTIVEFLUX_PATCH_FORT_H
#define FC2D_ACTIVEFLUX_PATCH_FORT_H

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

#define  ACTIVEFLUX_FORT_WRITE_FILE FCLAW_F77_FUNC(activeflux_fort_write_file,  \
                                                   ACTIVEFLUX_FORT_WRITE_FILE)
void     ACTIVEFLUX_FORT_WRITE_FILE(char* matname1,
                                 int* mx,        int* my,
                                 int* meqn,      int* mbc,
                                 double* xlower, double* ylower,
                                 double* dx,     double* dy,
                                 double q[],     double error[], double soln[],
                                 double *time,
                                 int* patch_num, int* level,
                                 int* blockno,   int* mpirank);

#define ACTIVEFLUX_FORT_HEADER_ASCII \
         FCLAW_F77_FUNC(activeflux_fort_header_ascii, \
                        ACTIVEFLUX_FORT_HEADER_ASCII)

void ACTIVEFLUX_FORT_HEADER_ASCII(char* matname1, char* matname2,
                               double* time, int* meqn, int* maux, 
                               int* ngrids);

#define ACTIVEFLUX_FORT_INTERPOLATE2FINE FCLAW_F77_FUNC( \
            activeflux_fort_interpolate2fine, ACTIVEFLUX_FORT_INTERPOLATE2FINE)
void ACTIVEFLUX_FORT_INTERPOLATE2FINE(const int* mx,const int* my,const int* mbc,
                                    const int* meqn, double qcoarse[], double qfine[],
                                    double areacoarse[], double areafine[],
                                    const int* igrid,const int* manifold);

#define ACTIVEFLUX_FORT_INTERPOLATE_FACE \
              FCLAW_F77_FUNC(activeflux_fort_interpolate_face, \
                             ACTIVEFLUX_FORT_INTERPOLATE_FACE)
void ACTIVEFLUX_FORT_INTERPOLATE_FACE(const int* mx, const int* my, 
                                               const int* mbc,const int* meqn,
                                               double qcoarse[],double qfine[],
                                               const int* idir, const int* iside,
                                               const int* num_neighbors,
                                               const int* refratio, const int* igrid,
                                               struct fclaw2d_patch_transform_data** 
                                               transform_cptr);

#define ACTIVEFLUX_FORT_INTERPOLATE_CORNER \
      FCLAW_F77_FUNC(activeflux_fort_interpolate_corner, \
                     ACTIVEFLUX_FORT_INTERPOLATE_CORNER)
void ACTIVEFLUX_FORT_INTERPOLATE_CORNER(const int* mx, const int* my, 
                                                 const int* mbc,const int* meqn, 
                                                 const int* a_refratio, double this_q[],
                                                 double neighbor_q[], const int* a_corner,
                                                 struct fclaw2d_patch_transform_data** 
                                                 transform_cptr);






#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif



