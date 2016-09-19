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

#ifndef FCLAW2D_REGRID_DEFAULT_FORT_H
#define FCLAW2D_REGRID_DEFAULT_FORT_H

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

typedef void (*fclaw2d_fort_tag4refinement_t)(const int* mx,const int* my,
                                              const int* mbc,const int* meqn,
                                              const double* xlower, const double* ylower,
                                              const double* dx, const double* dy,
                                              const int* blockno,
                                              double q[],
                                              const double* tag_threshold,
                                              const int* init_flag,
                                              int* tag_patch);

typedef void (*fclaw2d_fort_tag4coarsening_t)(const int* mx, const int* my,
                                              const int* mbc, const int* meqn,
                                              const double* xlower, const double* ylower,
                                              const double* dx, const double* dy,
                                              const int* blockno,
                                              double q0[],double q1[],
                                              double q2[],double q3[],
                                              const double* tag_threshold,
                                              int* tag_patch);

typedef void (*fclaw2d_fort_interpolate2fine_t)(const int* mx, const int* my,
                                                const int* mbc, const int* meqn,
                                                double qcoarse[], double qfine[],
                                                double areacoarse[], double areafine[],
                                                const int* igrid, const int* manifold);

typedef void (*fclaw2d_fort_average2coarse_t)(const int* mx, const int* my,
                                              const int* mbc, const int* meqn,
                                              double qcoarse[],double qfine[],
                                              double areacoarse[],double areafine[],
                                              const int* igrid, const int* manifold);





#define FCLAW2D_FORT_AVERAGE_AREA FCLAW_F77_FUNC(fclaw2d_fort_average_area, \
                                                FCLAW2D_FORT_AVERAGE_AREA)
void FCLAW2D_FORT_AVERAGE_AREA(const int* mx, const int* my,
                               const int* mbc,
                               double areacoarse[],double areafine[],
                               const int* igrid);





/* ------------------------------------------------------------------
   Generic headers for customized user routines.  There are no object
   files associated with these, other than what the user might provide.
   ------------------------------------------------------------------- */
#define TAG4REFINEMENT FCLAW_F77_FUNC(tag4refinement,TAG4REFINEMENT)
void TAG4REFINEMENT(const int* mx,const int* my,
                    const int* mbc,const int* meqn,
                    const double* xlower, const double* ylower,
                    const double* dx, const double* dy,
                    const int* blockno,
                    double q[],
                    const double* tag_threshold,
                    const int* init_flag,
                    int* tag_patch);

#define TAG4COARSENING FCLAW_F77_FUNC(tag4coarsening,TAG4COARSENING)
void TAG4COARSENING(const int* mx, const int* my,
                    const int* mbc, const int* meqn,
                    const double* xlower, const double* ylower,
                    const double* dx, const double* dy,
                    const int* blockno,
                    double q0[],double q1[],
                    double q2[],double q3[],
                    const double* tag_threshold,
                    int* tag_patch);

#define INTERPOLATE2FINE FCLAW_F77_FUNC(interpolate2fine, INTERPOLATE2FINE)
void INTERPOLATE2FINE(const int* mx,const int* my,const int* mbc,
                      const int* meqn, double qcoarse[], double qfine[],
                      double areacoarse[], double areafine[],
                      const int* igrid,const int* manifold);

#define AVERAGE2COARSE FCLAW_F77_FUNC(average2coarse, AVERAGE2COARSE)
void AVERAGE2COARSE(const int* mx,const int* my,const int* mbc,
                    const int* meqn,
                    double qcoarse[],double qfine[],
                    double areacoarse[],double areafine[],
                    const int* igrid, const int* manifold);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
