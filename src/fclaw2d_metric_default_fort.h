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

#ifndef FCLAW2D_METRIC_DEFAULT_FORT_H
#define FCLAW2D_METRIC_DEFAULT_FORT_H

#include "forestclaw2d.h"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/* Typedef'd functions that user can set to whatever they want */
typedef void (*fclaw2d_fort_setup_mesh_t)(const int* mx, const int* my, const int* mbc,
                                         const double* xlower, const double* ylower,
                                         const double* dx, const double* dy,
                                         int* blockno,
                                         double xp[], double yp[], double zp[],
                                         double xd[], double yd[], double zd[]);

typedef void (*fclaw2d_fort_compute_area_t)(const int* mx, const int* my, const int* mbc,
                                            const double* dx, const double* dy,
                                            const double* xlower, const double* ylower,
                                            const int* blockno, double area[],
                                            const int* level, const int* maxlevel,
                                            const int* refratio, const int* ghost_only);

typedef void (*fclaw2d_fort_compute_normals_t)(const int* mx, const int* my, const int* mbc,
                                               double xp[], double yp[], double zp[],
                                               double xd[], double yd[], double zd[],
                                               double xnormals[], double ynormals[]);

typedef void (*fclaw2d_fort_compute_tangents_t)(const int* mx, const int* my, const int* mbc,
                                                double xd[], double yd[], double zd[],
                                                double xtangents[], double ytangents[],
                                                double edge_lengths[]);

typedef void (*fclaw2d_fort_compute_surf_normals_t)(const int* mx, const int* my, const int* mbc,
                                                    double xnormals[],double ynormals[],
                                                    double edge_lengths[],
                                                    double curvature[],
                                                    double surfnormals[], double area[]);


typedef double (*fclaw2d_fort_aux_func_t)(double* xc, double *yc);


#define FCLAW2D_FORT_AVERAGE_AREA FCLAW_F77_FUNC(fclaw2d_fort_average_area, \
                                               FCLAW2D_FORT_AVERAGE_AREA) 
void FCLAW2D_FORT_AVERAGE_AREA(const int* mx, const int* my,
                                 const int* mbc,
                                 double areacoarse[],double areafine[],
                                 const int* igrid);

#define FCLAW2D_FORT_SETUP_MESH \
    FCLAW_F77_FUNC(fclaw2d_fort_setup_mesh,FCLAW2D_FORT_SETUP_MESH)
void FCLAW2D_FORT_SETUP_MESH(const int* mx, const int* my, const int* mbc,
                             const double* xlower, const double* ylower,
                             const double* dx, const double* dy,
                             int* blockno,
                             double xp[], double yp[], double zp[],
                             double xd[], double yd[], double zd[]);

#define FCLAW2D_FORT_COMPUTE_AREA \
    FCLAW_F77_FUNC(fclaw2d_fort_compute_area,FCLAW2D_FORT_COMPUTE_AREA)
void FCLAW2D_FORT_COMPUTE_AREA(const int* mx, const int* my, const int* mbc,
                               const double* dx, const double* dy,
                               const double* xlower, const double* ylower,
                               const int* blockno, double area[],
                               const int* quadsize, double quadstore[],
                               const int* ghost_only);

#define FCLAW2D_FORT_INTEGRATE_EXACT FCLAW_F77_FUNC(fclaw2d_fort_integrate_exact, \
                                                    FCLAW2D_FORT_INTEGRATE_EXACT)
void FCLAW2D_FORT_INTEGRATE_EXACT(int* mx,int* my,int* mbc, double* dx,double* dy,
                                  double* xlower, double* ylower, int* blockno,
                                  double area[], fclaw2d_fort_aux_func_t *f,
                                  double favg[], int* compute_avg, int* ghost_only);

#define FCLAW2D_FORT_COMPUTE_NORMALS \
    FCLAW_F77_FUNC(fclaw2d_fort_compute_normals,FCLAW2D_FORT_COMPUTE_NORMALS)
void FCLAW2D_FORT_COMPUTE_NORMALS(const int* mx, const int* my, const int* mbc,
                                  double xp[], double yp[], double zp[],
                                  double xd[], double yd[], double zd[],
                                  double xnormals[], double ynormals[]);

#define FCLAW2D_FORT_COMPUTE_TANGENTS \
    FCLAW_F77_FUNC(fclaw2d_fort_compute_tangents,FCLAW2D_FORT_COMPUTE_TANGENTS)
void FCLAW2D_FORT_COMPUTE_TANGENTS(const int* mx, const int* my, const int* mbc,
                                   double xd[], double yd[], double zd[],
                                   double xtangents[], double ytangents[],
                                   double edge_lengths[]);

#define FCLAW2D_FORT_COMPUTE_SURF_NORMALS \
    FCLAW_F77_FUNC(fclaw2d_fort_compute_surf_normals,FCLAW2D_FORT_COMPUTE_SURF_NORMALS)
void FCLAW2D_FORT_COMPUTE_SURF_NORMALS(const int* mx, const int* my, const int* mbc,
                                       double xnormals[],double ynormals[],
                                       double edge_lengths[],
                                       double curvature[],
                                       double surfnormals[], double area[]);

#define ONE FCLAW_F77_FUNC(one,ONE)
double ONE(double* xc, double *yc);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
