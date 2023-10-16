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

#ifndef SWIRL_USER_H
#define SWIRL_USER_H

#include "../all/advection_user.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct user_options
{
    int example;
    int use_claw3d;

    double alpha;

    double maxelev;    /* Extruded height */

    double *center;
    const char* center_string;

    int claw_version;
    int is_registered;

} user_options_t;

#ifdef P8HACK

void swirl_link_solvers(fclaw_global_t *glob);

#endif /* P8HACK */

/* ------------------------------------- Options ---------------------------------------*/

user_options_t* swirl_options_register (fclaw_app_t * app,
                                        const char *configfile);

void swirl_options_store (fclaw3d_global_t* glob, user_options_t* user);

const user_options_t* swirl_get_options(fclaw3d_global_t* glob);

#ifdef P8HACK

#define SWIRL_SETAUX_MANIFOLD FCLAW_F77_FUNC(swirl_setaux_manifold, \
                                               SWIRL_SETAUX_MANIFOLD)

void SWIRL_SETAUX_MANIFOLD(const int* mbc,
                            const int* mx, const int* my, const int*mz,
                            const double* xlower, const double* ylower,
                            const double* zlower,
                            const double* dx, const double* dy,
                            const double *dz,
                            const int* maux, double aux[],
                            const int* blockno,
                            double xrot[], double yrot[], double zrot[],
                            double volume[],double faceareas[]);


#define SWIRL_SET_VELOCITY_MANIFOLD FCLAW_F77_FUNC(swirl_set_velocity_manifold, \
                                               SWIRL_SET_VELOCITY_MANIFOLD)


void SWIRL_SET_VELOCITY_MANIFOLD(const int* mx, const int* my, const int*mz,
                                 const int* mbc,
                                 const double* dx, const double* dy,
                                 const double *dz,
                                 const double *t, 
                                 const int *blockno,
                                 const double* xlower, const double* ylower,
                                 const double* zlower,
                                 double xrot[], double yrot[], double zrot[],
                                 double faceareas[],
                                 double aux[], const int* maux);

void swirl_map_extrude(fclaw2d_map_context_t *cont, double maxelev);

#endif /* P8HACK */

#ifdef __cplusplus
}
#endif

#endif
