/*
  Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun
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

#ifndef OVERPRESSURE_USER_H
#define OVERPRESSURE_USER_H

#include <fclaw2d_include_all.h>

#include <fclaw3dx_clawpatch_options.h>
#include <fclaw3dx_clawpatch.h>

#include <fc3d_clawpack46.h>
#include <fc3d_clawpack46_options.h>

#include "../all/euler3d_user.h"


#ifdef __cplusplus
extern "C"
{
#endif

typedef struct user_options
{
    int example;
    int mapping; 
    int init_choice; 

    double gamma;
    double x0;
    double y0;
    double z0;
    double r0;
    double rhoin;
    double rhoout;
    double pin;
    double pout;

    double min_z;
    double mid_z;
    double max_z;
    double scale_bump;

    const char *latitude_string;
    double *latitude;

    const char *longitude_string;
    double *longitude;

    double maxelev;

    int claw_version;

    int is_registered;
} user_options_t;


void overpressure_link_solvers(fclaw2d_global_t *glob);

user_options_t* overpressure_options_register (fclaw_app_t * app,
                                          const char *configfile);

void overpressure_options_store (fclaw2d_global_t* glob, user_options_t* user);

user_options_t* overpressure_get_options(fclaw2d_global_t* glob);

void overpressure_map_extrude(fclaw2d_map_context_t* cont,
                              const double maxelev,
                              const int mapping,
                              const double minz,
                              const double maxz,
                              const double midz,
                              const double scale_bump);


/* ------------------------------- Fortran routines ----------------------------------- */

#ifdef __cplusplus
}
#endif

#endif
