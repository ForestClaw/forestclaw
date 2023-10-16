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

#ifndef LATLONG_USER_H
#define LATLONG_USER_H

#include "../all/advection_user.h"

#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix syntax highlighting */
#endif

typedef struct user_options
{
    int example;
    int claw_version;

    double revs_per_second;
    double maxelev;

    const char *latitude_string;
    double *latitude;

    const char *longitude_string;
    double *longitude;

    int is_registered;
}
user_options_t;

void latlong_link_solvers(fclaw_global_t *glob);

user_options_t*  latlong_options_register (fclaw_app_t * app,
                                           const char *configfile);\

void latlong_options_store (fclaw_global_t* glob, user_options_t* user);

const user_options_t* latlong_get_options(fclaw_global_t* glob);

/* Local mapping that is user defined. */
void latlong_map_extrude (fclaw2d_map_context_t* cont,
                          const double maxelev);




#ifdef __cplusplus
}
#endif

#endif /* LATLONG_USER_H */
