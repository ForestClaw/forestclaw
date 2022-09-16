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

#ifndef DISK_USER_H
#define DISK_USER_H

#include "../all/advection_user.h"

#ifdef __cplusplus
extern "C"
{
#endif


typedef struct user_options
{
    int example;
    double alpha;
    int claw_version;

    double maxelev;   

    int is_registered;

} user_options_t;


user_options_t* disk_options_register (fclaw_app_t * app,
                                       const char *configfile);

void disk_options_store (fclaw2d_global_t* glob, user_options_t* user_opt);

const user_options_t* disk_get_options(fclaw2d_global_t* glob);

void disk_link_solvers(fclaw2d_global_t *glob);

void disk_global_post_process(fclaw_options_t *fclaw_opt,
                              fclaw3dx_clawpatch_options_t *clawpatch_opt,
                              user_options_t *user_opt);

void disk_map_extrude(fclaw2d_map_context_t *cont, 
                      const double maxelev);



#ifdef __cplusplus
}
#endif

#endif
