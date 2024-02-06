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

#ifndef INTERFACE_USER_H
#define INTERFACE_USER_H

#include <fclaw_include_all.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef struct user_options
{
    int example;
    double rhol;
    double cl;
    double rhor;
    double cr;

    int claw_version;

    int is_registered;

} user_options_t;

user_options_t* interface_options_register (fclaw_app_t * app,
                            const char *configfile);

void interface_options_store (fclaw_global_t* glob, user_options_t* user);

const user_options_t* interface_get_options(fclaw_global_t* glob);

void interface_problem_setup(fclaw_global_t* glob);


void interface_link_solvers(fclaw_global_t *glob);

fclaw_map_context_t* fclaw_map_new_nomap();

/* ----------------------------- Fortran code ----------------------------------------- */

#define INTERFACE_SETPROB FCLAW_F77_FUNC(interface_setprob,INTERFACE_SETPROB)
void INTERFACE_SETPROB(const double *rhol, const double* cl, 
                       const double* rhor, const double *cr);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
