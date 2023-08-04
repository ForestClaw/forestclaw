/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#ifndef SLOSH_USER_H
#define SLOSH_USER_H

#include <fclaw2d_include_all.h>

#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>

#include <fc2d_geoclaw.h>
#include <fc2d_geoclaw_options.h>


#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix Syntax highlighting */
#endif


typedef struct slosh_user_options
{
  int is_registered;
} slosh_user_options_t;

void slosh_link_solvers(fclaw2d_global_t *glob);

slosh_user_options_t* slosh_options_register (fclaw_app_t * app,
                                        const char *section,
                                        const char *configfile);

void slosh_options_store (fclaw2d_global_t* glob, slosh_user_options_t* user);

slosh_user_options_t* slosh_get_options(fclaw2d_global_t* glob);

fclaw2d_domain_t* slosh_create_domain(sc_MPI_Comm mpicomm, fclaw_options_t* gparms);

void slosh_run_program(fclaw2d_global_t* glob);

#ifdef __cplusplus
}
#endif

#endif
