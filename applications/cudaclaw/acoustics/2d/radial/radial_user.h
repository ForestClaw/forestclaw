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

#ifndef RADIAL_USER_H
#define RADIAL_USER_H

#include <fc2d_cudaclaw.h>


#include <fclaw2d_include_all.h>

#include <fclaw2d_clawpatch.h>

#include <fc2d_cudaclaw.h>
#include <cudaclaw_user_fort.h>

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
    double rho;
    double bulk;

    int claw_version;
    int is_registered;

} user_options_t;

      
user_options_t* radial_options_register (fclaw_app_t * app, const char *configfile);

void radial_options_store (fclaw2d_global_t* glob, user_options_t* user);

user_options_t* radial_get_options(fclaw2d_global_t* glob);

void radial_link_solvers(fclaw2d_global_t *glob);

#if 0
fclaw2d_map_context_t* fclaw2d_map_new_nomap();

fclaw2d_map_context_t* fclaw2d_map_new_pillowdisk5(const double scale[],
                                                   const double shift[],
                                                   const double rotate[],
                                                   const double alpha);
#endif                                                   
/* --------------------------------------- Cuda ----------------------------------------*/

void setprob_cuda();

void radial_assign_rpn2(cudaclaw_cuda_rpn2_t *rpn2);
void radial_assign_rpt2(cudaclaw_cuda_rpt2_t *rpt2);

#ifdef __cplusplus
}
#endif

#endif
