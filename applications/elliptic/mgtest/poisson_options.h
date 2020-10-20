/*
Copyright (c) 2019 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#ifndef POISSON_OPTIONS_H
#define POISSON_OPTIONS_H

#include <fclaw2d_include_all.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/* ------------------------------------- Options ---------------------------------------*/

#if 0
typedef enum {
    STARPATCH = 0,  // ThunderEgg solver : Uses FFT
    FIVEPOINT,      // Laplacian (no beta)
    VARPOISSON      // Variable coefficient (uses beta)
} poisson_operator_types;
#endif

typedef struct poisson_options
{
    /* Put any user options here */
    int example;
    int beta_choice;

    double alpha;
    double x0;
    double y0; 

    double a; 
    double b;

    double eps_disk;
    
#if 0    
    int patch_operator;
    sc_keyvalue_t *kv_patch_operator;
#endif    

    int m_polar;    // number of polar flowers

    double *x0_polar;
    const char* x0_polar_string;

    double *y0_polar;
    const char* y0_polar_string;

    double *r0_polar;
    const char* r0_polar_string;

    double *r1_polar;
    const char* r1_polar_string;

    int *n_polar;
    const char* n_polar_string;

    int is_registered;

} poisson_options_t;


poisson_options_t* poisson_options_register (fclaw_app_t * app,
                                           const char *configfile);

void poisson_options_store(fclaw2d_global_t* glob, poisson_options_t* user);

const poisson_options_t* poisson_get_options(fclaw2d_global_t* glob);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
