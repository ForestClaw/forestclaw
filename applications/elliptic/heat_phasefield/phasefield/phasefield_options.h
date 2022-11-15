/*
Copyright (c) 2019-2020 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#ifndef PHASEFIELD_OPTIONS_H
#define PHASEFIELD_OPTIONS_H

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
} phasefield_operator_types;
#endif

typedef struct phasefield_options
{
    /* Put any user options here */
    int example;

    double S;
    double alpha;
    double m;
    double xi;

    /* Initial conditions */
    double k; 
    double gamma; 
    double r0; 

    double x0;
    double y0;

    int is_registered;

} phasefield_options_t;


phasefield_options_t* phasefield_options_register (fclaw_app_t * app,
                                                   const char *section,
                                                   const char *configfile);

void phasefield_options_store(fclaw2d_global_t* glob, phasefield_options_t* user);

const phasefield_options_t* phasefield_get_options(fclaw2d_global_t* glob);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
