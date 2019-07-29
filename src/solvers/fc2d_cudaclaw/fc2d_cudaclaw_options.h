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

#ifndef FCLAW2D_CUDACLAW_OPTIONS_H
#define FCLAW2D_CUDACLAW_OPTIONS_H

#include <fclaw_base.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

struct fclaw2d_global;
typedef struct fc2d_cudaclaw_options fc2d_cudaclaw_options_t;

struct fc2d_cudaclaw_options
{
    int mwaves;

    const char *order_string;
    int *order;

    int *mthlim;
    const char *mthlim_string;

    int *mthbc;
    const char *mthbc_string;

    int method[7];
    int mcapa;
    int src_term;
    int use_fwaves;

    /* Output */
    int ascii_out;
    int vtk_out;

    int buffer_len;

    int is_registered;
};

fclaw_exit_type_t fc2d_cudaclaw_postprocess (fc2d_cudaclaw_options_t *
                                               clawopt);

fclaw_exit_type_t fc2d_cudaclaw_check (fc2d_cudaclaw_options_t * clawopt);

void fc2d_cudaclaw_reset (fc2d_cudaclaw_options_t * clawopt);

fc2d_cudaclaw_options_t*  fc2d_cudaclaw_options_register (fclaw_app_t * app,
                                                              const char *configfile);

void fc2d_cudaclaw_package_register(fclaw_app_t* app,
                                      fc2d_cudaclaw_options_t* clawopt);

fc2d_cudaclaw_options_t* fc2d_cudaclaw_get_options(struct fclaw2d_global *glob);

void fc2d_cudaclaw_options_store (struct fclaw2d_global* glob, fc2d_cudaclaw_options_t* clawopt);

void fc2d_cudaclaw_output(struct fclaw2d_global *glob, int iframe);

int cudaclaw_check_parameters(int mwaves);
void cudaclaw_set_method_parameters(int order[], int mthlim[], int mwaves);



#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
