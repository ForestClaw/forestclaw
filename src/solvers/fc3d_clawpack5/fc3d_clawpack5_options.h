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

#ifndef FC3D_CLAWPACK5_OPTIONS_H
#define FC3D_CLAWPACK5_OPTIONS_H

#include <fclaw_options.h>
#include <fclaw2d_base.h>
#include <fclaw_clawpatch3_options.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/* Only one copy of clawpack5_options for each run */
typedef struct fc3d_clawpack5_options
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
}
fc3d_clawpack5_options_t;

fclaw_exit_type_t fc3d_clawpack5_postprocess (fc3d_clawpack5_options_t *
                                               clawopt);
fclaw_exit_type_t fc3d_clawpack5_check (fc3d_clawpack5_options_t * clawopt);
void fc3d_clawpack5_reset (fc3d_clawpack5_options_t * clawopt);

/**
 * @brief Register options in SC
 * 
 * @param a the app context
 * @param section the section name
 * @param configfile the config file
 * @return fc3d_clawpack5_options_t* a newly allocated options struct
 */
fc3d_clawpack5_options_t *fc3d_clawpack5_options_register (fclaw_app_t *app,
                                                           const char *section,
                                                           const char *configfile);

#define SET_AMR_MODULE FCLAW_F77_FUNC(set_amr_module,SET_AMR_MODULE)
void SET_AMR_MODULE(const int* mwaves_in, const int* mcapa_in,
                    const int mthlim_in[], const int method_in[]);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
