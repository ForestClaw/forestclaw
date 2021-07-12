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

/** \file
 *
 * Routines for handling general ForestClaw input options.
 *
 */

#ifndef FCLAW2D_CLAWPATCH_OPTIONS_H
#define FCLAW2D_CLAWPATCH_OPTIONS_H

#include <fclaw_base.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif


/* Criteria for tagging patches */
#define  FCLAW_REFINE_CRITERIA_VALUE          0
#define  FCLAW_REFINE_CRITERIA_DIFFERENCE     1
#define  FCLAW_REFINE_CRITERIA_MINMAX         2
#define  FCLAW_REFINE_CRITERIA_GRADIENT       3
#define  FCLAW_REFINE_CRITERIA_USER           4


struct fclaw2d_global;

typedef struct fclaw2d_clawpatch_options fclaw2d_clawpatch_options_t;


struct fclaw2d_clawpatch_options
{
    /* These are constant for all clawpatch's */
    int mx;
    int my;
    int maux;
    int mbc;

    int meqn;
    int rhs_fields;

    int refinement_criteria;
    sc_keyvalue_t *kv_refinement_criteria;

    /* Advanced options */
    int interp_stencil_width;
    int ghost_patch_pack_aux;

    int is_registered;

};

fclaw2d_clawpatch_options_t *
fclaw2d_clawpatch_options_register(fclaw_app_t* app, const char* configfile);

void fclaw2d_clawpatch_options_store (struct fclaw2d_global *glob, 
                                      fclaw2d_clawpatch_options_t* clawpatch_options);

fclaw2d_clawpatch_options_t* fclaw2d_clawpatch_get_options(struct fclaw2d_global* glob);


void fclaw2d_clawpatch_set_refinement_criteria(int r);

int fclaw2d_clawpatch_get_refinement_criteria();



#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* FCLAW2D_CLAWPATCH_OPTIONS_H */
