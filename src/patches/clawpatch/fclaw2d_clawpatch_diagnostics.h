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

#ifndef FCLAW2D_CLAWPATCH_DIAGNOSTICS_H
#define FCLAW2D_CLAWPATCH_DIAGNOSTICS_H

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

struct fclaw2d_global;
struct fclaw2d_patch;

typedef struct {
    double* local_error;  /* meqn x 3 array of errors on a patch */
    double* global_error; /* meqn x 3 array of errors after gather */
    double area;
    double *mass;
    double *mass0;  /* Mass at initial time */
    double *c_kahan;  
} error_info_t;

/************************************* typedefs ******************************/

/********************************* virtualized functions *****************************/
void fclaw2d_clawpatch_diagnostics_initialize(struct fclaw2d_global *glob,
                                              void** patch_acc);

void fclaw2d_clawpatch_diagnostics_compute(struct fclaw2d_global *glob,
                                           void* patch_acc);

void fclaw2d_clawpatch_diagnostics_gather(struct fclaw2d_global *glob,
                                          void* patch_acc, int init_flag);

void fclaw2d_clawpatch_diagnostics_reset(struct fclaw2d_global *glob,
                                         void* patch_acc);

void fclaw2d_clawpatch_diagnostics_finalize(struct fclaw2d_global *glob,
                                            void** patch_acc);

void fclaw2d_clawpatch_diagnostics_vtable_initialize();


void fclaw2d_clawpatch_diagnostics_cons_default(struct fclaw2d_global *glob,
                                                struct fclaw2d_patch *patch,
                                                int blockno,
                                                int patchno,
                                                void *error_data);

void fclaw2d_clawpatch_diagnostics_error_default(struct fclaw2d_global *glob,
                                                 struct fclaw2d_patch *patch,
                                                 int blockno,
                                                 int patchno,
                                                 void *error_data);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
