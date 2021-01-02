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

#ifndef PILLOWSPHERE_H
#define PILLOWSPHERE_H

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

struct fclaw2d_glob;
struct fclaw2d_patch;
struct fclaw2d_patch_transform;

typedef struct fclaw2d_clawpatch_pillow_vtable fclaw2d_clawpatch_pillow_vtable_t;


/* ----------------------------- Fortran typedefs ------------------------------------- */

typedef void (*pillow_fort_copy_block_corner_t)(int* mx, int* my,
                                                int* mbc, int* meqn,
                                                double qthis[], 
                                                double qneighbor[], 
                                                int* icorner,
                                                int* iblock);

typedef void  (*pillow_fort_average_block_corner_t)(int* mx, int* my, int* mbc,
                                                    int* meqn, 
                                                    int* refratio, 
                                                    double qcoarse[],
                                                    double qfine[], 
                                                    double areacoarse[], 
                                                    double areafine[],
                                                    int* a_coarse_corner,
                                                    int* blockno);

typedef void  (*pillow_fort_interpolate_block_corner_t)(int* mx, int* my, int* mbc,
                                                        int* meqn, int* refratio,
                                                        double qcoarse[],
                                                        double qfine[], 
                                                        int* icoarse_corner,
                                                        int* blockno);

/* ----------------------------- Use pillow sphere ------------------------------------ */

void fclaw2d_clawpatch_use_pillowsphere();

/* --------------------------------- Virtual table ------------------------------------ */

void fclaw2d_clawpatch_pillow_vtable_initialize(int claw_version);

struct fclaw2d_clawpatch_pillow_vtable
{
    /* Block corners */
    pillow_fort_copy_block_corner_t           fort_copy_block_corner;
    pillow_fort_average_block_corner_t        fort_average_block_corner;
    pillow_fort_interpolate_block_corner_t    fort_interpolate_block_corner;

    int is_set;
};


/* ------------------------------- Public access functions ---------------------------- */

fclaw2d_clawpatch_pillow_vtable_t* fclaw2d_clawpatch_pillow_vt();

int fclaw2d_clawpatch_pillow_vtable_is_set();


#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif    