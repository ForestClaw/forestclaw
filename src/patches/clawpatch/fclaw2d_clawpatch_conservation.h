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

#ifndef FCLAW2D_CLAWPATCH_CONSERVATION_H
#define FCLAW2D_CLAWPATCH_CONSERVATION_H

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

struct fclaw2d_patch_transform_data;
struct fclaw2d_global;
struct fclaw2d_patch;

typedef struct fclaw2d_clawpatch_registers fclaw2d_clawpatch_registers_t;

struct fclaw2d_clawpatch_registers
{
	/* Two 1d arrays stored on each of four faces */
	double *edge_fluxes[4];

	/* Scaling factors */
	double *edgelengths[4];
	double *area[4];

	/* 1d arrays stored on left/right faces (fp/fm) and top/bottom faces (gp/gm) */
	double *fp[2];
	double *fm[2];
	double *gp[2];
	double *gm[2];
};

struct fclaw2d_global;
struct fclaw2d_patch;

void fclaw2d_clawpatch_time_sync_new(struct fclaw2d_global* glob,
                                     struct fclaw2d_patch* this_patch,
                                     int blockno,int patchno,
                                     fclaw2d_clawpatch_registers_t **registers);

void fclaw2d_clawpatch_time_sync_delete(fclaw2d_clawpatch_registers_t **registers);

void fclaw2d_clawpatch_time_sync_setup(struct fclaw2d_global* glob,
                                       struct fclaw2d_patch* this_patch,
                                       int blockno,int patchno);

void fclaw2d_clawpatch_time_sync_f2c(struct fclaw2d_global* glob,
                                     struct fclaw2d_patch* coarse_patch,
                                     struct fclaw2d_patch* fine_patch,
                                     int idir,
                                     int igrid,
                                     int iface_coarse,
                                     int time_interp,
                                     struct fclaw2d_patch_transform_data
                                     *transform_data);
	
void fclaw2d_clawpatch_time_sync_samesize(struct fclaw2d_global* glob,
                                          struct fclaw2d_patch* this_patch,
                                          struct fclaw2d_patch* neighbor_patch,
                                          int this_iface,int idir,
                                          struct fclaw2d_patch_transform_data
                                          *transform_data);

void fclaw2d_clawpatch_time_sync_reset(struct fclaw2d_global* glob,
                                       struct fclaw2d_patch *this_patch,
                                       int coarse_level,
                                       int reset_mode);

#if 0
void fclaw2d_clawpatch_update_cons_metric(struct fclaw2d_global* glob,
										  struct fclaw2d_patch* this_patch,
										  int blockno,int patchno);
#endif


#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW2D_CLAWPATCH_CONSERVATION_H */

