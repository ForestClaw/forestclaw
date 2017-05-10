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

#ifndef FCLAW2D_CLAWPATCH_OUTPUT_H
#define FCLAW2D_CLAWPATCH_OUTPUT_H

#include <fclaw2d_global.h>
#include <fclaw2d_patch.h>
#include <fclaw2d_clawpatch_vtk.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/* ----------------------------------------------------------------
	Ascii output
   --------------------------------------------------------------- */
typedef void  (*fclaw2d_fort_header_ascii_t)(char* matname1,char* matname2,
                                             double* time, int* meqn, int* maux,
                                             int* ngrids);

/* Write out data */
typedef void (*fclaw2d_fort_output_ascii_t)(char* matname1,
                                          int* mx,        int* my,
                                          int* meqn,      int* mbc,
                                          double* xlower, double* ylower,
                                          double* dx,     double* dy,
                                          double q[],
                                          int* patch_num, int* level,
                                          int* blockno,   int* mpirank);

void cb_clawpatch_output_ascii (fclaw2d_domain_t * domain,
                                fclaw2d_patch_t * this_patch,
                                int this_block_idx, int this_patch_idx,
                                void *user);

void fclaw2d_clawpatch_output_ascii(fclaw2d_global_t* glob,int iframe);


/* ----------------------------------------------------------------
	VTK output
   --------------------------------------------------------------- */
int
fclaw2d_vtk_write_file (fclaw2d_global_t * glob, const char *basename,
                        int mx, int my, int meqn,
                        double vtkspace, int vtkwrite,
                        fclaw2d_vtk_patch_data_t coordinate_cb,
                        fclaw2d_vtk_patch_data_t value_cb);

void fclaw2d_clawpatch_output_vtk (fclaw2d_global_t* glob, int iframe);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
