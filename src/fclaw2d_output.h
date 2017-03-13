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

#ifndef FCLAW2D_OUTPUT_H
#define FCLAW2D_OUTPUT_H

#include <fclaw2d_vtable.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

void fclaw2d_output_frame (fclaw2d_global_t * glob, int iframe);

/* --------------------------------------------------------------------
 * Write a one-stop function that creates a big VTK file in parallel.
 * The file extension is basename.vtu.  It contains per-patch data for
 * mpirank, blockno, and global patchno, and the meqn-vector data.
 * Function will be refactored into our callback structure eventually.
 * -------------------------------------------------------------------- */

void fclaw2d_output_write_vtk (fclaw2d_global_t *glob, const char *basename);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
