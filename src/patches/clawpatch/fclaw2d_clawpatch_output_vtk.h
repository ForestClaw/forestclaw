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

#ifndef FCLAW2D_CLAWPATCH_OUTPUT_VTK_H
#define FCLAW2D_CLAWPATCH_OUTPUT_VTK_H


#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

struct fclaw2d_global;
struct fclaw2d_patch;

/** Callback to access/compute patch data for visualization.
 * \param [in,out] a    The callback should write into this memory.
 *                      For coordinate locations, it holds
 *                      (my + 1) * (mx + 1) * 3 floats.
 *                      For patch data, it holds my * mx * meqn floats.
 *                      The vector index changes fastest, then mx, then my
 *                      slowest.
 */
typedef void (*fclaw2d_vtk_patch_data_t) (struct fclaw2d_global * glob,
                                          struct fclaw2d_patch * this_patch,
                                          int this_block_idx, int this_patch_idx,
                                          char *a);
	
/** Write a file in VTK format for the whole domain in parallel.
 * \param [in] vtkspace     Relative width of visual separation of patches.
 *                          Between 0. (none) and 1. (patch width becomes 0).
 * \param [in] vtkwrite     Mode of writing:  (unused; uses 1 by default)
 *                          0 for MPI_File_write_all (faster),
 *                          1 for MPI_File_write (less memory usage).
 * \return          0 if successful, negative otherwise.
 *                  Collective with identical value on all ranks.
 */
int
fclaw2d_vtk_write_file (struct fclaw2d_global * glob, const char *basename,
                        int mx, int my,
#if FCLAW2D_PATCHDIM == 3
                        int mz,
#endif
                        int meqn,
                        double vtkspace, int vtkwrite,
                        fclaw2d_vtk_patch_data_t coordinate_cb,
                        fclaw2d_vtk_patch_data_t value_cb);


void fclaw2d_clawpatch_output_vtk (struct fclaw2d_global* glob, int iframe);


#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW2D_VTK_H */
