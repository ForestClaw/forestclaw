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

#ifndef FCLAW_CLAWPATCH_OUTPUT_HDF_H
#define FCLAW_CLAWPATCH_OUTPUT_HDF_H

/** 
 * @file
 * Routines for hdf5 output 
 */

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

struct fclaw_global;
struct fclaw_patch;

/** 
 * Callback to access/compute patch data for visualization.
 * @param[in] global the global context
 * @param[in] this_patch the patch context
 * @param[in] this_block_idx the block index
 * @param[in] this_patch_ids the patch index
 * @param[in,out] a    The callback should write into this memory.
 *                     For coordinate locations, it holds
 *                     (my + 1) * (mx + 1) * 3 floats.
 *                     For patch data, it holds my * mx * meqn floats.
 *                     The vector index changes fastest, then mx, then my
 *                     slowest.
 */
typedef void (*fclaw_hdf_patch_data_t) (struct fclaw_global * glob,
                                        struct fclaw_patch * this_patch,
                                        int this_block_idx, int this_patch_idx,
                                        char *a);
	
/**
 * @brief Output vtu file
 * 
 * @param glob the global context
 * @param filename the filename to output to
 * @param coordinate_cb callback for coordinate data, set to NULL for default callback
 * @param value_cb callback cell data, set to NULL for default callback
 */
void fclaw_clawpatch_output_hdf5_to_file (struct fclaw_global* glob, 
                                         const char* filename,
                                         fclaw_hdf_patch_data_t coordinate_cb,
                                         fclaw_hdf_patch_data_t value_cb);


/**
 * @brief Output vtu file
 * 
 * @param glob the global context
 * @param iframe the the frame index
 */
void fclaw_clawpatch_output_hdf5 (struct fclaw_global* glob, int iframe);


#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW2D_VTK_H */
