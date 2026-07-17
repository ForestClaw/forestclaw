/*
Copyright (c) 2012-2026 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#ifndef FCLAW_CLAWPATCH_VTK_VTABLE_H
#define FCLAW_CLAWPATCH_VTK_VTABLE_H

/** 
 * @file
 * Routines for vtk output 
 */

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

#include <fclaw_base.h>

struct fclaw_global;
struct fclaw_patch;

/**
 * @brief Callback context for VTK patch callbacks
 */
typedef struct fclaw_vtk_cb_context
{
    int fits32; /* whether the mesh connectivity fits in 32 bits */
    int cells_per_patch; /* number of cells per patch */
} fclaw_vtk_cb_context_t;

/** 
 * Callback to access/compute patch data for visualization.
 * @param[in] global the global context
 * @param[in] patch the patch context
 * @param[in] blockno the block index
 * @param[in] patchno the patch index
 * @param[in,out] context the callback context
 * @param[in,out] buffer The callback should write into this memory.
 */
typedef void (*fclaw_vtk_patch_cb_t) (struct fclaw_global * glob,
                                      struct fclaw_patch * patch,
                                      int blockno, int patchno,
                                      fclaw_vtk_cb_context_t *context,
                                      char *buffer);

/** 
 * Callback to compute the number of elements (e.g. points or cells) for a patch.
 * @param[in] global the global context
 * @param[in] patch the patch context
 * @param[in] blockno the block index
 * @param[in] patchno the patch index
 * @param[in] context the callback context
 * @return the number of elements (e.g. points or cells) for this patch
 */
typedef size_t (*fclaw_vtk_patch_elements_cb_t) (struct fclaw_global * glob,
                                                 struct fclaw_patch * patch,
                                                 int blockno, int patchno,
                                                 fclaw_vtk_cb_context_t *context);

typedef enum {
    FCLAW_VTK_UINT8,
    FCLAW_VTK_INT32,
    FCLAW_VTK_INT64,
    FCLAW_VTK_UINT64,
    FCLAW_VTK_FLOAT32,
    FCLAW_VTK_FLOAT64,
    FCLAW_VTK_INT32_OR_64 /* use int32 if connectivity fits in 32 bits, otherwise use int64 */
} fclaw_vtk_entry_type_t;

typedef struct fclaw_clawpatch_vtk_vtable_entry
{
    const char* name;
    fclaw_vtk_entry_type_t type;
    int number_of_components; 
    size_t elements_per_patch; /* number of elements (e.g. points or cells) per patch */
    fclaw_vtk_patch_elements_cb_t elements_in_patch; /* callback to compute number of elements for a patch */
    fclaw_vtk_patch_cb_t callback;
} fclaw_clawpatch_vtk_vtable_entry_t; 
	
typedef struct fclaw_clawpatch_vtk_vtable
{
    sc_list_t* field_entries; /* entries for the field section */
    sc_list_t* point_entries; /* entries for the point section */
    sc_list_t* cell_entries;  /* entries for the cell section */
    sc_list_t* celldata_entries; /* entries for the celldata section */
} fclaw_clawpatch_vtk_vtable_t;

/**
 * @brief get the fclaw2d vtable
 * 
 * @param glob the global context
 */
fclaw_clawpatch_vtk_vtable_t* fclaw_clawpatch_vtk_vtable(struct fclaw_global* glob);

/**
 * @brief Initialize fclaw2d vtable
 * 
 * @param glob the global context
 */
void fclaw_clawpatch_vtk_vtable_initialize(struct fclaw_global* glob);
#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW_CLAWPATCH_VTK_VTABLE_H */
