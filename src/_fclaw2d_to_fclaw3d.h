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

//Metric terms
#define fclaw2d_metric_vtable_t              fclaw3d_metric_vtable_t
#define fclaw2d_metric_vt                    fclaw3d_metric_vt
#define fclaw2d_metric_patch_scalar          fclaw3d_metric_patch_scalar
#define fclaw2d_metric_patch_new             fclaw3d_metric_patch_new
#define fclaw2d_metric_patch_delete          fclaw3d_metric_patch_delete
#define fclaw2d_metric_patch_build           fclaw3d_metric_patch_build
#define fclaw2d_metric_patch_build_from_fine fclaw3d_metric_patch_build_from_fine
#define fclaw2d_metric_vtable_initialize     fclaw3d_metric_vtable_initialize
#define fclaw2d_metric_patch_get_area        fclaw3d_metric_patch_get_volume
#define fclaw2d_metric_get_metric_patch      fclaw3d_metric_get_metric_patch
#define fclaw2d_metric_patch_compute_area    fclaw3d_metric_patch_compute_volume
#define fclaw2d_metric_vtable_initialize     fclaw3d_metric_vtable_initialize
#define fclaw2d_metric_patch_nodes_size      fclaw3d_metric_patch_nodes_size

// static names
#define metric_average_area_from_fine        metric_average_volume_from_fine

// fclaw3d_metric.hpp
#define fclaw2d_metric_patch_t               fclaw3d_metric_patch_t

//For Clawpatch
#define fclaw2d_clawpatch_get_metric_patch   fclaw3d_clawpatch_get_metric_patch
#define fclaw2d_clawpatch_metric_scalar      fclaw3d_clawpatch_metric_scalar
#define fclaw2d_clawpatch_metric_vector      fclaw3d_clawpatch_metric_basis
#define fclaw2d_clawpatch_metric_data        fclaw3d_clawpatch_metric_data
#define fclaw2d_clawpatch_get_area           fclaw3d_clawpatch_get_volume
#define clawpatch_get_area                   clawpatch_get_volume


#if 0
// ForestClaw files
#define fclaw_map_context_t                fclaw_map_context_t
#define fclaw2d_map_destroy                  fclaw2d_map_destroy
#define fclaw2d_map_query                    fclaw2d_map_query
#endif





