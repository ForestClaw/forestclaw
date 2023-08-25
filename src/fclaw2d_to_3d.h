/*
Copyright (c) 2012-2023 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#ifndef FCLAW2D_TO_3D_H
#define FCLAW2D_TO_3D_H

#ifdef FORESTCLAW2D_H
#error "The include files forestclaw2d.h and fclaw2d_to_3d.h cannot be combined"
#endif

#include <p4est_to_p8est.h>

/* redefine macros */
#define FCLAW2D_SPACEDIM                FCLAW3D_SPACEDIM
#define FCLAW2D_NUMFACES                FCLAW3D_NUMFACES
#define FCLAW2D_NUMCORNERS              FCLAW3D_NUMCORNERS
#define FCLAW2D_NUMSIBLINGS             FCLAW3D_NUMSIBLINGS
#define FCLAW2D_NUMFACENEIGHBORS        FCLAW3D_NUMFACENEIGHBORS
/* not redefining REFINEFACTOR, which should be dimension-independent */
#define FCLAW2D_FILE_USER_STRING_BYTES  FCLAW3D_FILE_USER_STRING_BYTES
#define FCLAW2D_FILE_MAX_BLOCK_SIZE     FCLAW3D_FILE_MAX_BLOCK_SIZE
#define FCLAW2D_FILE_ERR_SUCCESS        FCLAW3D_FILE_ERR_SUCCESS
#define FCLAW2D_FILE_ERR_FORMAT         FCLAW3D_FILE_ERR_FORMAT
#define FCLAW2D_FILE_ERR_NOT_IMPLEMENTED FCLAW3D_FILE_ERR_NOT_IMPLEMENTED

/* redefine typedefs */
#define fclaw2d_patch_flags_t           fclaw3d_patch_flags_t
#define fclaw2d_file_context            fclaw3d_file_context
#define fclaw2d_file_context_t          fclaw3d_file_context_t
#define fclaw2d_domain_exchange_t       fclaw3d_domain_exchange_t
#define fclaw2d_domain_indirect         fclaw3d_domain_indirect
#define fclaw2d_domain_indirect_t       fclaw3d_domain_indirect_t
#define fclaw2d_domain_indirect_t       fclaw3d_domain_indirect_t
#define fclaw2d_domain_exchange_t       fclaw3d_domain_exchange_t
#define fclaw_patch_d2_t                fclaw_patch_d3_t
#define fclaw_block_d2_t                fclaw_block_d3_t

/* redefine enums */
#define FCLAW2D_PATCH_CHILDID           FCLAW3D_PATCH_CHILDID
#define FCLAW2D_PATCH_FIRST_SIBLING     FCLAW3D_PATCH_FIRST_SIBLING
#define FCLAW2D_PATCH_ON_PARALLEL_BOUNDARY  FCLAW3D_PATCH_ON_PARALLEL_BOUNDARY
#define FCLAW2D_PATCH_IS_GHOST          FCLAW3D_PATCH_IS_GHOST
#define FCLAW2D_PATCH_ON_BLOCK_FACE_0   FCLAW3D_PATCH_ON_BLOCK_FACE_0
#define FCLAW2D_PATCH_ON_BLOCK_FACE_1   FCLAW3D_PATCH_ON_BLOCK_FACE_1
#define FCLAW2D_PATCH_ON_BLOCK_FACE_2   FCLAW3D_PATCH_ON_BLOCK_FACE_2
#define FCLAW2D_PATCH_ON_BLOCK_FACE_3   FCLAW3D_PATCH_ON_BLOCK_FACE_3

/* redefine variables */
#define fclaw2d_patch_block_face_flags  fclaw3d_patch_block_face_flags
#define fclaw2d_smallest_h              fclaw3d_smallest_h

/* redefine functions */
#define fclaw2d_check_initial_level     fclaw3d_check_initial_level
#define fclaw2d_file_open_write         fclaw3d_file_open_write
#define fclaw2d_file_write_block        fclaw3d_file_write_block
#define fclaw2d_file_write_array        fclaw3d_file_write_array
#define fclaw2d_file_open_read          fclaw3d_file_open_read
#define fclaw2d_file_read_block         fclaw3d_file_read_block
#define fclaw2d_file_read_array         fclaw3d_file_read_array
#define fclaw2d_file_error_string       fclaw3d_file_error_string
#define fclaw2d_file_close              fclaw3d_file_close
#define fclaw2d_domain_new_unitsquare   fclaw3d_domain_new_unitcube
#define fclaw2d_domain_new_brick        fclaw3d_domain_new_brick
#define fclaw2d_domain_new_conn         fclaw3d_domain_new_conn
#define fclaw2d_domain_corner_faces     fclaw3d_domain_corner_faces
#define fclaw2d_domain_indirect_begin   fclaw3d_domain_indirect_begin
#define fclaw2d_domain_indirect_neighbors fclaw3d_domain_indirect_neighbors
#define fclaw2d_domain_indirect_end     fclaw3d_domain_indirect_end
#define fclaw2d_domain_indirect_destroy fclaw3d_domain_indirect_destroy
#define fclaw2d_patch_corner_dimension  fclaw3d_patch_corner_dimension
#define fclaw2d_patch_childid           fclaw3d_patch_childid
#define fclaw2d_patch_is_first_sibling  fclaw3d_patch_is_first_sibling
#define fclaw2d_patch_is_ghost          fclaw3d_patch_is_ghost
#define fclaw2d_patch_boundary_type     fclaw3d_patch_boundary_type
#define fclaw2d_patch_normal_match      fclaw3d_patch_normal_match
#define fclaw2d_patch_face_neighbors    fclaw3d_patch_face_neighbors
#define fclaw2d_patch_face_swap         fclaw3d_patch_face_swap
#define fclaw2d_patch_face_transformation   fclaw3d_patch_face_transformation
#define fclaw2d_patch_face_transformation_block fclaw3d_patch_face_transformation_block
#define fclaw2d_patch_face_transformation_intra fclaw3d_patch_face_transformation_intra
#define fclaw2d_patch_face_transformation_valid fclaw3d_patch_face_transformation_valid
#define fclaw2d_patch_transform_face    fclaw3d_patch_transform_face
#define fclaw2d_patch_transform_face2   fclaw3d_patch_transform_face2
#define fclaw2d_patch_corner_neighbors  fclaw3d_patch_corner_neighbors
#define fclaw2d_patch_corner_swap       fclaw3d_patch_corner_swap
#define fclaw2d_patch_transform_corner  fclaw3d_patch_transform_corner
#define fclaw2d_patch_transform_corner2 fclaw3d_patch_transform_corner2
#define fclaw2d_domain_set_refinement   fclaw3d_domain_set_refinement
#define fclaw2d_patch_mark_refine       fclaw3d_patch_mark_refine
#define fclaw2d_patch_mark_coarsen      fclaw3d_patch_mark_coarsen
#define fclaw2d_domain_iterate_adapted  fclaw3d_domain_iterate_adapted
#define fclaw2d_domain_allocate_before_partition    fclaw3d_domain_allocate_before_partition
#define fclaw2d_domain_retrieve_after_partition     fclaw3d_domain_retrieve_after_partition
#define fclaw2d_domain_iterate_partitioned  fclaw3d_domain_iterate_partitioned
#define fclaw2d_domain_free_after_partition fclaw3d_domain_free_after_partition
#define fclaw2d_domain_allocate_before_exchange fclaw3d_domain_allocate_before_exchange
#define fclaw2d_domain_free_after_exchange  fclaw3d_domain_free_after_exchange
#define fclaw2d_domain_ghost_exchange   fclaw3d_domain_ghost_exchange
#define fclaw2d_domain_ghost_exchange_begin fclaw3d_domain_ghost_exchange_begin
#define fclaw2d_domain_ghost_exchange_end   fclaw3d_domain_ghost_exchange_end
#define fclaw2d_domain_free_after_exchange  fclaw3d_domain_free_after_exchange
#define fclaw2d_domain_indirect_destroy     fclaw3d_domain_indirect_destroy
#define fclaw2d_domain_indirect_end         fclaw3d_domain_indirect_end
#define fclaw2d_domain_indirect_begin       fclaw3d_domain_indirect_begin
#define fclaw2d_map_store                   fclaw3d_map_store
#define fclaw2d_map_get                     fclaw3d_map_get

#define fclaw2d_exchange_setup          fclaw3d_exchange_setup
#define fclaw2d_exchange_delete         fclaw3d_exchange_delete
#define fclaw2d_exchange_ghost_patches_begin fclaw3d_exchange_ghost_patches_begin
#define fclaw2d_exchange_ghost_patches_end fclaw3d_exchange_ghost_patches_end
#define fclaw2d_domain_serialization_enter  fclaw3d_domain_serialization_enter
#define fclaw2d_domain_serialization_leave  fclaw3d_domain_serialization_leave
#define fclaw2d_domain_is_meta          fclaw3d_domain_is_meta
#define fclaw2d_domain_init_meta        fclaw3d_domain_init_meta
#define fclaw2d_domain_destroy          fclaw3d_domain_destroy
#define fclaw2d_domain_adapt            fclaw3d_domain_adapt
#define fclaw2d_domain_partition        fclaw3d_domain_partition
#define fclaw2d_domain_partition_unchanged  fclaw3d_domain_partition_unchanged
#define fclaw2d_domain_complete         fclaw3d_domain_complete
#define fclaw2d_domain_write_vtk        fclaw3d_domain_write_vtk
#define fclaw2d_domain_list_levels      fclaw3d_domain_list_levels
#define fclaw2d_domain_list_neighbors   fclaw3d_domain_list_neighbors
#define fclaw2d_domain_list_adapted     fclaw3d_domain_list_adapted
#define fclaw2d_domain_search_points    fclaw3d_domain_search_points
#define fclaw2d_domain_iterate_cb       fclaw3d_domain_iterate_cb
#define fclaw2d_domain_iterate_level_mthread fclaw3d_domain_iterate_level_mthread

/* translations for maps */
#define fclaw2d_map_context_t           fclaw3d_map_context_t
#define fclaw2d_map_destroy_t           fclaw3d_map_destroy_t
#define fclaw2d_map_destroy             fclaw3d_map_destroy
#define fclaw2d_map_new_nomap           fclaw3d_map_new_nomap

/* translations for the global data structure */
#define fclaw2d_iterate_patch_cb        fclaw3d_iterate_patch_cb
#define fclaw2d_iterate_family_cb       fclaw3d_iterate_family_cb
#define fclaw2d_domain_integrate_rays   fclaw3d_domain_integrate_rays
#define fclaw2d_overlap_exchange        fclaw3d_overlap_exchange
#define fclaw2d_global_store_map        fclaw3d_global_store_map

/* translations not found in p4est */
#ifndef p4est_wrap_new_unitsquare
#define p4est_wrap_new_unitsquare       p8est_wrap_new_unitcube
#endif

#endif /* !FCLAW2D_TO_3D_H */
