#!/usr/bin/env python3

fclaw_1to2 = {

    #===== ./src/forestclaw2d.h =====
    "<forestclaw2d.h>" : "<forestclaw.h>",
    "fclaw2d_domain_t" : "fclaw_domain_t",
    "fclaw2d_block_t"  : "fclaw_block_t",
    "fclaw2d_patch_t"  : "fclaw_patch_t",
    "fclaw2d_patch"  : "fclaw_patch",

    "fclaw2d_block"    : "fclaw_block",

    "fclaw2d_domain_persist"   : "fclaw_domain_persist",
    "fclaw2d_domain_persist_t" : "fclaw_domain_persist_t",

    "fclaw2d_domain"    : "fclaw_domain",

    "fclaw2d_domain_attribute_add"    : "fclaw_domain_attribute_add",
    "fclaw2d_domain_attribute_access" : "fclaw_domain_attribute_access",
    "fclaw2d_domain_attribute_remove" : "fclaw_domain_attribute_remove",

    "fclaw2d_domain_dimension"        : "fclaw_domain_dimension",
    "fclaw2d_domain_num_faces"        : "fclaw_domain_num_faces",
    "fclaw2d_domain_num_corners"      : "fclaw_domain_num_corners",
    "fclaw2d_domain_num_face_corners" : "fclaw_domain_num_face_corners",
    "fclaw2d_domain_num_orientations" : "fclaw_domain_num_orientations",
    "fclaw2d_domain_corner_faces"     : "fclaw_domain_corner_faces",
    "fclaw2d_patch_corner_dimension"  : "fclaw_patch_corner_dimension",

    "fclaw2d_patch_childid"           : "fclaw_patch_childid",
    "fclaw2d_patch_is_first_sibling"  : "fclaw_patch_is_first_sibling",
    "fclaw2d_patch_is_ghost"          : "fclaw_patch_is_ghost",
    "fclaw2d_patch_callback_t"        : "fclaw_patch_callback_t",
    "fclaw2d_domain_iterate_level"    : "fclaw_domain_iterate_level",
    "fclaw2d_domain_iterate_patches"  : "fclaw_domain_iterate_patches",
    "fclaw2d_domain_iterate_families" : "fclaw_domain_iterate_families",
    "fclaw2d_patch_boundary_type"     : "fclaw_patch_boundary_type",
    "fclaw2d_patch_normal_match"      : "fclaw_patch_normal_match",

    "fclaw2d_face_neighbor"           : "fclaw_face_neighbor",
    "FCLAW2D_PATCH_BOUNDARY"          : "FCLAW_PATCH_BOUNDARY",
    "FCLAW2D_PATCH_HALFSIZE"          : "FCLAW_PATCH_HALFSIZE",
    "FCLAW2D_PATCH_SAMESIZE"          : "FCLAW_PATCH_SAMESIZE",
    "FCLAW2D_PATCH_DOUBLESIZE"        : "FCLAW_PATCH_DOUBLESIZE",
    "fclaw2d_patch_relation_t"        : "fclaw_patch_relation_t",

    "fclaw2d_patch_face_neighbors"            : "fclaw_patch_face_neighbors",
    "fclaw2d_patch_face_swap"                 : "fclaw_patch_face_swap",
    "fclaw2d_patch_face_transformation"       : "fclaw_patch_face_transformation",
    "fclaw2d_patch_face_transformation_block" : "fclaw_patch_face_transformation_block",
    "fclaw2d_patch_face_transformation_intra" : "fclaw_patch_face_transformation_intra",
    "fclaw2d_patch_face_transformation_valid" : "fclaw_patch_face_transformation_valid",
    "fclaw2d_patch_transform_face"            : "fclaw_patch_2d_transform_face",
    "fclaw2d_patch_transform_face2"           : "fclaw_patch_2d_transform_face2",
    "fclaw2d_patch_corner_neighbors"          : "fclaw_patch_corner_neighbors",
    "fclaw2d_patch_corner_swap"               : "fclaw_patch_corner_swap",
    "fclaw2d_patch_transform_corner"          : "fclaw_patch_2d_transform_corner",
    "fclaw2d_patch_transform_corner2"         : "fclaw_patch_2d_transform_corner2",

    "fclaw2d_domain_set_refinement"           : "fclaw_domain_set_refinement",
    "fclaw2d_patch_mark_refine"               : "fclaw_patch_mark_refine",
    "fclaw2d_patch_mark_coarsen"              : "fclaw_patch_mark_coarsen",
    "fclaw2d_match_callback_t"                : "fclaw_match_callback_t",
    "fclaw2d_domain_iterate_adapted"          : "fclaw_domain_iterate_adapted",

    "fclaw2d_domain_allocate_before_partition": "fclaw_domain_allocate_before_partition",
    "fclaw2d_domain_retrieve_after_partition" : "fclaw_domain_retrieve_after_partition",
    "fclaw2d_transfer_callback_t"             : "fclaw_transfer_callback_t",
    "fclaw2d_domain_iterate_partitioned"      : "fclaw_domain_iterate_partitioned",
    "fclaw2d_domain_free_after_partition"     : "fclaw_domain_free_after_partition",

    "fclaw2d_domain_exchange"                 : "fclaw_domain_exchange",
    "fclaw2d_domain_exchange_t"               : "fclaw_domain_exchange_t",

    "fclaw2d_domain_allocate_before_exchange" : "fclaw_domain_allocate_before_exchange",
    "fclaw2d_domain_ghost_exchange"           : "fclaw_domain_ghost_exchange",
    "fclaw2d_domain_ghost_exchange_begin"     : "fclaw_domain_ghost_exchange_begin",
    "fclaw2d_domain_ghost_exchange_end"       : "fclaw_domain_ghost_exchange_end", 
    "fclaw2d_domain_free_after_exchange"      : "fclaw_domain_free_after_exchange",
    "fclaw2d_domain_indirect_t"               : "fclaw_domain_indirect_t",
    "fclaw2d_domain_indirect_begin"           : "fclaw_domain_indirect_begin",
    "fclaw2d_domain_indirect_end"             : "fclaw_domain_indirect_end",
    "fclaw2d_domain_indirect_neighbors"       : "fclaw_domain_indirect_neighbors",
    "fclaw2d_domain_indirect_destroy"         : "fclaw_domain_indirect_destroy",
    "fclaw2d_domain_global_maximum"           : "fclaw_domain_global_maximum",
    "fclaw2d_domain_global_sum"               : "fclaw_domain_global_sum",
    "fclaw2d_domain_barrier"                  : "fclaw_domain_barrier",
    "fclaw2d_domain_serialization_enter"      : "fclaw_domain_serialization_enter",
    "fclaw2d_domain_serialization_leave"      : "fclaw_domain_serialization_leave",

    "fclaw2d_domain_is_meta"                  : "fclaw_domain_is_meta",

    #===== ./src/forestclaw2d_convenience.h =====
    "<fclaw2d_convenience.h>"       : "<fclaw_convenience.h>",
    "fclaw2d_domain_new_unitsquare" : "fclaw_domain_new_unitsquare",
    "fclaw2d_domain_new_torus"      : "fclaw_domain_new_2d_torus",
    "fclaw2d_domain_new_twosphere"  : "fclaw_domain_new_2d_twosphere",
    "fclaw2d_domain_new_cubedsphere": "fclaw_domain_new_2d_cubedsphere",
    "fclaw2d_domain_new_disk"       : "fclaw_domain_new_2d_disk",
    "fclaw2d_domain_new_brick"      : "fclaw_domain_new_2d_brick",
    "fclaw2d_domain_destroy"        : "fclaw_domain_destroy",

    "fclaw2d_domain_adapt"               : "fclaw_domain_adapt",
    "fclaw2d_domain_partition"           : "fclaw_domain_partition",
    "fclaw2d_domain_partition_unchanged" : "fclaw_domain_partition_unchanged",
    "fclaw2d_domain_complete"            : "fclaw_domain_complete",
    "fclaw2d_domain_write_vtk"           : "fclaw_domain_write_vtk",
    "fclaw2d_domain_list_levels"         : "fclaw_domain_list_levels",
    "fclaw2d_domain_list_neighbors"      : "fclaw_domain_list_neighbors",
    "fclaw2d_domain_list_adapted"        : "fclaw_domain_list_adapted",
    "fclaw2d_domain_search_points"       : "fclaw_domain_search_points",

    "fclaw2d_integrate_ray_t"            : "fclaw_integrate_ray_t",
    "fclaw2d_domain_integrate_rays"      : "fclaw_domain_integrate_rays",
    "fclaw2d_interpolate_point_t"        : "fclaw_interpolate_point_t",
    "fclaw2d_overlap_exchange"           : "fclaw_overlap_exchange",

    #===== ./src/fclaw2d_advance.h =====
    "<fclaw2d_advance.h>"        : "<fclaw_advance.h>",
    "fclaw2d_advance_all_levels" : "fclaw_advance_all_levels",

    #===== ./src/fclaw2d_block.h =====
    "<fclaw2d_block.h>"                : "<fclaw_block.h>",
    "fclaw2d_block_get_block_boundary" : "fclaw_block_get_block_boundary",

    #===== ./src/fclaw2d_corner_neighbors.h =====
    "fclaw2d_corner_neighbors.h" : "fclaw_corner_neighbors.h",

    #===== ./src/fclaw2d_diagnostics.h =====
    "<fclaw2d_diagnostics.h>"               : "<fclaw_diagnostics.h>",
    "fclaw2d_diagnostics_vtable"            : "fclaw_diagnostics_vtable",
    "fclaw2d_diagnostics_vtable_t"          : "fclaw_diagnostics_vtable_t",
    "fclaw2d_diagnostics_accumulator"       : "fclaw_diagnostics_accumulator",
    "fclaw2d_diagnostics_accumulator_t"     : "fclaw_diagnostics_accumulator_t",
    "fclaw2d_diagnostics_initialize_t"      : "fclaw_diagnostics_initialize_t",
    "fclaw2d_diagnostics_compute_t"         : "fclaw_diagnostics_compute_t",
    "fclaw2d_diagnostics_gather_t"          : "fclaw_diagnostics_gather_t",
    "fclaw2d_diagnostics_reset_t"           : "fclaw_diagnostics_reset_t",
    "fclaw2d_diagnostics_finalize_t"        : "fclaw_diagnostics_finalize_t",
    "fclaw2d_diagnostics_vt"                : "fclaw_diagnostics_vt",
    "fclaw2d_diagnostics_vtable_initialize" : "fclaw_diagnostics_vtable_initialize",
    "fclaw2d_domain_global_minimum"         : "fclaw_domain_global_minimum",
    "fclaw2d_diagnostics_initialize"        : "fclaw_diagnostics_initialize",
    "fclaw2d_diagnostics_gather"            : "fclaw_diagnostics_gather",
    "fclaw2d_diagnostics_reset"             : "fclaw_diagnostics_reset",
    "fclaw2d_diagnostics_finalize"          : "fclaw_diagnostics_finalize",

    #===== ./src/fclaw2d_domain.h =====
    "<fclaw2d_domain.h>"                   : "<fclaw_domain.h>",
    "fclaw2d_domain_data"                  : "fclaw_domain_data",
    "fclaw2d_domain_data_t"                : "fclaw_domain_data_t",
    "fclaw2d_domain_data_new"              : "you_can_safely_remove_this_call",
    "fclaw2d_domain_data_delete"           : "you_can_safely_remove_this_call",
    "fclaw2d_domain_setup"                 : "fclaw_domain_setup",
    "fclaw2d_domain_reset"                 : "fclaw_domain_reset",
    "fclaw2d_domain_get_data"              : "fclaw_domain_get_data",
    "fclaw2d_domain_iterate_level_mthread" : "fclaw_domain_iterate_level_mthread",

    #===== ./src/fclaw2d_elliptic_solver.h =====
    "<fclaw2d_elliptic_solver.h>"        : "<fclaw_elliptic_solver.h>",
    "fclaw2d_elliptic_setup_t"           : "fclaw_elliptic_setup_t",
    "fclaw2d_elliptic_rhs_t"             : "fclaw_elliptic_rhs_t",
    "fclaw2d_elliptic_solve_t"           : "fclaw_elliptic_solve_t",
    "fclaw2d_elliptic_physical_bc_t"     : "fclaw_elliptic_physical_bc_t",
    "fclaw2d_elliptic_vtable"            : "fclaw_elliptic_vtable",
    "fclaw2d_elliptic_vtable_t"          : "fclaw_elliptic_vtable_t",
    "fclaw2d_elliptic_vtable_initialize" : "fclaw_elliptic_vtable_initialize",
    "fclaw2d_elliptic_solve"             : "fclaw_elliptic_solve",
    "fclaw2d_elliptic_vt"                : "fclaw_elliptic_vt",

    #===== ./src/fclaw2d_exchange.h =====
    "<fclaw2d_exchange.h>"                 : "<fclaw_exchange.h>",
    "fclaw2d_exchange_setup"               : "fclaw_exchange_setup",
    "fclaw2d_exchange_delete"              : "fclaw_exchange_delete",
    "fclaw2d_exchange_ghost_patches_begin" : "fclaw_exchange_ghost_patches_begin",
    "fclaw2d_exchange_ghost_patches_end"   : "fclaw_exchange_ghost_patches_end",

    #===== ./src/fclaw2d_face_neighbors.h =====
    "<fclaw2d_face_neighbors.h>"  : "<fclaw_face_neighbors.h>",
    "fclaw2d_face_neighbor_ghost" : "fclaw_face_neighbor_ghost",

    #===== ./src/fclaw2d_farraybox.hpp =====
    "<fclaw2d_farraybox.hpp>"      : "<fclaw_farraybox.hpp>",
    "fclaw2d_farraybox_set_to_nan" : "fclaw_farraybox_set_to_nan",

    #===== ./src/fclaw2d_forestclaw.h =====
    "<fclaw2d_forestclaw.h>"     : "<fclaw_forestclaw.h>",
    "fclaw2d_problem_setup"      : "fclaw_problem_setup",
    "fclaw2d_vtables_initialize" : "fclaw_vtables_initialize",
    "fclaw2d_initialize"         : "fclaw_initialize",
    "fclaw2d_run"                : "fclaw_run",
    "fclaw2d_finalize"           : "fclaw_finalize",
    
    #===== ./src/fclaw2d_ghost_fill.h =====
    "<fclaw2d_ghost_fill.h>"             : "<fclaw_ghost_fill.h>",
    "fclaw2d_ghost_fill_parallel_mode"   : "fclaw_ghost_fill_parallel_mode",
    "FCLAW2D_BOUNDARY_INTERIOR_ONLY"     : "FCLAW_BOUNDARY_INTERIOR_ONLY",
    "FCLAW2D_BOUNDARY_GHOST_ONLY"        : "FCLAW_BOUNDARY_GHOST_ONLY",
    "FCLAW2D_BOUNDARY_ALL"               : "FCLAW_BOUNDARY_ALL",
    "FCLAW2D_BOUNDARY_LOCAL_ALL"         : "FCLAW_BOUNDARY_LOCAL_ALL",
    "fclaw2d_ghost_fill_parallel_mode_t" : "fclaw_ghost_fill_parallel_mode_t",
    "fclaw2d_exchange_type"              : "fclaw_exchange_type",
    "FCLAW2D_COPY"                       : "FCLAW_COPY",
    "FCLAW2D_AVERAGE"                    : "FCLAW_AVERAGE",
    "FCLAW2D_INTERPOLATE"                : "FCLAW_INTERPOLATE",
    "FCLAW2D_TIME_SYNC_F2C"              : "FCLAW_TIME_SYNC_F2C",
    "FCLAW2D_TIME_SYNC_SAMESIZE"         : "FCLAW_TIME_SYNC_SAMESIZE",
    "fclaw2d_exchange_type_t"            : "fclaw_exchange_type_t",
    "fclaw2d_grid_type"                  : "fclaw_grid_type",
    "FCLAW2D_IS_COARSE"                  : "FCLAW_IS_COARSE",
    "FCLAW2D_IS_FINE"                    : "FCLAW_IS_FINE",
    "fclaw2d_grid_type_t"                : "fclaw_grid_type_t",
    "fclaw2d_exchange_info"              : "fclaw_exchange_info",
    "fclaw2d_exchange_info_t"            : "fclaw_exchange_info_t",
    "fclaw2d_ghost_update"               : "fclaw_ghost_update",
    "fclaw2d_ghost_update_async"         : "fclaw_ghost_update_async",
    "fclaw2d_ghost_update_nonasync"      : "fclaw_ghost_update_nonasync",
    "fclaw2d_face_neighbor_ghost"        : "fclaw_face_neighbor_ghost",

    #===== ./src/fclaw2d_global.h =====
    "<fclaw2d_global.h>"                   : "<fclaw_global.h>",
    "fclaw2d_iterate_patch_cb"             : "fclaw_iterate_patch_cb",
    "fclaw2d_iterate_family_cb"            : "fclaw_iterate_family_cb",
    "fclaw2d_global_t"                     : "fclaw_global_t",
    "fclaw2d_global_iterate_t"             : "fclaw_global_iterate_t",
    "fclaw2d_global"                       : "fclaw_global",
    "fclaw2d_global_iterate"               : "fclaw_global_iterate",
    "fclaw2d_global_new"                   : "fclaw_global_new",
    "fclaw2d_global_new_comm"              : "fclaw_global_new_comm",
    "fclaw2d_global_destroy"               : "fclaw_global_destroy",
    "fclaw2d_global_store_domain"          : "fclaw_global_store_domain",
    "fclaw2d_global_store_map"             : "fclaw_map_store",
    "fclaw2d_global_iterate_level"         : "fclaw_global_iterate_level",
    "fclaw2d_global_iterate_patches"       : "fclaw_global_iterate_patches",
    "fclaw2d_global_iterate_families"      : "fclaw_global_iterate_families",
    "fclaw2d_global_iterate_adapted"       : "fclaw_global_iterate_adapted",
    "fclaw2d_global_iterate_level_mthread" : "fclaw_global_iterate_level_mthread",
    "fclaw2d_global_iterate_partitioned"   : "fclaw_global_iterate_partitioned",
    "fclaw2d_global_options_store"         : "fclaw_global_options_store",
    "fclaw2d_global_get_options"           : "fclaw_global_get_options",
    "fclaw2d_global_set_global"            : "fclaw_global_set_static",
    "fclaw2d_global_unset_global"          : "fclaw_global_clear_static",
    "fclaw2d_global_get_global"            : "fclaw_global_get_static_global",
    "fclaw2d_set_global_context"           : "fclaw_set_global_context",
    "fclaw2d_clear_global_context"         : "fclaw_clear_global_context",

    #===== ./src/fclaw2d_include_all.h =====
    "<fclaw2d_include_all.h>" : "<fclaw_include_all.h>",

    #===== ./src/fclaw2d_options.h =====
    "<fclaw2d_options.h>"   : "<fclaw_options.h>",
    "fclaw2d_options_store" : "fclaw_options_store",
    "fclaw2d_get_options"   : "fclaw_get_options",

    #===== ./src/fclaw2d_output.h =====
    "<fclaw2d_output.h>"              : "<fclaw_output.h>",
    "fclaw2d_output_frame"            : "fclaw_output_frame",

    #===== ./src/fclaw2d_partition.h =====
    "<fclaw2d_partition.h>"    : "<fclaw_partition.h>",
    "fclaw2d_partition_domain" : "fclaw_partition_domain",

    #===== ./src/fclaw2d_patch.h =====
    "<fclaw2d_patch.h>"                     : "<fclaw_patch.h>",

    "fclaw2d_patch_vtable_t"                : "fclaw_patch_vtable_t",
    "fclaw2d_patch_data_t"                  : "fclaw_patch_data_t",
    "fclaw2d_patch_transform_data_t"        : "fclaw_patch_transform_data_t",
    "FCLAW2D_BUILD_FOR_GHOST_AREA_COMPUTED" : "FCLAW_BUILD_FOR_GHOST_AREA_COMPUTED",
    "FCLAW2D_BUILD_FOR_GHOST_AREA_PACKED"   : "FCLAW_BUILD_FOR_GHOST_AREA_PACKED",
    "FCLAW2D_BUILD_FOR_UPDATE"              : "FCLAW_BUILD_FOR_UPDATE",
    "FCLAW2D_BUILD_CUSTOM"                  : "FCLAW_BUILD_CUSTOM",
    "fclaw2d_build_mode_t"                  : "fclaw_build_mode_t",
    "fclaw2d_patch_data"                    : "fclaw_patch_data",
    "fclaw2d_patch_transform_data"          : "fclaw_patch_transform_data",

    "fclaw2d_patch_reset_data"                : "fclaw_patch_reset_data",
    "fclaw2d_patch_data_delete"               : "fclaw_patch_data_delete",
    "fclaw2d_patch_build"                     : "fclaw_patch_build",
    "fclaw2d_patch_build_from_fine"           : "fclaw_patch_build_from_fine",
    "fclaw2d_patch_initialize"                : "fclaw_patch_initialize",
    "fclaw2d_patch_physical_bc"               : "fclaw_patch_physical_bc",
    "fclaw2d_patch_single_step_update"        : "fclaw_patch_single_step_update",
    "fclaw2d_patch_set_rhs"                   : "fclaw_patch_set_rhs",
    "fclaw2d_patch_restore_step"              : "fclaw_patch_restore_step",
    "fclaw2d_patch_save_step"                 : "fclaw_patch_save_step",
    "fclaw2d_patch_setup_timeinterp"          : "fclaw_patch_setup_timeinterp",
    "fclaw2d_patch_copy_face"                 : "fclaw_patch_copy_face",
    "fclaw2d_patch_average_face"              : "fclaw_patch_average_face",
    "fclaw2d_patch_interpolate_face"          : "fclaw_patch_interpolate_face",
    "fclaw2d_patch_copy_corner"               : "fclaw_patch_copy_corner",
    "fclaw2d_patch_average_corner"            : "fclaw_patch_average_corner",
    "fclaw2d_patch_interpolate_corner"        : "fclaw_patch_interpolate_corner",
    "fclaw2d_patch_create_user_data"          : "fclaw_patch_create_user_data",
    "fclaw2d_patch_destroy_user_data"         : "fclaw_patch_destroy_user_data",
    "fclaw2d_patch_transform_init_data"       : "fclaw_patch_transform_init_data",
    "fclaw2d_patch_transform_blockface"       : "fclaw_patch_transform_blockface",
    "fclaw2d_patch_transform_blockface_intra" : "fclaw_patch_transform_blockface_intra",
    "fclaw2d_patch_tag4refinement"            : "fclaw_patch_tag4refinement",
    "fclaw2d_patch_tag4coarsening"            : "fclaw_patch_tag4coarsening",
    "fclaw2d_patch_interpolate2fine"          : "fclaw_patch_interpolate2fine",
    "fclaw2d_patch_average2coarse"            : "fclaw_patch_average2coarse",
    "fclaw2d_patch_ghost_packsize"            : "fclaw_patch_ghost_packsize",
    "fclaw2d_patch_local_ghost_alloc"         : "fclaw_patch_local_ghost_alloc",
    "fclaw2d_patch_local_ghost_free"          : "fclaw_patch_local_ghost_free",
    "fclaw2d_patch_local_ghost_pack"          : "fclaw_patch_local_ghost_pack",
    "fclaw2d_patch_remote_ghost_build"        : "fclaw_patch_remote_ghost_build",
    "fclaw2d_patch_remote_ghost_unpack"       : "fclaw_patch_remote_ghost_unpack",
    "fclaw2d_patch_remote_ghost_delete"       : "fclaw_patch_remote_ghost_delete",
    "fclaw2d_patch_partition_pack"            : "fclaw_patch_partition_pack",
    "fclaw2d_patch_partition_unpack"          : "fclaw_patch_partition_unpack",
    "fclaw2d_patch_partition_packsize"        : "fclaw_patch_partition_packsize",
    "fclaw2d_patch_time_sync_f2c"             : "fclaw_patch_time_sync_f2c",
    "fclaw2d_patch_time_sync_samesize"        : "fclaw_patch_time_sync_samesize",
    "fclaw2d_patch_time_sync_reset"           : "fclaw_patch_time_sync_reset",

    "fclaw2d_patch_new_t"                     : "fclaw_patch_new_t",
    "fclaw2d_patch_delete_t"                  : "fclaw_patch_delete_t",
    "fclaw2d_patch_build_t"                   : "fclaw_patch_build_t",
    "fclaw2d_patch_build_from_fine_t"         : "fclaw_patch_build_from_fine_t",
    "fclaw2d_patch_setup_t"                   : "fclaw_patch_setup_t",
    "fclaw2d_patch_initialize_t"              : "fclaw_patch_initialize_t",
    "fclaw2d_patch_physical_bc_t"             : "fclaw_patch_physical_bc_t",
    "fclaw2d_patch_single_step_update_t"      : "fclaw_patch_single_step_update_t",
    "fclaw2d_patch_rhs_t"                     : "fclaw_patch_rhs_t",
    "fclaw2d_patch_setup_timeinterp_t"        : "fclaw_patch_setup_timeinterp_t",
    "fclaw2d_patch_restore_step_t"            : "fclaw_patch_restore_step_t",
    "fclaw2d_patch_save_step_t"               : "fclaw_patch_save_step_t",
    "fclaw2d_patch_copy_face_t"               : "fclaw_patch_copy_face_t",
    "fclaw2d_patch_average_face_t"            : "fclaw_patch_average_face_t",
    "fclaw2d_patch_interpolate_face_t"        : "fclaw_patch_interpolate_face_t",
    "fclaw2d_patch_copy_corner_t"             : "fclaw_patch_copy_corner_t",
    "fclaw2d_patch_average_corner_t"          : "fclaw_patch_average_corner_t",
    "fclaw2d_patch_interpolate_corner_t"      : "fclaw_patch_interpolate_corner_t",
    "fclaw2d_patch_transform_init_data_t"     : "fclaw_patch_transform_init_data_t",
    "fclaw2d_patch_transform_blockface_t"     : "fclaw_patch_transform_blockface_t",
    "fclaw2d_patch_transform_blockface_intra_t" : "fclaw_patch_transform_blockface_intra_t",
    "fclaw2d_patch_tag4refinement_t"          : "fclaw_patch_tag4refinement_t",
    "fclaw2d_patch_tag4coarsening_t"          : "fclaw_patch_tag4coarsening_t",
    "fclaw2d_patch_interpolate2fine_t"        : "fclaw_patch_interpolate2fine_t",
    "fclaw2d_patch_average2coarse_t"          : "fclaw_patch_average2coarse_t",
    "fclaw2d_patch_ghost_packsize_t"          : "fclaw_patch_ghost_packsize_t",
    "fclaw2d_patch_local_ghost_pack_t"        : "fclaw_patch_local_ghost_pack_t",
    "fclaw2d_patch_local_ghost_alloc_t"       : "fclaw_patch_local_ghost_alloc_t",
    "fclaw2d_patch_local_ghost_free_t"        : "fclaw_patch_local_ghost_free_t",
    "fclaw2d_patch_remote_ghost_build_t"      : "fclaw_patch_remote_ghost_build_t",
    "fclaw2d_patch_remote_ghost_setup_t"      : "fclaw_patch_remote_ghost_setup_t",
    "fclaw2d_patch_remote_ghost_unpack_t"     : "fclaw_patch_remote_ghost_unpack_t",
    "fclaw2d_patch_remote_ghost_delete_t"     : "fclaw_patch_remote_ghost_delete_t",
    "fclaw2d_patch_partition_packsize_t"      : "fclaw_patch_partition_packsize_t",
    "fclaw2d_patch_partition_pack_t"          : "fclaw_patch_partition_pack_t",
    "fclaw2d_patch_partition_unpack_t"        : "fclaw_patch_partition_unpack_t",
    "fclaw2d_patch_time_sync_f2c_t"           : "fclaw_patch_time_sync_f2c_t",
    "fclaw2d_patch_time_sync_samesize_t"      : "fclaw_patch_time_sync_samesize_t",
    "fclaw2d_patch_time_sync_reset_t"         : "fclaw_patch_time_sync_reset_t",
    "fclaw2d_patch_create_user_data_t"        : "fclaw_patch_create_user_data_t",
    "fclaw2d_patch_destroy_user_data_t"       : "fclaw_patch_destroy_user_data_t",
    "fclaw2d_patch_metric_patch_t"            : "fclaw_patch_metric_patch_t",

    "fclaw2d_patch_vtable"                    : "fclaw_patch_vtable",
    "fclaw2d_patch_vt"                        : "fclaw_patch_vt",
    "fclaw2d_patch_vtable_initialize"         : "fclaw_patch_vtable_initialize",

    "fclaw2d_patch_get_info"                  : "fclaw_patch_get_info",
    "fclaw2d_patch_get_info2"                 : "fclaw_patch_get_info2",
    "fclaw2d_patch_get_user_patch"            : "fclaw_patch_get_user_patch",
    "fclaw2d_patch_get_patch_data"            : "fclaw_patch_get_patch_data",
    "fclaw2d_patch_get_user_data"             : "fclaw_patch_get_user_data",
    "fclaw2d_patch_metric_patch"              : "fclaw_patch_metric_patch",
    "fclaw2d_patch_get_blockno"               : "fclaw_patch_get_blockno",
    "fclaw2d_patch_get_patchno"               : "fclaw_patch_get_patchno",
    "fclaw2d_patch_user_data"                 : "fclaw_patch_user_data",
    "fclaw2d_patch_set_user_data"             : "fclaw_patch_set_user_data",

    "fclaw2d_patch_on_parallel_boundary"      : "fclaw_patch_on_parallel_boundary",
    "fclaw2d_patch_set_face_type"             : "fclaw_patch_set_face_type",
    "fclaw2d_patch_set_corner_type"           : "fclaw_patch_set_corner_type",
    "fclaw2d_patch_set_missing_corner"        : "fclaw_patch_set_missing_corner",
    "fclaw2d_patch_get_face_type"             : "fclaw_patch_get_face_type",
    "fclaw2d_patch_get_corner_type"           : "fclaw_patch_get_corner_type",
    "fclaw2d_patch_corner_is_missing"         : "fclaw_patch_corner_is_missing",
    "fclaw2d_patch_neighbors_set"             : "fclaw_patch_neighbors_set",
    "fclaw2d_patch_neighbors_reset"           : "fclaw_patch_neighbors_reset",
    "fclaw2d_patch_neighbor_type_set"         : "fclaw_patch_neighbor_type_set",
    "fclaw2d_patch_has_finegrid_neighbors"    : "fclaw_patch_has_finegrid_neighbors",
    "fclaw2d_patch_on_coarsefine_interface"   : "fclaw_patch_on_coarsefine_interface",
    "fclaw2d_patch_block_corner_count"        : "fclaw_patch_block_corner_count",
    "fclaw2d_patch_set_block_corner_count"    : "fclaw_patch_set_block_corner_count",

    #===== ./src/fclaw2d_physical_bc.h =====
    "<fclaw2d_physical_bc.h>"      : "<fclaw_physical_bc.h>",

    "fclaw2d_physical_time_info"   : "fclaw_physical_time_info",
    "fclaw2d_physical_time_info_t" : "fclaw_physical_time_info_t",
    "cb_fclaw2d_physical_set_bc"   : "cb_fclaw_physical_set_bc",
    "fclaw2d_physical_get_bc"      : "fclaw_physical_get_bc",
    "fclaw2d_physical_set_bc"      : "fclaw_physical_set_bc",
    "fclaw2d_physical_bc_default"  : "fclaw_physical_bc_default",

    #===== ./src/fclaw2d_regrid.h =====
    "<fclaw2d_regrid.h>"                     : "<fclaw_regrid.h>",
    "cb_fclaw2d_regrid_tag4refinement"       : "cb_fclaw_regrid_tag4refinement",
    "cb_fclaw2d_regrid_repopulate"           : "cb_fclaw_regrid_repopulate",
    "fclaw2d_regrid_set_neighbor_types"      : "fclaw_regrid_set_neighbor_types",
    "fclaw2d_regrid"                         : "fclaw_regrid",
    "fclaw2d_after_regrid"                   : "fclaw_after_regrid",
    "fclaw2d_regrid_set_neighbor_types"      : "fclaw_regrid_set_neighbor_types",

    #===== ./src/fclaw2d_time_sync.h =====
    "<fclaw2d_time_sync.h>"                  : "<fclaw_time_sync.h>",
    "fclaw2d_time_sync_type"                 : "fclaw_time_sync_type",
    "FCLAW2D_TIME_SYNC_RESET_F2C"            : "FCLAW_TIME_SYNC_RESET_F2C",
    "FCLAW2D_TIME_SYNC_RESET_SAMESIZE"       : "FCLAW_TIME_SYNC_RESET_SAMESIZE",
    "FCLAW2D_TIME_SYNC_RESET_PHYS"           : "FCLAW_TIME_SYNC_RESET_PHYS",             
    "fclaw2d_time_sync_type_t"               : "fclaw_time_sync_type_t",
    "fclaw2d_time_sync"                      : "fclaw_time_sync",

    #===== ./src/fclaw2d_timeinterp.h =====
    "<fclaw2d_timeinterp.h>"                 : "<fclaw_timeinterp.h>",
    "fclaw2d_timeinterp"                     : "fclaw_timeinterp",

    #===== ./src/fclaw_timer.h =====
    "FCLAW2D_TIMER_NONE"                     : "FCLAW_TIMER_NONE",
    "FCLAW2D_TIMER_INIT"                     : "FCLAW_TIMER_INIT",
    "FCLAW2D_TIMER_ADVANCE"                  : "FCLAW_TIMER_ADVANCE",
    "FCLAW2D_TIMER_ELLIPTIC_SOLVE"           : "FCLAW_TIMER_ELLIPTIC_SOLVE",
    "FCLAW2D_TIMER_GHOSTFILL"                : "FCLAW_TIMER_GHOSTFILL",
    "FCLAW2D_TIMER_REGRID"                   : "FCLAW_TIMER_REGRID",
    "FCLAW2D_TIMER_DIAGNOSTICS"              : "FCLAW_TIMER_DIAGNOSTICS",
    "FCLAW2D_TIMER_OUTPUT"                   : "FCLAW_TIMER_OUTPUT",
    "FCLAW2D_TIMER_GHOSTPATCH_COMM"          : "FCLAW_TIMER_GHOSTPATCH_COMM",
    "FCLAW2D_TIMER_ADAPT_COMM"               : "FCLAW_TIMER_ADAPT_COMM",
    "FCLAW2D_TIMER_PARTITION_COMM"           : "FCLAW_TIMER_PARTITION_COMM",
    "FCLAW2D_TIMER_DIAGNOSTICS_COMM"         : "FCLAW_TIMER_DIAGNOSTICS_COMM",
    "FCLAW2D_TIMER_CFL_COMM"                 : "FCLAW_TIMER_CFL_COMM",
    "FCLAW2D_TIMER_WALLTIME"                 : "FCLAW_TIMER_WALLTIME",
    "FCLAW2D_TIMER_UNACCOUNTED"              : "FCLAW_TIMER_UNACCOUNTED",
    "FCLAW2D_TIMER_ADVANCE_STEPS_COUNTER"    : "FCLAW_TIMER_ADVANCE_STEPS_COUNTER",
    "FCLAW2D_TIMER_ELLIPTIC_GRIDS_COUNTER"   : "FCLAW_TIMER_ELLIPTIC_GRIDS_COUNTER",
    "FCLAW2D_TIMER_GRIDS_PER_PROC"           : "FCLAW_TIMER_GRIDS_PER_PROC",
    "FCLAW2D_TIMER_GRIDS_INTERIOR"           : "FCLAW_TIMER_GRIDS_INTERIOR",
    "FCLAW2D_TIMER_GRIDS_LOCAL_BOUNDARY"     : "FCLAW_TIMER_GRIDS_LOCAL_BOUNDARY",
    "FCLAW2D_TIMER_GRIDS_REMOTE_BOUNDARY"    : "FCLAW_TIMER_GRIDS_REMOTE_BOUNDARY",
    "FCLAW2D_TIMER_REGRID_BUILD"             : "FCLAW_TIMER_REGRID_BUILD",
    "FCLAW2D_TIMER_REGRID_TAGGING"           : "FCLAW_TIMER_REGRID_TAGGING",
    "FCLAW2D_TIMER_TIMESYNC"                 : "FCLAW_TIMER_TIMESYNC",
    "FCLAW2D_TIMER_GHOSTPATCH_BUILD"         : "FCLAW_TIMER_GHOSTPATCH_BUILD",
    "FCLAW2D_TIMER_PARTITION"                : "FCLAW_TIMER_PARTITION",
    "FCLAW2D_TIMER_PARTITION_BUILD"          : "FCLAW_TIMER_PARTITION_BUILD",
    "FCLAW2D_TIMER_ADVANCE_STEP2"            : "FCLAW_TIMER_ADVANCE_STEP2",
    "FCLAW2D_TIMER_ADVANCE_B4STEP2"          : "FCLAW_TIMER_ADVANCE_B4STEP2",
    "FCLAW2D_TIMER_GHOSTFILL_COPY"           : "FCLAW_TIMER_GHOSTFILL_COPY",
    "FCLAW2D_TIMER_GHOSTFILL_AVERAGE"        : "FCLAW_TIMER_GHOSTFILL_AVERAGE",
    "FCLAW2D_TIMER_GHOSTFILL_INTERP"         : "FCLAW_TIMER_GHOSTFILL_INTERP",
    "FCLAW2D_TIMER_GHOSTFILL_PHYSBC"         : "FCLAW_TIMER_GHOSTFILL_PHYSBC",
    "FCLAW2D_TIMER_GHOSTFILL_STEP1"          : "FCLAW_TIMER_GHOSTFILL_STEP1",
    "FCLAW2D_TIMER_GHOSTFILL_STEP2"          : "FCLAW_TIMER_GHOSTFILL_STEP2",
    "FCLAW2D_TIMER_GHOSTFILL_STEP3"          : "FCLAW_TIMER_GHOSTFILL_STEP3",
    "FCLAW2D_TIMER_NEIGHBOR_SEARCH"          : "FCLAW_TIMER_NEIGHBOR_SEARCH",
    "FCLAW2D_TIMER_LOCAL_COMM"               : "FCLAW_TIMER_LOCAL_COMM",
    "FCLAW2D_TIMER_GLOBAL_COMM"              : "FCLAW_TIMER_GLOBAL_COMM",
    "FCLAW2D_TIMER_CUDA_ALLOCATE"            : "FCLAW_TIMER_CUDA_ALLOCATE",
    "FCLAW2D_TIMER_CUDA_MEMCOPY_H2H"         : "FCLAW_TIMER_CUDA_MEMCOPY_H2H",
    "FCLAW2D_TIMER_CUDA_MEMCOPY_H2D"         : "FCLAW_TIMER_CUDA_MEMCOPY_H2D",
    "FCLAW2D_TIMER_CUDA_MEMCOPY_D2H"         : "FCLAW_TIMER_CUDA_MEMCOPY_D2H",
    "FCLAW2D_TIMER_EXTRA1"                   : "FCLAW_TIMER_EXTRA1",
    "FCLAW2D_TIMER_EXTRA2"                   : "FCLAW_TIMER_EXTRA2",
    "FCLAW2D_TIMER_EXTRA3"                   : "FCLAW_TIMER_EXTRA3",
    "FCLAW2D_TIMER_EXTRA4"                   : "FCLAW_TIMER_EXTRA4",
    "FCLAW2D_TIMER_COUNT"                    : "FCLAW_TIMER_COUNT",
    "fclaw2d_timer_names_t"                  : "fclaw_timer_names_t",
    "fclaw2d_timer_t"                        : "fclaw_timer_t",
    "fclaw2d_timer_wtime"                    : "fclaw_timer_wtime",
    "fclaw2d_timer_init"                     : "fclaw_timer_init",
    "fclaw2d_timer_start"                    : "fclaw_timer_start",
    "fclaw2d_timer_stop"                     : "fclaw_timer_stop",
    "fclaw2d_timer_start_threadsafe"         : "fclaw_timer_start_threadsafe",
    "fclaw2d_timer_stop_threadsafe"          : "fclaw_timer_stop_threadsafe",

    #===== ./src/fclaw2d_update_single_step.h =====
    "<fclaw2d_update_single_step.h>"         : "<fclaw_update_single_step.h>",
    "fclaw2d_single_step_buffer_data"        : "fclaw_single_step_buffer_data",
    "fclaw2d_single_step_buffer_data_t"      : "fclaw_single_step_buffer_data_t",
    "fclaw2d_single_step_data"               : "fclaw_single_step_data",
    "fclaw2d_single_step_data_t"             : "fclaw_single_step_data_t",
    "fclaw2d_update_single_step"             : "fclaw_update_single_step",

    #===== ./src/fclaw2d_vtable.h =====
    "<fclaw2d_vtable.h>"                     : "<fclaw_vtable.h>",
    "fclaw2d_problem_setup_t"                : "fclaw_problem_setup_t",
    "fclaw2d_output_frame_t"                 : "fclaw_output_frame_t",
    "fclaw2d_after_regrid_t"                 : "fclaw_after_regrid_t",
    "fclaw2d_vtable"                         : "fclaw_vtable",
    "fclaw2d_vtable_t"                       : "fclaw_vtable_t",
    "fclaw2d_vt"                             : "fclaw_vt",
    "fclaw2d_vtable_initialize"              : "fclaw_vtable_initialize",
    "fclaw2d_after_regrid"                   : "fclaw_after_regrid",
                                            
    #===== ./src/patches/clawpatch/fclaw2d_clawpatch.h =====
    "<fclaw2d_clawpatch.h>"                  : "<fclaw_clawpatch.h>",
    "fclaw2d_clawpatch_vtable_t"             : "fclaw_clawpatch_vtable_t",
    "fclaw2d_clawpatch_vtable_initialize"    : "fclaw_clawpatch_vtable_initialize",
    "fclaw2d_clawpatch_vt"                   : "fclaw_clawpatch_vt",
    "fclaw2d_clawpatch_vtable"               : "fclaw_clawpatch_vtable",
    "fclaw2d_clawpatch_save_current_step"    : "fclaw_clawpatch_save_current_step",
    "fclaw2d_clawpatch_grid_data"            : "fclaw_clawpatch_2d_grid_data",
    "fclaw2d_clawpatch_metric_scalar"        : "fclaw_clawpatch_2d_metric_scalar",
    "fclaw2d_clawpatch_metric_vector"        : "fclaw_clawpatch_2d_metric_vector",
    "fclaw2d_clawpatch_metric_data"          : "fclaw_clawpatch_2d_metric_data",
    "fclaw2d_clawpatch_metric_data2"         : "fclaw_clawpatch_2d_metric_data2",
    "fclaw2d_clawpatch_get_area"             : "fclaw_clawpatch_get_2d_area",
    "fclaw2d_clawpatch_soln_data"            : "fclaw_clawpatch_soln_data",
    "fclaw2d_clawpatch_aux_data"             : "fclaw_clawpatch_aux_data",
    "fclaw2d_clawpatch_rhs_data"             : "fclaw_clawpatch_rhs_data",
    "fclaw2d_clawpatch_elliptic_error_data"  : "fclaw_clawpatch_elliptic_error_data",
    "fclaw2d_clawpatch_elliptic_soln_data"   : "fclaw_clawpatch_elliptic_soln_data",
    "fclaw2d_clawpatch_get_q"                : "fclaw_clawpatch_get_q",
    "fclaw2d_clawpatch_get_error"            : "fclaw_clawpatch_get_error",
    "fclaw2d_clawpatch_get_exactsoln"        : "fclaw_clawpatch_get_exactsoln",
    "fclaw2d_clawpatch_size"                 : "fclaw_clawpatch_size",
    "fclaw2d_clawpatch_get_user_data"        : "fclaw_clawpatch_get_user_data",
    "fclaw2d_clawpatch_set_user_data"        : "fclaw_clawpatch_set_user_data",
    "fclaw2d_clawpatch_get_solver_data"      : "fclaw_clawpatch_get_solver_data",
    "fclaw2d_clawpatch_set_solver_data"      : "fclaw_clawpatch_set_solver_data",
    "fclaw2d_clawpatch_timesync_data"        : "fclaw_clawpatch_timesync_data",
    "fclaw2d_clawpatch_get_q_timesync"       : "fclaw_clawpatch_get_q_timesync",
    "fclaw2d_clawpatch_get_registers"        : "fclaw_clawpatch_get_2d_registers",

    #===== ./src/patches/clawpatch/fclaw2d_clawpatch.hpp =====
    "<fclaw2d_clawpatch.hpp>"                : "<fclaw_clawpatch.hpp>",
    "fclaw2d_clawpatch_t"                    : "fclaw_clawpatch_t",
    "fclaw2d_clawpatch_get_clawpatch"        : "fclaw_clawpatch_get_clawpatch",
    "fclaw2d_clawpatch_get_metric_patch"     : "fclaw_clawpatch_get_2d_metric_patch",

    #===== ./src/patches/clawpatch/fclaw2d_clawpatch_diagnostics.h =====
    "<fclaw2d_clawpatch_diagnostics.h>"               : "<fclaw_clawpatch_diagnostics.h>",
    "fclaw2d_clawpatch_diagnostics_initialize"        : "fclaw_clawpatch_diagnostics_initialize",
    "fclaw2d_clawpatch_diagnostics_compute"           : "fclaw_clawpatch_diagnostics_compute",
    "fclaw2d_clawpatch_diagnostics_gather"            : "fclaw_clawpatch_diagnostics_gather",
    "fclaw2d_clawpatch_diagnostics_reset"             : "fclaw_clawpatch_diagnostics_reset",
    "fclaw2d_clawpatch_diagnostics_finalize"          : "fclaw_clawpatch_diagnostics_finalize",
    "fclaw2d_clawpatch_diagnostics_vtable_initialize" : "fclaw_clawpatch_diagnostics_vtable_initialize",
    "fclaw2d_clawpatch_diagnostics_cons_default"      : "fclaw_clawpatch_diagnostics_cons_default",
    "fclaw2d_clawpatch_diagnostics_error_default"     : "fclaw_clawpatch_diagnostics_error_default",

    #===== ./src/patches/clawpatch/fclaw2d_clawpatch_options.h =====
    "<fclaw2d_clawpatch_options.h>"           : "<fclaw_clawpatch_options.h>",
    "fclaw2d_clawpatch_options_t"             : "fclaw_clawpatch_options_t",
    "fclaw2d_clawpatch_options"               : "fclaw_clawpatch_options",
    "fclaw2d_clawpatch_options_register"      : "fclaw_clawpatch_2d_options_register",
    "fclaw2d_clawpatch_options_store"         : "fclaw_clawpatch_options_store",
    "fclaw2d_clawpatch_get_options"           : "fclaw_clawpatch_get_options",

    #===== ./src/patches/clawpatch/fclaw2d_clawpatch_output_ascii.h =====
    "<fclaw2d_clawpatch_output_ascii.h>"      : "<fclaw_clawpatch_output_ascii.h>",
    "fclaw2d_clawpatch_output_ascii"          : "fclaw_clawpatch_output_ascii",
    "fclaw2d_clawpatch_time_header_ascii"     : "fclaw_clawpatch_time_header_ascii",

    #===== ./src/patches/clawpatch/fclaw2d_clawpatch_output_vtk.h =====
    "<fclaw2d_clawpatch_output_vtk.h>"        : "<fclaw_clawpatch_output_vtk.h>",
    "fclaw2d_vtk_patch_data_t"                : "fclaw_vtk_patch_data_t",
    "fclaw2d_vtk_write_file"                  : "fclaw_vtk_write_2d_file",
    "fclaw2d_clawpatch_output_vtk"            : "fclaw_clawpatch_output_vtk",

    #===== ./src/patches/clawpatch/fclaw2d_clawpatch_pillow.h =====
    "<fclaw2d_clawpatch_pillow.h>"               : "<fclaw_clawpatch_pillow.h>",
                               
    "fclaw2d_clawpatch_pillow_vtable_t"          : "fclaw_clawpatch_pillow_vtable_t",
    "fclaw2d_clawpatch_use_pillowsphere"         : "fclaw_clawpatch_use_pillowsphere",
    "fclaw2d_clawpatch_pillow_vtable_initialize" : "fclaw_clawpatch_pillow_vtable_initialize",
    "fclaw2d_clawpatch_pillow_vtable"            : "fclaw_clawpatch_pillow_vtable",
    "fclaw2d_clawpatch_pillow_vt"                : "fclaw_clawpatch_pillow_vt",

    #===== ./src/patches/clawpatch/fclaw3dx_clawpatch.h =====
    "<fclaw3dx_clawpatch.h>"                    : "<fclaw_clawpatch.h>",
    "fclaw3dx_clawpatch_vtable_t"               : "fclaw_clawpatch_vtable_t",
    "fclaw3dx_clawpatch_vtable_initialize"      : "fclaw_clawpatch_vtable_initialize",
    "fclaw3dx_clawpatch_vt"                     : "fclaw_clawpatch_vt",
    "fclaw3dx_clawpatch_vtable"                 : "fclaw_clawpatch_vtable",
    "fclaw3dx_clawpatch_save_current_step"      : "fclaw_clawpatch_save_current_step",
    "fclaw3dx_clawpatch_grid_data"              : "fclaw_clawpatch_3d_grid_data",
    "fclaw3d_clawpatch_grid_data"               : "fclaw_clawpatch_3d_grid_data",
    "fclaw3d_clawpatch_get_volume"              : "fclaw_clawpatch_get_3d_volume",
    "fclaw3d_clawpatch_metric_scalar"           : "fclaw_clawpatch_3d_metric_scalar",
    "fclaw3d_clawpatch_metric_basis"            : "fclaw_clawpatch_3d_metric_basis",
    "fclaw3d_clawpatch_mesh_data"               : "fclaw_clawpatch_3d_mesh_data",
    "fclaw3dx_clawpatch_soln_data"              : "fclaw_clawpatch_soln_data",
    "fclaw3dx_clawpatch_aux_data"               : "fclaw_clawpatch_aux_data",
    "fclaw3dx_clawpatch_rhs_data"               : "fclaw_clawpatch_rhs_data",
    "fclaw3dx_clawpatch_elliptic_error_data"    : "fclaw_clawpatch_elliptic_error_data",
    "fclaw3dx_clawpatch_elliptic_soln_data"     : "fclaw_clawpatch_elliptic_soln_data",
    "fclaw3dx_clawpatch_get_q"                  : "fclaw_clawpatch_get_q",
    "fclaw3dx_clawpatch_get_error"              : "fclaw_clawpatch_get_error",
    "fclaw3dx_clawpatch_get_exactsoln"          : "fclaw_clawpatch_get_exactsoln",
    "fclaw3dx_clawpatch_size"                   : "fclaw_clawpatch_size",
    "fclaw3dx_clawpatch_get_user_data"          : "fclaw_clawpatch_get_user_data",
    "fclaw3dx_clawpatch_set_user_data"          : "fclaw_clawpatch_set_user_data",
    "fclaw3dx_clawpatch_get_solver_data"        : "fclaw_clawpatch_get_solver_data",
    "fclaw3dx_clawpatch_set_solver_data"        : "fclaw_clawpatch_set_solver_data",
    "fclaw3dx_clawpatch_timesync_data"          : "fclaw_clawpatch_timesync_data",
    "fclaw3dx_clawpatch_get_q_timesync"         : "fclaw_clawpatch_get_q_timesync",

    #===== ./src/patches/clawpatch/fclaw3dx_clawpatch.hpp =====
    "<fclaw3dx_clawpatch.hpp>"                  : "<fclaw_clawpatch.hpp>",
    "fclaw3dx_clawpatch_t"                      : "fclaw_clawpatch_t",
    "fclaw3dx_clawpatch_get_clawpatch"          : "fclaw_clawpatch_get_clawpatch",
    "fclaw3dx_clawpatch_get_metric_patch"       : "fclaw_clawpatch_get_3d_metric_patch",

    #===== ./src/patches/clawpatch/fclaw3dx_clawpatch_diagnostics.h =====
    "<fclaw3dx_clawpatch_diagnostics.h>"               : "<fclaw_clawpatch_diagnostics.h>",
    "fclaw3dx_clawpatch_error_info_t"                  : "error_info_t",
    "fclaw3dx_clawpatch_diagnostics_initialize"        : "fclaw_clawpatch_diagnostics_initialize",
    "fclaw3dx_clawpatch_diagnostics_compute"           : "fclaw_clawpatch_diagnostics_compute",
    "fclaw3dx_clawpatch_diagnostics_gather"            : "fclaw_clawpatch_diagnostics_gather",
    "fclaw3dx_clawpatch_diagnostics_reset"             : "fclaw_clawpatch_diagnostics_reset",
    "fclaw3dx_clawpatch_diagnostics_finalize"          : "fclaw_clawpatch_diagnostics_finalize",
    "fclaw3dx_clawpatch_diagnostics_vtable_initialize" : "fclaw_clawpatch_diagnostics_vtable_initialize",
    "fclaw3dx_clawpatch_diagnostics_cons_default"      : "fclaw_clawpatch_diagnostics_cons_default",
    "fclaw3dx_clawpatch_diagnostics_error_default"     : "fclaw_clawpatch_diagnostics_error_default",

    #===== ./src/patches/clawpatch/fclaw3dx_clawpatch_options.h =====
    "<fclaw3dx_clawpatch_options.h>"             : "<fclaw_clawpatch_options.h>",                           
    "fclaw3dx_clawpatch_options_t"             : "fclaw_clawpatch_options_t",
    "fclaw3dx_clawpatch_options"               : "fclaw_clawpatch_options",
    "fclaw3dx_clawpatch_options_t"             : "fclaw_clawpatch_options_t",
    "fclaw3dx_clawpatch_options_register"      : "fclaw_clawpatch_3d_options_register",
    "fclaw3dx_clawpatch_options_store"         : "fclaw_clawpatch_options_store",
    "fclaw3dx_clawpatch_get_options"           : "fclaw_clawpatch_get_options",

    #===== ./src/patches/clawpatch/fclaw3dx_clawpatch_output_ascii.h =====
    "<fclaw3dx_clawpatch_output_ascii.h>"      : "<fclaw_clawpatch_output_ascii.h>",
    "fclaw3dx_clawpatch_output_ascii_cb"       : "fclaw_clawpatch_output_ascii_cb",
    "fclaw3dx_clawpatch_output_ascii"          : "fclaw_clawpatch_output_ascii",
    "fclaw3dx_clawpatch_time_header_ascii"     : "fclaw_clawpatch_time_header_ascii",

    #===== ./src/patches/clawpatch/fclaw3dx_clawpatch_output_vtk.h =====
    "<fclaw3dx_clawpatch_output_vtk.h>"        : "<fclaw_clawpatch_output_vtk.h>",
    "fclaw3dx_vtk_patch_data_t"                : "fclaw_vtk_patch_data_t",
    "fclaw3dx_vtk_write_file"                  : "fclaw_vtk_write_3d_file",
    "fclaw3dx_clawpatch_output_vtk"            : "fclaw_clawpatch_output_vtk",

    #===== ./src/patches/clawpatch/fclaw3dx_clawpatch_pillow.h =====
    "<fclaw3dx_clawpatch_pillow.h>"               : "<fclaw_clawpatch_pillow.h>",
    "fclaw3dx_clawpatch_pillow_vtable_t"          : "fclaw_clawpatch_pillow_vtable_t",
    "fclaw3dx_clawpatch_use_pillowsphere"         : "fclaw_clawpatch_use_pillowsphere",
    "fclaw3dx_clawpatch_pillow_vtable_initialize" : "fclaw_clawpatch_pillow_vtable_initialize",
    "fclaw3dx_clawpatch_pillow_vtable"            : "fclaw_clawpatch_pillow_vtable",
    "fclaw3dx_clawpatch_pillow_vt"                : "fclaw_clawpatch_pillow_vt",

    #===== ./src/patches/clawpatch/fclaw3dx_clawpatch_fort.h =====
    "FCLAW3DX_CLAWPATCH_GET_REFINEMENT_CRITERIA" : "FCLAW3D_CLAWPATCH_GET_REFINEMENT_CRITERIA",
    "FCLAW3DX_CLAWPATCH_TAG_CRITERIA"            : "FCLAW3D_CLAWPATCH_TAG_CRITERIA",
    "FCLAW3DX_CLAWPATCH_EXCEEDS_THRESHOLD"       : "FCLAW3D_CLAWPATCH_EXCEEDS_THRESHOLD",
    "FCLAW3DX_CLAWPATCH_VALUE_EXCEEDS_TH"        : "FCLAW3D_CLAWPATCH_VALUE_EXCEEDS_TH",
    "FCLAW3DX_CLAWPATCH_DIFFERENCE_EXCEEDS_TH"   : "FCLAW3D_CLAWPATCH_DIFFERENCE_EXCEEDS_TH",
    "FCLAW3DX_CLAWPATCH_MINMAX_EXCEEDS_TH"       : "FCLAW3D_CLAWPATCH_MINMAX_EXCEEDS_TH",
    "FCLAW3DX_CLAWPATCH_GRADIENT_EXCEEDS_TH"     : "FCLAW3D_CLAWPATCH_GRADIENT_EXCEEDS_TH",
    "FCLAW3DX_USER_EXCEEDS_TH"                   : "FCLAW3D_USER_EXCEEDS_TH",
    "FCLAW3DX_USER_TAG4REFINEMENT"               : "FCLAW3D_USER_TAG4REFINEMENT",
    "FCLAW3DX_USER_TAG4COARSENING"               : "FCLAW3D_USER_TAG4COARSENING",
    "FCLAW3DX_USER_INTERPOLATE2FINE"             : "FCLAW3D_USER_INTERPOLATE2FINE",
    "FCLAW3DX_USER_AVERAGE2COARSE"               : "FCLAW3D_USER_AVERAGE2COARSE",
    "fclaw3dx_clawpatch_get_refinement_criteria" : "fclaw3d_clawpatch_get_refinement_criteria",
    "fclaw3dx_clawpatch_tag_criteria"            : "fclaw3d_clawpatch_tag_criteria",
    "fclaw3dx_clawpatch_exceeds_threshold"       : "fclaw3d_clawpatch_exceeds_threshold",
    "fclaw3dx_clawpatch_value_exceeds_th"        : "fclaw3d_clawpatch_value_exceeds_th",
    "fclaw3dx_clawpatch_difference_exceeds_th"   : "fclaw3d_clawpatch_difference_exceeds_th",
    "fclaw3dx_clawpatch_minmax_exceeds_th"       : "fclaw3d_clawpatch_minmax_exceeds_th",
    "fclaw3dx_clawpatch_gradient_exceeds_th"     : "fclaw3d_clawpatch_gradient_exceeds_th",
    "fclaw3dx_user_exceeds_th"                   : "fclaw3d_user_exceeds_th",
    "fclaw3dx_user_tag4refinement"               : "fclaw3d_user_tag4refinement",
    "fclaw3dx_user_tag4coarsening"               : "fclaw3d_user_tag4coarsening",
    "fclaw3dx_user_interpolate2fine"             : "fclaw3d_user_interpolate2fine",
    "fclaw3dx_user_average2coarse"               : "fclaw3d_user_average2coarse",

    #===== ./src/patches/clawpatch/fclaw3dx_clawpatch_transform.h =====
    "<fclaw3dx_clawpatch_transform.h>"             : "<fclaw2d_clawpatch_transform.h>",
    "fclaw3dx_clawpatch_transform_init_data"       : "fclaw2d_clawpatch_transform_init_data",
    "fclaw3dx_clawpatch_face_transformation"       : "fclaw2d_clawpatch_face_transformation",
    "fclaw3dx_clawpatch_face_transformation_intra" : "fclaw2d_clawpatch_face_transformation_intra",
    "FCLAW3DX_CLAWPATCH_TRANSFORM_FACE"            : "FCLAW2D_CLAWPATCH_TRANSFORM_FACE",
    "FCLAW3DX_CLAWPATCH_TRANSFORM_FACE_HALF"       : "FCLAW2D_CLAWPATCH_TRANSFORM_FACE_HALF",
    "FCLAW3DX_CLAWPATCH_TRANSFORM_CORNER"          : "FCLAW2D_CLAWPATCH_TRANSFORM_CORNER",
    "FCLAW3DX_CLAWPATCH_TRANSFORM_CORNER_HALF"     : "FCLAW2D_CLAWPATCH_TRANSFORM_CORNER_HALF",
    "fclaw3dx_clawpatch_transform_face"            : "fclaw2d_clawpatch_transform_face",
    "fclaw3dx_clawpatch_transform_face_half"       : "fclaw2d_clawpatch_transform_face_half",
    "fclaw3dx_clawpatch_transform_corner"          : "fclaw2d_clawpatch_transform_corner",
    "fclaw3dx_clawpatch_transform_corner_half"     : "fclaw2d_clawpatch_transform_corner_half",

    #===== ./src/patches/clawpatch/fclaw3dx_clawpatch_fort.h =====
    "<fclaw3dx_clawpatch_fort.h>"                  : "<fclaw3d_clawpatch_fort.h>",

    #===== ./src/patches/clawpatch/fclaw3dx_clawpatch46_fort.h =====
    "<fclaw3dx_clawpatch46_fort.h>"                : "<fclaw3d_clawpatch46_fort.h>",

    #===== ./src/fclaw2d_rays.h =====
    "<fclaw2d_rays.h>"                             : "<fclaw_rays.h>",
    "fclaw2d_ray_vtable"                           : "fclaw_ray_vtable",
    "fclaw2d_ray_vtable_t"                         : "fclaw_ray_vtable_t",
    "fclaw2d_ray"                                  : "fclaw_ray",
    "fclaw2d_ray_t"                                : "fclaw_ray_t",
    "fclaw2d_ray_allocate_and_define_t"            : "fclaw_ray_allocate_and_define_t",
    "fclaw2d_ray_deallocate_t"                     : "fclaw_ray_deallocate_t",
    "fclaw2d_ray_create_files_t"                   : "fclaw_ray_create_files_t",
    "fclaw2d_ray_normalize_t"                      : "fclaw_ray_normalize_t",
    "fclaw2d_ray_update_t"                         : "fclaw_ray_update_t",
    "fclaw2d_ray_print_t"                          : "fclaw_ray_print_t",
    "fclaw2d_ray_destroy_buffer_data_t"            : "fclaw_ray_destroy_buffer_data_t",
    "fclaw2d_ray_vtable"                           : "fclaw_ray_vtable",
    "fclaw2d_ray_vtable_t"                         : "fclaw_ray_vtable_t",
    "fclaw2d_ray_allocate_and_define"              : "fclaw_ray_allocate_and_define",
    "fclaw2d_ray_deallocate"                       : "fclaw_ray_deallocate",
    "fclaw2d_ray_set_ray"                          : "fclaw_ray_set_ray",
    "fclaw2d_ray_get_ray"                          : "fclaw_ray_get_ray",
    "fclaw2d_ray_vt"                               : "fclaw_ray_vt",
    "fclaw2d_ray_vtable_initialize"                : "fclaw_ray_vtable_initialize",
    "fclaw2d_ray_allocate_rays"                    : "fclaw_ray_allocate_rays",
    "fclaw2d_ray_deallocate_rays"                  : "fclaw_ray_deallocate_rays",

    #===== ./src/fclaw2d_map.c =====
    "fclaw2d_map_context"                          : "fclaw_map_context",
    "fclaw2d_map_context_t"                        : "fclaw_map_context_t",
    "fclaw2d_map_store_t"                          : "fclaw_map_store_t",
    "fclaw2d_get_map"                              : "fclaw_get_map",
    "fclaw2d_map_query_t"                          : "fclaw_map_query_t",
    "fclaw2d_map_c2m_t"                            : "fclaw_map_2d_c2m_t",
    "fclaw2d_map_c2m_basis_t"                      : "fclaw_map_2d_c2m_basis_t",
    "fclaw3dx_map_c2m_t"                           : "fclaw_map_3d_c2m_t",
    "fclaw3dx_map_c2m_basis_t"                     : "fclaw_map_3d_c2m_basis_t",
    "fclaw2d_map_destroy_t"                        : "fclaw_map_destroy_t",
    "FCLAW2D_MAP_QUERY"                            : "FCLAW_MAP_QUERY",
    "fclaw2d_map_query"                            : "fclaw_map_query",
    "FCLAW2D_MAP_C2M"                              : "FCLAW_MAP_2D_C2M",
    "fclaw2d_map_c2m"                              : "fclaw_map_2d_c2m",
    "FCLAW2D_MAP_C2M_BASIS"                        : "FCLAW_MAP_2D_C2M_BASIS",
    "fclaw2d_map_c2m_basis"                        : "fclaw_map_2d_c2m_basis",
    "FCLAW3D_MAP_C2M"                              : "FCLAW_MAP_3D_C2M",
    "fclaw3d_map_c2m"                              : "fclaw_map_3d_c2m",
    "FCLAW3D_MAP_C2M_BASIS"                        : "FCLAW_MAP_3D_C2M_BASIS",
    "fclaw3d_map_c2m_basis"                        : "fclaw_map_3d_c2m_basis",
    "FCLAW2D_MAP_BRICK2C"                          : "FCLAW_MAP_2D_BRICK2C",
    "fclaw2d_map_brick2c"                          : "fclaw_map_2d_brick2c",
    "fclaw2d_map_destroy"                          : "fclaw_map_destroy"
}

import glob
import argparse
try:
    import pygments
except ModuleNotFoundError:
    print("The 'pygments' module is not installed.")
    print("This is needed for tokenizing the source code.")
    print("You can install it by running 'pip install pygments' or 'pip3 install pygments'.")
    exit(1)  # Exit the script with an error code
import pygments.lexers
from pygments.token import Token

def replace_identifiers_and_includes(filepath, code, identifier_map):
    lexer = pygments.lexers.get_lexer_for_filename(filepath, stripnl=False, ensurenl=False)
    tokens = lexer.get_tokens(code)
    new_code = ''
    
    changes = False

    for ttype, value in tokens:
        if ttype in Token.Comment.PreprocFile:
            new_value = identifier_map.get(value, value)
            changes = changes or (new_value != value)
            new_code += new_value
        elif ttype in Token.Name:
            new_value = identifier_map.get(value, value)
            changes = changes or (new_value != value)
            new_code += new_value
        else:
            new_code += value
            
    return new_code, changes

def process_directory(root_dir, identifier_map):
    # Use glob to find all C++ files recursively.
    for filepath in [f for f in glob.glob(f"{root_dir}/**/*", recursive=True) if f.lower().endswith(('.c', '.h', '.cpp', '.hpp', '.f', '.f90', '.cu'))]:
        # Read the existing code from the file
        with open(filepath, 'r') as f:
            old_code = f.read()
        
        # Replace the identifiers
        new_code, changes = replace_identifiers_and_includes(filepath, old_code, identifier_map)

        if changes:
            print(f"Upadated   {filepath}...")
            # Write the new code back to the file
            with open(filepath, 'w') as f:
                f.write(new_code)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A simple script to demonstrate argparse.")
    parser.add_argument("-d", "--directory", help="directory to recursively update files in", required=True)
    args = parser.parse_args()
    print("Processing files...")
    process_directory(args.directory, fclaw_1to2)
    print("Done!")