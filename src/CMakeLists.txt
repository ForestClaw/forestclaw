# we use auxiliary objects for each code language to avoid compiler option mingling

add_library(forestclaw_c OBJECT fclaw_base.c fclaw_options.c fclaw_gauges.c fclaw_package.c fclaw_math.c fclaw_timer.c fclaw_mpi.c fclaw2d_block.c fclaw2d_defs.c fclaw2d_options.c fclaw2d_global.c fclaw2d_forestclaw.c fclaw2d_vtable.c fclaw2d_map.c fclaw2d_map_query.c fclaw2d_map_brick.c fclaw2d_convenience.c fclaw2d_exchange.c fclaw2d_advance.c fclaw2d_elliptic_solver.c fclaw2d_physical_bc.c fclaw2d_ghost_fill.c fclaw2d_output.c fclaw2d_run.c fclaw2d_diagnostics.c fclaw2d_update_single_step.c fclaw2d_domain.c fclaw2d_regrid.c fclaw2d_patch.c fclaw2d_partition.c fclaw2d_timeinterp.c fclaw2d_initialize.c fclaw2d_finalize.c fclaw2d_time_sync.c fclaw2d_corner_neighbors.c fclaw2d_face_neighbors.c fclaw2d_farraybox.cpp fclaw2d_output_tikz.c forestclaw2d.c fp_exception_glibc_extension.c
mappings/fclaw2d_map_nomap.c mappings/fclaw2d_map_nomap_brick.c)
target_include_directories(forestclaw_c PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include patches/clawpatch)
target_link_libraries(forestclaw_c PRIVATE P4EST::P4EST SC::SC)
if(MPI_FOUND)
  target_link_libraries(forestclaw_c PRIVATE MPI::MPI_C)
endif()

add_library(forestclaw_f OBJECT mappings/context.f mappings/scale_rot_map.f mappings/map_query.f mappings/cart/mapc2m_cart.f mappings/bilinear/mapc2m_bilinear.f mappings/annulus/mapc2m_annulus.f mappings/identity/mapc2m_identity.f mappings/pillowdisk/mapc2m_pillowdisk.f mappings/squareddisk/mapc2m_squareddisk.f mappings/pillowdisk5/mapc2m_pillowdisk5.f mappings/cubedsphere/mapc2m_cubedsphere.f mappings/pillowsphere/mapc2m_pillowsphere.f mappings/fivepatch/mapc2m_fivepatch.f mappings/torus/mapc2m_torus.f mappings/torus/mapc2m_twisted_torus.f mappings/latlong/mapc2m_latlong.f
fortran_source2d/cellave2.f fortran_source2d/fdisc.f fortran_source2d/philim.f)
set_target_properties(forestclaw_f PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
target_include_directories(forestclaw_f INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)

target_sources(forestclaw PRIVATE $<TARGET_OBJECTS:forestclaw_f> $<TARGET_OBJECTS:forestclaw_c>)
