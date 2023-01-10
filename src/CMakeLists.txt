# we use auxiliary objects for each code language to avoid compiler option mingling

# -- c library
add_library(forestclaw_c OBJECT 
  fclaw_base.c 
  fclaw_options.c 
  fclaw_filesystem.cpp
  fclaw_gauges.c 
  fclaw_package.c 
  fclaw_packing.c
  fclaw_pointer_map.cpp
  fclaw_math.c
  fclaw_timer.c
  fclaw_mpi.c
  fclaw2d_block.c
  fclaw2d_options.c
  fclaw2d_global.c
  fclaw2d_forestclaw.c
  fclaw2d_vtable.c
  fclaw2d_map.c
  fclaw2d_map_query.c
  fclaw2d_map_brick.c
  fclaw2d_convenience.c
  fclaw2d_exchange.c
  fclaw2d_advance.c
  fclaw2d_elliptic_solver.c
  fclaw2d_physical_bc.c
  fclaw2d_ghost_fill.c
  fclaw2d_output.c
  fclaw2d_run.c
  fclaw2d_diagnostics.c
  fclaw2d_rays.c
  fclaw2d_update_single_step.c
  fclaw2d_domain.c
  fclaw2d_regrid.c
  fclaw2d_patch.c
  fclaw2d_partition.c
  fclaw2d_timeinterp.c
  fclaw2d_initialize.c
  fclaw2d_finalize.c
  fclaw2d_time_sync.c
  fclaw2d_corner_neighbors.c
  fclaw2d_face_neighbors.c
  fclaw2d_farraybox.cpp
  fclaw2d_output_tikz.c
  forestclaw2d.c
  forestclaw3d.c
  fclaw3d_convenience.c
  fclaw3d_domain.c
  fclaw3d_global.c
  fclaw3d_map.c
  fclaw_domain.c
  fclaw_global.c
  fp_exception_glibc_extension.c
  mappings/fclaw2d_map_nomap.c
  mappings/fclaw2d_map_nomap_brick.c
)
target_include_directories(forestclaw_c PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include patches/clawpatch)
target_link_libraries(forestclaw_c PRIVATE P4EST::P4EST SC::SC)
if(mpi)
  target_link_libraries(forestclaw_c PRIVATE MPI::MPI_C)
endif(mpi)

# -- fortran library
add_library(forestclaw_f OBJECT
  mappings/context.f
  mappings/scale_rot_map.f
  mappings/map_query.f
  mappings/cart/mapc2m_cart.f
  mappings/bilinear/mapc2m_bilinear.f
  mappings/annulus/mapc2m_annulus.f
  mappings/identity/mapc2m_identity.f
  mappings/pillowdisk/mapc2m_pillowdisk.f
  mappings/squareddisk/mapc2m_squareddisk.f
  mappings/pillowdisk5/mapc2m_pillowdisk5.f
  mappings/cubedsphere/mapc2m_cubedsphere.f
  mappings/pillowsphere/mapc2m_pillowsphere.f
  mappings/fivepatch/mapc2m_fivepatch.f
  mappings/torus/mapc2m_torus.f
  mappings/torus/mapc2m_twisted_torus.f
  mappings/latlong/mapc2m_latlong.f
  fortran_source2d/fclaw_dopri5.f
  fortran_source2d/cellave2.f
  fortran_source2d/fdisc.f
  fortran_source2d/philim.f
  fortran_source2d/cellave3.f
)

# -- put together in final library

add_library(forestclaw)
set_target_properties(forestclaw PROPERTIES EXPORT_NAME FORESTCLAW)

target_sources(forestclaw PRIVATE $<TARGET_OBJECTS:forestclaw_f> $<TARGET_OBJECTS:forestclaw_c>)

target_link_libraries(forestclaw PRIVATE ZLIB::ZLIB)
target_link_libraries(forestclaw PUBLIC P4EST::P4EST SC::SC)
if(mpi)
  target_link_libraries(forestclaw PUBLIC MPI::MPI_C INTERFACE MPI::MPI_CXX)
endif(mpi)

target_include_directories(forestclaw
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>)

# imported target, for use from FetchContent
add_library(FORESTCLAW::FORESTCLAW INTERFACE IMPORTED GLOBAL)
target_link_libraries(FORESTCLAW::FORESTCLAW INTERFACE forestclaw)

install(FILES
  fclaw_base.h 
	fclaw_timer.h 
	fclaw_package.h 
  fclaw_packing.h
	fclaw_pointer_map.h 
	fclaw_filesystem.h 
	fclaw_options.h 
	fclaw_gauges.h 
	fclaw_mpi.h 
	fclaw_math.h 
	forestclaw2d.h 
	fp_exception_glibc_extension.h 
  fclaw2d_include_all.h 
	fclaw2d_defs.h 
	fclaw2d_options.h 
	fclaw2d_advance.h 
	fclaw2d_elliptic_solver.h 
	fclaw2d_farraybox.hpp 
	fclaw2d_global.h
	fclaw2d_forestclaw.h
	fclaw2d_domain.h
	fclaw2d_block.h
	fclaw2d_patch.h
	fclaw2d_vtable.h
	fclaw2d_output.h
	fclaw2d_time_sync.h
	fclaw2d_update_single_step.h
	fclaw2d_physical_bc.h
	fclaw2d_rays.h
	fclaw2d_ghost_fill.h
	fclaw2d_corner_neighbors.h
	fclaw2d_face_neighbors.h
	fclaw2d_map.h
	fclaw2d_timeinterp.h
	fclaw2d_regrid.h
	fclaw2d_map_brick.h
	fclaw2d_partition.h
	fclaw2d_exchange.h
	fclaw2d_convenience.h
	fclaw2d_map_query.h
	fclaw2d_map_query_defs.h
	fclaw2d_diagnostics.h
  fclaw2d_to_3d.h
  forestclaw3d.h
  fclaw3d_convenience.h
  fclaw3d_defs.h
  fclaw3d_domain.h
  fclaw3d_global.h
  fclaw3d_map.h
  fclaw_patch.h
  fclaw_domain.h
  fclaw_global.h
  DESTINATION include
)

install(TARGETS forestclaw
  EXPORT ${PROJECT_NAME}-targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)

if(BUILD_TESTING)
  add_executable(forestclaw.TEST
      fclaw_gauges.h.TEST.cpp
      fclaw_packing.h.TEST.cpp
      fclaw_pointer_map.h.TEST.cpp
      fclaw2d_elliptic_solver.h.TEST.cpp
      fclaw2d_diagnostics.h.TEST.cpp
      fclaw2d_global.h.TEST.cpp
      fclaw2d_options.h.TEST.cpp
      fclaw2d_patch.h.TEST.cpp
      fclaw2d_vtable.h.TEST.cpp
  )

  target_link_libraries(forestclaw.TEST testutils forestclaw)
  register_unit_tests(forestclaw.TEST)
endif(BUILD_TESTING)
