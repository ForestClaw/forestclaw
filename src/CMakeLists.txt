# we use auxiliary objects for each code language to avoid compiler option mingling

# -- c library
add_library(forestclaw_c OBJECT 
  fclaw_advance.c
  fclaw_base.c 
  fclaw_block.c
  fclaw_convenience.c
  fclaw_corner_neighbors.c
  fclaw_diagnostics.c
  fclaw_domain.c
  fclaw_edge_neighbors.c
  fclaw_elliptic_solver.c
  fclaw_exchange.c
  fclaw_face_neighbors.c
  fclaw_farraybox.cpp
  fclaw_filesystem.cpp
  fclaw_finalize.c
  fclaw_forestclaw.c
  fclaw_gauges.c 
  fclaw_ghost_fill.c
  fclaw_global.c
  fclaw_initialize.c
  fclaw_math.c
  fclaw_mpi.c
  fclaw_options.c 
  fclaw_output.c
  fclaw_patch.c
  fclaw_package.c 
  fclaw_physical_bc.c
  fclaw_partition.c
  fclaw_pointer_map.c
  fclaw_regrid.c
  fclaw_run.c
  fclaw_time_sync.c
  fclaw_timer.c
  fclaw_timeinterp.c
  fclaw_vtable.c

  fclaw2d_convenience.c
  fclaw2d_exchange.c
  fclaw2d_map.c
  fclaw2d_map_brick.c
  fclaw2d_map_query.c
  fclaw2d_output_tikz.c
  fclaw2d_rays.c
  fclaw2d_update_single_step.c

  forestclaw.c
  forestclaw2d.c
  forestclaw3d.c

  fclaw3d_convenience.c
  fclaw3d_exchange.c
  fclaw3d_map.c

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

  fclaw_advance.h 
  fclaw_base.h 
  fclaw_block.h
  fclaw_convenience.h
  fclaw_corner_neighbors.h
  fclaw_diagnostics.h
  fclaw_domain.h
  fclaw_edge_neighbors.h
  fclaw_elliptic_solver.h 
  fclaw_exchange.h
  fclaw_face_neighbors.h
  fclaw_farraybox.hpp 
  fclaw_filesystem.h 
  fclaw_forestclaw.h
  fclaw_gauges.h 
  fclaw_ghost_fill.h
  fclaw_global.h 
  fclaw_include_all.h 
  fclaw_options.h 
  fclaw_output.h
  fclaw_package.h 
  fclaw_partition.h
  fclaw_patch.h
  fclaw_physical_bc.h
  fclaw_pointer_map.h 
  fclaw_regrid.h
  fclaw_timer.h 
  fclaw_time_sync.h
  fclaw_timeinterp.h
  fclaw_mpi.h 
  fclaw_math.h 
  fclaw_vtable.h

  fp_exception_glibc_extension.h 

  fclaw2d_convenience.h
  fclaw2d_defs.h 
  fclaw2d_domain.h
  fclaw2d_map.h
  fclaw2d_map_brick.h
  fclaw2d_map_query.h
  fclaw2d_map_query_defs.h
  fclaw2d_rays.h
  fclaw2d_to_3d.h
  fclaw2d_update_single_step.h

  fclaw3d_convenience.h
  fclaw3d_domain.h
  fclaw3d_defs.h
  fclaw3d_include_all.h 
  fclaw3d_map.h

  forestclaw.h
  forestclaw2d.h 
  forestclaw3d.h
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
      fclaw_pointer_map.h.TEST.cpp
      fclaw_elliptic_solver.h.TEST.cpp
      fclaw_diagnostics.h.TEST.cpp
      fclaw_global.h.TEST.cpp
      fclaw_options.h.TEST.cpp
      fclaw_patch.h.TEST.cpp
      fclaw_vtable.h.TEST.cpp
  )

  target_link_libraries(forestclaw.TEST testutils forestclaw)
  register_unit_tests(forestclaw.TEST)
endif(BUILD_TESTING)
