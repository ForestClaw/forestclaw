# -- fortran library
add_library(geoclaw_f OBJECT
    amrlib_source/amr_module.f90
    geolib_source/utility_module.f90
    geolib_source/geoclaw_module.f90
    amrlib_source/regions_module.f90
    geolib_source/topo_module.f90
    geolib_source/qinit_module.f90
    geolib_source/fixedgrids_module.f90
    geolib_source/refinement_module.f90
    geolib_source/surge/model_storm_module.f90
    geolib_source/surge/data_storm_module.f90
    geolib_source/surge/storm_module.f90
    geolib_source/friction_module.f90
    riemann_source/rpn2_geoclaw.f
    riemann_source/rpt2_geoclaw.f
    riemann_source/geoclaw_riemann_utils.f
    amrlib_source/opendatafile.f
    amrlib_source/check4nans.f90
    amrlib_source/quick_sort1.f
    geolib_source/topointegral.f
    geolib_source/topo_update.f90
    geolib_source/cellgridintegrate2.f
    geolib_source/allowflag.f
    geolib_source/bilinearintegral.f
    fortran_source/philim.f
    fortran_source/geoclaw_set_modules_fort.f90
    fortran_source/geoclaw_step2_wrap_fort.f
    fortran_source/geoclaw_flux2fw_fort.f
    fortran_source/geoclaw_limiter_fort.f
    fortran_source/geoclaw_block_fort.f
    fortran_source/geoclaw_gauges_fort.f90
    fortran_source/geoclaw_topo_update_fort.f
    fortran_source/geoclaw_bc2_fort.f
    fortran_source/geoclaw_step2_fort.f90
    fortran_source/geoclaw_setaux_fort.f90
    fortran_source/geoclaw_src2_fort.f90
    fortran_source/geoclaw_b4step2_fort.f90
    fortran_source/geoclaw_qinit_fort.f90
    fclaw2d_source/fc2d_geoclaw_copy_fort.f
    fclaw2d_source/fc2d_geoclaw_average_fort.f
    fclaw2d_source/fc2d_geoclaw_interpolate_fort.f
    fclaw2d_source/fc2d_geoclaw_output_ascii_fort.f
    fclaw2d_source/fc2d_geoclaw_tag4refinement_fort.f90
    fclaw2d_source/fc2d_geoclaw_tag4coarsening_fort.f90
    fclaw2d_source/fc2d_geoclaw_test_regions.f90
    fclaw2d_source/fc2d_geoclaw_flag2refine_fort.f90
    fclaw2d_source/fc2d_geoclaw_local_ghost_pack_fort.f
    fclaw2d_source/fc2d_geoclaw_local_ghost_pack_aux_fort.f
    fclaw2d_source/fc2d_geoclaw_diagnostics_fort.f
    fclaw2d_source/fc2d_geoclaw_timeinterp_fort.f
)

target_link_Libraries(geoclaw_f PRIVATE
    clawpatch
)

set_target_properties(geoclaw_f PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

# -- add library
add_library(geoclaw)
set_target_properties(geoclaw PROPERTIES EXPORT_NAME GEOCLAW)

target_sources(geoclaw PRIVATE
    $<TARGET_OBJECTS:geoclaw_f>
    fc2d_geoclaw.cpp
    fc2d_geoclaw_options.c
    fc2d_geoclaw_gauges_default.c
    fc2d_geoclaw_run.c
    fc2d_geoclaw_output_ascii.c
)

target_link_libraries(geoclaw PUBLIC forestclaw clawpatch)

target_include_directories(geoclaw PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
)

# -- install
install(FILES
  fc2d_geoclaw.h
  types.h
  fc2d_geoclaw_options.h
  fc2d_geoclaw_fort.h
  fc2d_geoclaw_gauges_default.h
  DESTINATION include
)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/include/amr_module.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include/data_storm_module.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include/fixedgrids_module.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include/friction_module.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include/geoclaw_module.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include/model_storm_module.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include/qinit_module.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include/refinement_module.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include/regions_module.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include/storm_module.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include/topo_module.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include/utility_module.mod
  DESTINATION include
)

install(TARGETS geoclaw
  EXPORT ${PROJECT_NAME}-targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)


# imported target, for use from FetchContent
add_library(FORESTCLAW::GEOCLAW INTERFACE IMPORTED GLOBAL)
target_link_libraries(FORESTCLAW::GEOCLAW INTERFACE geoclaw)

if(BUILD_TESTING)
  add_executable(fc2d_geoclaw.TEST
    fc2d_geoclaw.h.TEST.cpp
    fc2d_geoclaw_options.h.TEST.cpp
  )
  target_link_libraries(fc2d_geoclaw.TEST testutils geoclaw)
  register_unit_tests(fc2d_geoclaw.TEST)
endif(BUILD_TESTING)
