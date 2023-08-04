set(metric ${PROJECT_SOURCE_DIR}/src/patches/metric)

# -- fortran object library

add_library(clawpatch_f OBJECT
  fclaw2d_clawpatch_utils.f
  fclaw3dx_clawpatch_utils.f
  fort_user/fclaw2d_clawpatch_difference_exceeds_th.f90
  fort_user/fclaw2d_clawpatch_value_exceeds_th.f90
  fort_user/fclaw2d_clawpatch_minmax_exceeds_th.f90
  fort_user/fclaw2d_clawpatch_gradient_exceeds_th.f90
  fort3_user/fclaw3dx_clawpatch_difference_exceeds_th.f90
  fort3_user/fclaw3dx_clawpatch_value_exceeds_th.f90
  fort3_user/fclaw3dx_clawpatch_minmax_exceeds_th.f90
  fort3_user/fclaw3dx_clawpatch_gradient_exceeds_th.f90
  fort_4.6/fclaw2d_clawpatch46_average.f
  fort_4.6/fclaw2d_clawpatch46_copy.f
  fort_4.6/fclaw2d_clawpatch46_interpolate.f
  fort_4.6/fclaw2d_clawpatch46_local_ghost_pack.f
  fort_4.6/fclaw2d_clawpatch46_diagnostics.f
  fort_4.6/fclaw2d_clawpatch46_pillowsphere.f
  fort_4.6/fclaw2d_clawpatch46_output_ascii.f
  fort_4.6/fclaw2d_clawpatch46_tag4refinement.f
  fort_4.6/fclaw2d_clawpatch46_tag4coarsening.f
  fort_4.6/fclaw2d_clawpatch46_timeinterp.f
  fort_5.0/fclaw2d_clawpatch5_average.f
  fort_5.0/fclaw2d_clawpatch5_copy.f
  fort_5.0/fclaw2d_clawpatch5_interpolate.f
  fort_5.0/fclaw2d_clawpatch5_local_ghost_pack.f
  fort_5.0/fclaw2d_clawpatch5_diagnostics.f
  fort_5.0/fclaw2d_clawpatch5_pillowsphere.f
  fort_5.0/fclaw2d_clawpatch5_output_ascii.f
  fort_5.0/fclaw2d_clawpatch5_tag4refinement.f
  fort_5.0/fclaw2d_clawpatch5_tag4coarsening.f
  fort_5.0/fclaw2d_clawpatch5_timeinterp.f
  fort3_4.6/fclaw3dx_clawpatch46_tag4refinement.f90 
  fort3_4.6/fclaw3dx_clawpatch46_tag4coarsening.f90 
  fort3_4.6/fclaw3dx_clawpatch46_average.f90 
  fort3_4.6/fclaw3dx_clawpatch46_copy.f90 
  fort3_4.6/fclaw3dx_clawpatch46_interpolate.f90 
  fort3_4.6/fclaw3dx_clawpatch46_local_ghost_pack.f90 
  fort3_4.6/fclaw3dx_clawpatch46_output_ascii.f90 
  fort3_4.6/fclaw3dx_clawpatch46_timeinterp.f90  
  fort3_4.6/fclaw3dx_clawpatch46_pillowsphere.f90 
  fort3_4.6/fclaw3dx_clawpatch46_diagnostics.f90 
  ${metric}/fortran_source2d/fclaw2d_metric_default_fort.f90
  ${metric}/fortran_source3d/fclaw3d_metric_default_fort.f90
)

target_include_directories(clawpatch_f PRIVATE
  ${metric}/fortran_source2d
)

set_target_properties(clawpatch_f PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

target_include_directories(clawpatch_f INTERFACE
  ${CMAKE_CURRENT_BINARY_DIR}/include
)

# -- clawpatch library

add_library(clawpatch
  $<TARGET_OBJECTS:clawpatch_f>
  fclaw_clawpatch_options.c
  fclaw_clawpatch_diagnostics.c
  fclaw_clawpatch_diagnostics_default.c

  fclaw2d_clawpatch.cpp
  fclaw2d_clawpatch_pillow.c
  fclaw2d_clawpatch_transform.c
  fclaw2d_clawpatch_output_ascii.c
  fclaw2d_clawpatch_output_vtk.c
  fclaw2d_clawpatch_conservation.c

  fclaw3dx_clawpatch.cpp
  fclaw3dx_clawpatch_pillow.c
  fclaw3dx_clawpatch_transform.c
  fclaw3dx_clawpatch_output_ascii.c
  fclaw3dx_clawpatch_output_vtk.c

  fort_user/fclaw2d_clawpatch_tag_criteria.c
  fort3_user/fclaw3dx_clawpatch_tag_criteria.c
  ${metric}/fclaw2d_metric.cpp
  ${metric}/fclaw2d_metric_default.c
  ${metric}/fclaw3d_metric.cpp 
  ${metric}/fclaw3d_metric_default.c
)

target_include_directories(clawpatch
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${metric}/fortran_source2d>
  $<BUILD_INTERFACE:${metric}/fortran_source3d>
  $<BUILD_INTERFACE:${metric}>
  $<INSTALL_INTERFACE:include>
)

target_link_libraries(clawpatch PUBLIC forestclaw)

# -- install

install(FILES
  fclaw_clawpatch_enums.h

	fclaw_clawpatch.hpp
	fclaw_clawpatch_diagnostics.h
	fclaw_clawpatch_options.h

  fclaw2d_clawpatch.h
	fclaw2d_clawpatch_conservation.h
	fclaw2d_clawpatch_transform.h
	fclaw2d_clawpatch_pillow.h
	fclaw2d_clawpatch_fort.h
	fclaw2d_clawpatch46_fort.h
	fclaw2d_clawpatch5_fort.h
	fclaw2d_clawpatch_output_ascii.h
	fclaw2d_clawpatch_output_vtk.h

  fclaw3dx_clawpatch.h
	fclaw3dx_clawpatch_transform.h
	fclaw3dx_clawpatch_diagnostics.h
	fclaw3dx_clawpatch_pillow.h
	fclaw3dx_clawpatch_fort.h
	fclaw3dx_clawpatch46_fort.h
	fclaw3dx_clawpatch_output_ascii.h
	fclaw3dx_clawpatch_output_vtk.h

	${metric}/fclaw2d_metric.h
	${metric}/fclaw2d_metric.hpp
	${metric}/fclaw2d_metric_default_fort.h
	${metric}/fortran_source2d/fclaw2d_metric_terms.i 

  ${metric}/fclaw3d_metric.h
  ${metric}/fclaw3d_metric.hpp
  ${metric}/fclaw3d_metric_default_fort.h
  DESTINATION include
)

install(TARGETS clawpatch
  EXPORT ${PROJECT_NAME}-targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)

# -- export name
set_target_properties(clawpatch PROPERTIES EXPORT_NAME CLAWPATCH)

# -- imported target, for use from FetchContent
add_library(FORESTCLAW::CLAWPATCH INTERFACE IMPORTED GLOBAL)
target_link_libraries(FORESTCLAW::CLAWPATCH INTERFACE clawpatch)

if(BUILD_TESTING)
  add_executable(clawpatch.TEST
      fclaw_clawpatch_options.h.TEST.cpp
      fclaw_clawpatch_diagnostics.h.TEST.cpp
      fclaw2d_clawpatch.h.TEST.cpp
      fclaw2d_clawpatch_fort.h.TEST.cpp
      fclaw3dx_clawpatch.h.TEST.cpp
      ${metric}/fclaw2d_metric.h.TEST.cpp
  )
  target_link_libraries(clawpatch.TEST testutils clawpatch forestclaw)
  register_unit_tests(clawpatch.TEST)
endif(BUILD_TESTING)
