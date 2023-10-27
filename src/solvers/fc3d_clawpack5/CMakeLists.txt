# -- fortran library
add_library(clawpack5_3d_f OBJECT
	fclaw2d_source/fc3d_clawpack5_average.f
	fclaw2d_source/fc3d_clawpack5_copy.f
	fclaw2d_source/fc3d_clawpack5_interpolate.f
	fclaw2d_source/fc3d_clawpack5_ghostpack_qarea.f
	fclaw2d_source/fc3d_clawpack5_diagnostics.f
	fclaw2d_source/fc3d_clawpack5_tag4refinement.f
	fclaw2d_source/fc3d_clawpack5_tag4coarsening.f
	fclaw2d_source/fc3d_clawpack5_output_ascii.f
	fclaw2d_source/fc3d_clawpack5_timeinterp.f
)

set_target_properties(clawpack5_3d_f PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

# -- add library
add_library(clawpack5_3d STATIC)
set_target_properties(clawpack5_3d PROPERTIES EXPORT_NAME CLAWPACK5_3D)

target_sources(clawpack5_3d PRIVATE
    fc3d_clawpack5.cpp
    fc3d_clawpack5_options.c
	$<TARGET_OBJECTS:clawpack5_3d_f>
)

target_include_directories(clawpack5_3d
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

# -- install
install(FILES
	fc3d_clawpack5.h
	clawpack5_user_fort.h
	fc3d_clawpack5_options.h
  TYPE INCLUDE
)

install(TARGETS clawpack5_3d EXPORT ${PROJECT_NAME}-targets)

# imported target, for use from FetchContent
add_library(FORESTCLAW::CLAWPACK5_3D INTERFACE IMPORTED GLOBAL)
target_link_libraries(FORESTCLAW::CLAWPACK5_3D INTERFACE clawpack5_3d)
