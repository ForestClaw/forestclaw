# -- fortran library
add_library(clawpack5_f OBJECT
  fortran_source/clawpack5_amr_module.f90
  fortran_source/clawpack5_set_amr_module.f90
  fortran_source/clawpack5_time_sync.f
  fortran_source/clawpack5_inlinelimiter.f
  fortran_source/clawpack5_limiter.f
  fortran_source/clawpack5_set_capacity.f
  fortran_source/clawpack5_block.f
  fortran_source/clawpack5_bc2_default.f
  fortran_source/clawpack5_flux2.f
  fortran_source/clawpack5_step2.f90
  fortran_source/clawpack5_step2_wrap.f
)

set_target_properties(clawpack5_f PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

# -- add library
add_library(clawpack5)
set_target_properties(clawpack5 PROPERTIES EXPORT_NAME CLAWPACK5)

target_sources(clawpack5 PRIVATE
  fc2d_clawpack5.cpp
  fc2d_clawpack5_options.c
  $<TARGET_OBJECTS:clawpack5_f>
)

target_link_libraries(clawpack5 PUBLIC clawpatch)
target_include_directories(clawpack5
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

# -- install
install(FILES
	fc2d_clawpack5.h
	fc2d_clawpack5_options.h
	clawpack5_user_fort.h
	fc2d_clawpack5_fort.h
  TYPE INCLUDE
)

install(TARGETS clawpack5 EXPORT ${PROJECT_NAME}-targets)

# imported target, for use from parent project
add_library(FORESTCLAW::CLAWPACK5 INTERFACE IMPORTED GLOBAL)
target_link_libraries(FORESTCLAW::CLAWPACK5 INTERFACE clawpack5)

if(BUILD_TESTING)
  add_executable(fc2d_clawpack5.TEST
    fc2d_clawpack5.h.TEST.cpp
    fc2d_clawpack5_options.h.TEST.cpp
  )
  target_link_libraries(fc2d_clawpack5.TEST testutils clawpack5)
  register_unit_tests(fc2d_clawpack5.TEST)
endif(BUILD_TESTING)
