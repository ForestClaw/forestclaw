# -- fortran library
add_library(clawpack3_46_f OBJECT
  fortran_source/fc3d_clawpack46_bc3_default.f90 
  fortran_source/fc3d_clawpack46_inlinelimiter.f90 
  fortran_source/fc3d_clawpack46_limiter.f90 
  fortran_source/fc3d_clawpack46_set_capacity.f90
  fortran_source/fc3d_clawpack46_block.f90 
  fortran_source/fc3d_clawpack46_bc3_default.f90 
  fortran_source/fc3d_clawpack46_flux3.f
  fortran_source/fc3d_clawpack46_step3.f90 
  fortran_source/fc3d_clawpack46_step3_wrap.f90
)

target_include_directories(clawpack3_46_f INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)

# -- add library
add_library(clawpack3_46)
set_target_properties(clawpack3_46 PROPERTIES EXPORT_NAME CLAWPACK3_46)

target_sources(clawpack3_46 PRIVATE 
  fc3d_clawpack46.cpp 
  fc3d_clawpack46_options.c
  $<TARGET_OBJECTS:clawpack3_46_f>
)

target_link_libraries(clawpack3_46 PUBLIC clawpatch)

target_include_directories(clawpack3_46
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

# -- install

install(FILES
	fc3d_clawpack46.h
  fc3d_clawpack46_fort.h 
	fc3d_clawpack46_options.h
  fc3d_clawpack46_user_fort.h 
DESTINATION include
)

install(TARGETS clawpack3_46
  EXPORT ${PROJECT_NAME}Targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)


# -- imported target, for use from FetchContent
add_library(FORESTCLAW::CLAWPACK3_46 INTERFACE IMPORTED GLOBAL)
target_link_libraries(FORESTCLAW::CLAWPACK3_46 INTERFACE clawpack3_46)

if(BUILD_TESTING)
  add_executable(fc3d_clawpack46.TEST
    fc3d_clawpack46.h.TEST.cpp
    fc3d_clawpack46_options.h.TEST.cpp
  )
  target_link_libraries(fc3d_clawpack46.TEST testutils clawpack3_46)
  register_unit_tests(fc3d_clawpack46.TEST)
endif(BUILD_TESTING)