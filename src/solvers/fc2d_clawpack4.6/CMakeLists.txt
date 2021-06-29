# -- fortran library
add_library(clawpack4.6_f OBJECT
  fortran_source/clawpack46_time_sync.f 
  fortran_source/clawpack46_inlinelimiter.f 
  fortran_source/clawpack46_set_capacity.f
  fortran_source/clawpack46_block.f 
  fortran_source/clawpack46_bc2_default.f 
  fortran_source/clawpack46_flux2.f
  fortran_source/clawpack46_flux2fw.f 
  fortran_source/clawpack46_step2.f 
  fortran_source/clawpack46_step2_wrap.f
)

target_include_directories(clawpack4.6_f INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)

# -- add library
add_library(clawpack4.6)
set_target_properties(clawpack4.6 PROPERTIES EXPORT_NAME CLAWPACK4.6)

target_sources(clawpack4.6 PRIVATE 
  fc2d_clawpack46.cpp 
  fc2d_clawpack46_options.c
  $<TARGET_OBJECTS:clawpack4.6_f>
)

target_link_libraries(clawpack4.6 PUBLIC clawpatch)

target_include_directories(clawpack4.6
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

# -- install

install(FILES
	fc2d_clawpack46.h
	fc2d_clawpack46_options.h
	clawpack46_user_fort.h
	fc2d_clawpack46_fort.h 
  DESTINATION include
)

install(TARGETS clawpack4.6
  EXPORT ${PROJECT_NAME}Targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)


# -- imported target, for use from FetchContent
add_library(FORESTCLAW::CLAWPACK4.6 INTERFACE IMPORTED GLOBAL)
target_link_libraries(FORESTCLAW::CLAWPACK4.6 INTERFACE clawpack4.6)