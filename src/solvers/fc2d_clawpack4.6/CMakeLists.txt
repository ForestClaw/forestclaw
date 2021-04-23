set(f fortran_source/)
add_library(clawpack4.6_f OBJECT
  ${f}clawpack46_time_sync.f ${f}clawpack46_inlinelimiter.f ${f}clawpack46_set_capacity.f
  ${f}clawpack46_block.f ${f}clawpack46_bc2_default.f ${f}clawpack46_flux2.f
  ${f}clawpack46_flux2fw.f ${f}clawpack46_step2.f ${f}clawpack46_step2_wrap.f)
set_target_properties(clawpack4.6_f PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
target_include_directories(clawpack4.6_f INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)

target_sources(clawpack4.6 PRIVATE fc2d_clawpack46.cpp fc2d_clawpack46_options.c
  $<TARGET_OBJECTS:clawpack4.6_f>)
target_include_directories(clawpack4.6 PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_BINARY_DIR}/include)
