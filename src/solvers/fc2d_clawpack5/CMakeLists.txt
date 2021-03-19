add_library(clawpack5_f OBJECT
  fortran_source/clawpack5_amr_module.f90 
  fortran_source/clawpack5_set_amr_module.f90 
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

target_include_directories(clawpack5_f INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)

target_sources(clawpack5 PRIVATE 
  fc2d_clawpack5.cpp	
  fc2d_clawpack5_options.c
  $<TARGET_OBJECTS:clawpack5_f>
)

target_include_directories(clawpack5 PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_BINARY_DIR}/include
)
