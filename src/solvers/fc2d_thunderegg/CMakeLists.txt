
add_library(thunderegg_f OBJECT
	fort_4.6/thunderegg_apply_bc_default.f90 
	fort_4.6/thunderegg_eval_bc_default.f90
)

set_target_properties(thunderegg_f PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

target_include_directories(thunderegg_f INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)

target_sources(thunderegg PRIVATE 
  fc2d_thunderegg.cpp
  fc2d_thunderegg_options.c 
  fc2d_thunderegg_vector.cpp 
  fc2d_thunderegg_physical_bc.c 
  operators/fc2d_thunderegg_starpatch.cpp 
  operators/fc2d_thunderegg_fivepoint.cpp 
  operators/fc2d_thunderegg_varpoisson.cpp 
  operators/fc2d_thunderegg_heat.cpp 
  $<TARGET_OBJECTS:thunderegg_f>
)

#TODO clean this up
target_include_directories(thunderegg PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/operators
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_BINARY_DIR}/include
)
