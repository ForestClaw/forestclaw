
# -- fortran library
add_library(fc2d_thunderegg_f OBJECT
	fort_4.6/thunderegg_apply_bc_default.f90
	fort_4.6/thunderegg_eval_bc_default.f90
)

set_target_properties(fc2d_thunderegg_f PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

target_include_directories(fc2d_thunderegg_f INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)

# -- add library
add_library(fc2d_thunderegg)
set_target_properties(fc2d_thunderegg PROPERTIES EXPORT_NAME FC2D_THUNDEREGG)

target_sources(fc2d_thunderegg PRIVATE
  fc2d_thunderegg.cpp
  fc2d_thunderegg_options.c
  fc2d_thunderegg_vector.cpp
  fc2d_thunderegg_physical_bc.c
  operators/fc2d_thunderegg_starpatch.cpp
  operators/fc2d_thunderegg_fivepoint.cpp
  operators/fc2d_thunderegg_varpoisson.cpp
  operators/fc2d_thunderegg_heat.cpp
  $<TARGET_OBJECTS:fc2d_thunderegg_f>
)

target_link_libraries(fc2d_thunderegg PUBLIC forestclaw clawpatch ThunderEgg::ThunderEgg)

target_include_directories(fc2d_thunderegg PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/operators>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

# -- install
install(FILES
	fc2d_thunderegg.h
	fc2d_thunderegg_fort.h
	fc2d_thunderegg_options.h
	fc2d_thunderegg_physical_bc.h
	fc2d_thunderegg_vector.hpp
	operators/fc2d_thunderegg_starpatch.h
	operators/fc2d_thunderegg_fivepoint.h
	operators/fc2d_thunderegg_varpoisson.h
	operators/fc2d_thunderegg_heat.h
  DESTINATION include
)

install(TARGETS fc2d_thunderegg
  EXPORT ${PROJECT_NAME}-targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)


# imported target, for use from FetchContent
add_library(FORESTCLAW::FC2D_THUNDEREGG INTERFACE IMPORTED GLOBAL)
target_link_libraries(FORESTCLAW::FC2D_THUNDEREGG INTERFACE fc2d_thunderegg)

if(BUILD_TESTING)
  add_executable(fc2d_thunderegg.TEST
    fc2d_thunderegg.h.TEST.cpp
    fc2d_thunderegg_vector_TEST.cpp
  )
  target_link_libraries(fc2d_thunderegg.TEST testutils fc2d_thunderegg forestclaw)
  register_unit_tests(fc2d_thunderegg.TEST)
endif(BUILD_TESTING)
