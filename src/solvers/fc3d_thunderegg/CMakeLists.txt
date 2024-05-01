
# -- fortran library
add_library(fc3d_thunderegg_f OBJECT
	fort_4.6/thunderegg_apply_bc_default.f90
	fort_4.6/thunderegg_eval_bc_default.f90
)

set_target_properties(fc3d_thunderegg_f PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

target_include_directories(fc3d_thunderegg_f INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)

# -- add library
add_library(fc3d_thunderegg)
set_target_properties(fc3d_thunderegg PROPERTIES EXPORT_NAME FC3D_THUNDEREGG)

target_sources(fc3d_thunderegg PRIVATE
  fc3d_thunderegg.cpp
  fc3d_thunderegg_options.c
  fc3d_thunderegg_vector.cpp
  fc3d_thunderegg_physical_bc.c
  $<TARGET_OBJECTS:fc3d_thunderegg_f>
)

target_link_libraries(fc3d_thunderegg PUBLIC forestclaw clawpatch ThunderEgg::ThunderEgg)

target_include_directories(fc3d_thunderegg PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

# -- install
install(FILES
	fc3d_thunderegg.h
	fc3d_thunderegg_fort.h
	fc3d_thunderegg_options.h
	fc3d_thunderegg_physical_bc.h
	fc3d_thunderegg_vector.hpp
  TYPE INCLUDE
)

install(TARGETS fc3d_thunderegg
  EXPORT ${PROJECT_NAME}-targets)


# imported target, for use from FetchContent
add_library(FORESTCLAW::FC3D_THUNDEREGG INTERFACE IMPORTED GLOBAL)
target_link_libraries(FORESTCLAW::FC3D_THUNDEREGG INTERFACE fc3d_thunderegg)

if(BUILD_TESTING)
  add_executable(fc3d_thunderegg.TEST
    fc3d_thunderegg.h.TEST.cpp
    fc3d_thunderegg_options.h.TEST.cpp
    fc3d_thunderegg_vector_TEST.cpp
  )
  target_link_libraries(fc3d_thunderegg.TEST testutils fc3d_thunderegg forestclaw)
  register_unit_tests(fc3d_thunderegg.TEST)
endif(BUILD_TESTING)
