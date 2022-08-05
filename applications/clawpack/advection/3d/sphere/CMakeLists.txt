set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/3d/all)
set(all2d ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(sphere3dx_f OBJECT 
  psi.f 
  fdisc.f 
  setprob.f 
  user_4.6/qinit.f90 
  ${all}/claw3_setaux_manifold.f90
  ${all}/rpn3_vc_advection.f90 
  ${all}/rpt3_vc_advection.f90 
  ${all}/rptt3_vc_advection.f90 
)

add_executable(sphere3dx
  sphere.cpp 
  sphere_user.cpp 
  sphere_options.c 
  ${all}/advection_patch_setup_manifold.c
  ${all2d}/fclaw2d_map_cubedsphere.c 
  ${all2d}/fclaw2d_map_pillowsphere.c 
  $<TARGET_OBJECTS:sphere3dx_f>
)

target_include_directories(sphere3dx PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(sphere3dx PRIVATE 
 FORESTCLAW::CLAWPACK3_46
)

add_test(NAME clawpack_advection_3d_sphere COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_advection_3d_sphere PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")