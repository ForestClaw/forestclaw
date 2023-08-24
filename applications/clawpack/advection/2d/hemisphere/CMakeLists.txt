set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(hemisphere_f OBJECT 
  setprob.f 
  fdisc.f 
  psi.f 
  user_4.6/qinit.f 
  user_5.0/qinit.f90 
  ${all}/user46_setaux_manifold.f 
  ${all}/user5_setaux_manifold.f90 
  ${all}/clawpack46_rpn2adv_manifold.f 
  ${all}/clawpack46_rpt2adv_manifold.f 
  ${all}/clawpack5_rpn2adv_manifold.f90 
  ${all}/clawpack5_rpt2adv_manifold.f90 
)

add_executable(hemisphere 
  hemisphere.cpp 
  hemisphere_options.c 
  hemisphere_user.cpp 
  ${all}/advection_patch_setup_manifold.c
  ${all}/fclaw2d_map_pillowsphere.c 
  ${all}/fclaw2d_map_pillowsphere5.c 
  $<TARGET_OBJECTS:hemisphere_f>
)

target_include_directories(hemisphere PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(hemisphere PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_test(NAME clawpack_advection_2d_hemisphere COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_advection_2d_hemisphere PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")