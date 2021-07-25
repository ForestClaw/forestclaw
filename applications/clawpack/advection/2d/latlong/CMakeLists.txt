set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(latlong_f OBJECT 
  fdisc.f 
  setprob.f 
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

add_executable(latlong 
  latlong.cpp 
  latlong_options.c 
  latlong_user.cpp 
  ${all}/advection_patch_setup_manifold.c
  ${all}/fclaw2d_map_latlong.c 
  $<TARGET_OBJECTS:latlong_f>
)

target_include_directories(latlong PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(latlong PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_test(NAME clawpack_advection_2d_latlong COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_advection_2d_latlong PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")