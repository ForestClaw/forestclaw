set(all ${PROJECT_SOURCE_DIR}/clawpack/transport/2d/all)

add_library(gaussian_f OBJECT 
  gaussian_sum.f 
  user_4.6/qinit.f 
  user_5.0/qinit.f90   
  ${all}/transport_util.f 
  ${all}/transport_setprob.f 
  ${all}/velocity_wind.f 
  ${all}/user46_b4step2_manifold.f 
  ${all}/user5_b4step2_manifold.f90 
  ${all}/user46_setaux_manifold.f 
  ${all}/user5_setaux_manifold.f90 
  ${all}/clawpack46_rpn2adv_manifold.f 
  ${all}/clawpack46_rpt2adv_manifold.f 
  ${all}/clawpack5_rpn2adv_manifold.f90 
  ${all}/clawpack5_rpt2adv_manifold.f90
)

add_executable(gaussian 
  gaussian.cpp 
  gaussian_user.cpp 
  gaussian_options.c 
  ${all}/transport_user.cpp
  ${all}/fclaw2d_map_cubedsphere.c 
  ${all}/fclaw2d_map_pillowsphere.c 
  $<TARGET_OBJECTS:gaussian_f>
)

target_include_directories(gaussian PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(gaussian PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_test(NAME clawpack_transport_gaussian COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_transport_gaussian PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")

